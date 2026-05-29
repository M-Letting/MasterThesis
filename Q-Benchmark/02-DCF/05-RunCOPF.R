options(warn = -1)

suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(arrow, warn.conflicts = FALSE))
suppressMessages(library(VIM, warn.conflicts = FALSE))
# Requires the proteoformLocationMapping branch:
# remotes::install_github("CCprofiler/CCprofiler", ref = "proteoformLocationMapping")
suppressMessages(library(CCprofiler, warn.conflicts = FALSE))

setwd(here::here("Q-Benchmark"))

score_cutoff <- 0.1
adj_pval_cutoff <- 0.05

# Run toggles
run_proteomaker <- TRUE
run_proteomaker_na <- FALSE
run_proteomaker_na_imputed <- TRUE
run_complex_raw <- FALSE
run_complex_imputed <- FALSE

output_dir <- "./02-DCF/00-Data/03-Results/COPF"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

################################################################################
## Helpers                                                                     ##
################################################################################

format_time_diff <- function(td) {
  td <- as.numeric(td, units = "secs")
  if (td < 60) {
    return(paste(round(td, 2), "seconds"))
  }
  if (td < 3600) {
    return(paste(round(td / 60, 2), "minutes"))
  }
  paste(round(td / 3600, 2), "hours")
}

# Run the full CCprofiler COPF pipeline on a long-format data.table.
# long_data    : data.table with 'peptide_id', 'protein_id', 'filename', 'intensity' columns
# fraction_annot: data.table with 'filename' and 'fraction_number'
# trace_annot  : data.table with 'peptide_id' and 'protein_id' columns (for annotateTraces)
# label        : string used in messages and PDF name (if ever enabled)
# Returns the trace_annotation data.table, or NULL on failure.
run_copf_pipeline <- function(long_data, fraction_annot, trace_annot, label) {
  # Drop NA-intensity rows so importPCPdata's dcast(fill=0) treats them as
  # missing combinations and fills with 0, matching the reference feather format.
  long_data <- long_data[!is.na(intensity)]

  # Collapse duplicate (peptide_id, protein_id, filename) rows by averaging
  # intensities. Duplicates occur when the same peptide sequence appears in
  # multiple rows of the source data (e.g., same peptide mapped to the complex
  # more than once). Without deduplication, importPCPdata's dcast falls back to
  # counting occurrences (all 1s), wiping out all intensity variance.
  n_before_dedup <- nrow(long_data)
  long_data <- long_data[,
    .(intensity = mean(intensity)),
    by = .(peptide_id, protein_id, filename)
  ]
  n_removed <- n_before_dedup - nrow(long_data)
  if (n_removed > 0) {
    cat(
      "  Collapsed",
      n_removed,
      "duplicate (peptide, fraction) rows by averaging\n"
    )
  }

  traces <- tryCatch(
    CCprofiler::importPCPdata(
      input_data = long_data,
      fraction_annotation = fraction_annot
    ),
    error = function(e) {
      message("importPCPdata failed for ", label, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(traces)) {
    return(NULL)
  }

  # Diagnostic: check ID overlap between traces and trace_annot
  trace_ids <- traces$trace_annotation$id
  annot_ids <- trace_annot$peptide_id
  n_match <- sum(trace_ids %in% annot_ids)
  cat(
    "  [diag] traces ids (n=",
    length(trace_ids),
    ") sample:",
    paste(head(trace_ids, 2), collapse = " | "),
    "\n"
  )
  cat(
    "  [diag] trace_annot ids (n=",
    length(annot_ids),
    ") sample:",
    paste(head(annot_ids, 2), collapse = " | "),
    "\n"
  )
  cat("  [diag] matching ids:", n_match, "\n")

  traces <- tryCatch(
    CCprofiler::annotateTraces(
      traces,
      trace_annot,
      traces_id_column = "id",
      trace_annotation_id_column = "peptide_id",
      uniprot_mass_format = FALSE
    ),
    error = function(e) {
      message("annotateTraces failed for ", label, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(traces)) {
    return(NULL)
  }

  int_mat <- CCprofiler::getIntensityMatrix(traces)
  zerovar <- apply(int_mat, 1, var, na.rm = TRUE) > 0
  traces <- subset(traces, trace_subset_ids = names(zerovar[zerovar]))

  n_before <- nrow(traces$trace_annotation)
  traces <- CCprofiler::filterSinglePeptideHits(traces)
  n_after <- nrow(traces$trace_annotation)

  cat(
    "  Peptides remaining after zero-variance + single-peptide filter:",
    n_after,
    "(removed",
    n_before - n_after,
    ")\n"
  )

  if (n_after == 0) {
    message("No multi-peptide proteins remaining for ", label)
    return(NULL)
  }

  traces_corr <- tryCatch(
    CCprofiler::calculateGeneCorrMatrices(traces),
    error = function(e) {
      message(
        "calculateGeneCorrMatrices failed for ",
        label,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  if (is.null(traces_corr)) {
    return(NULL)
  }

  traces_clustered <- CCprofiler::clusterPeptides(
    traces_corr,
    method = "average",
    plot = FALSE,
    PDF = FALSE,
    name = paste0("COPF_", label)
  )

  traces_cut <- CCprofiler::cutClustersInNreal(
    traces_clustered,
    clusterN = 2,
    min_peptides_per_cluster = 2
  )

  traces_scored <- tryCatch(
    CCprofiler::calculateProteoformScore(traces_cut),
    error = function(e) {
      message("calculateProteoformScore failed for ", label, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(traces_scored)) {
    return(NULL)
  }

  # Guard: calculateProteoformScore does not add proteoform_score_pval_adj when
  # all proteins remain in a single cluster (the merge never runs).
  if (!"proteoform_score_pval_adj" %in% names(traces_scored$trace_annotation)) {
    traces_scored$trace_annotation[, proteoform_score_pval_adj := NA_real_]
  }

  # Use fully permissive thresholds so every protein with a two-cluster solution
  # contributes to the output. The actual proteoform_score_pval_adj flows into
  # peptide_adj_pval for minority-cluster peptides, making the external benchmark
  # threshold sweep the real sensitivity/specificity control knob.
  # The module-level score_cutoff / adj_pval_cutoff are intentionally not used here.
  traces_pf <- CCprofiler::annotateTracesWithProteoforms(
    traces_scored,
    score_cutoff = 0,
    adj_pval_cutoff = 1
  )

  res <- as.data.table(traces_pf$trace_annotation)
  setnames(res, "id", "peptide_id")

  # peptide_adj_pval: for proteins with n_proteoforms > 1, the minority cluster
  # (fewer members, excluding singletons) is the discordant complexoform/proteoform;
  # assign proteoform_score_pval_adj to those peptides, 1.0 to all others.
  res[, n_in_cluster := .N, by = .(protein_id, cluster)]
  res[, peptide_adj_pval := 1.0]
  res[
    n_proteoforms > 1 & cluster != 100,
    peptide_adj_pval := ifelse(
      n_in_cluster == min(n_in_cluster),
      proteoform_score_pval_adj,
      1.0
    ),
    by = protein_id
  ]
  res[, n_in_cluster := NULL]

  res
}

################################################################################
## ProteoMaker datasets                                                        ##
################################################################################

pm_datasets <- list(
  "LowNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerLowNA_wide_processed.csv",
  "LowNoNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerLowNoNA_wide_processed.csv",
  "MedNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerMedNA_wide_processed.csv",
  "MedNoNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerMedNoNA_wide_processed.csv",
  "HighNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerHighNA_wide_processed.csv",
  "HighNoNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerHighNoNA_wide_processed.csv"
)

if (!run_proteomaker_na) {
  pm_datasets <- pm_datasets[grepl("NoNA$", names(pm_datasets))]
}

if (run_proteomaker) {
  for (name in names(pm_datasets)) {
    file_path <- pm_datasets[[name]]
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }

    out_path <- file.path(
      output_dir,
      paste0("COPF_ProteoMaker", name, "_results.csv")
    )
    if (file.exists(out_path)) {
      cat("Skipping", name, "- result already exists\n")
      next
    }

    cat(paste0(rep("#", 80), collapse = ""), "\n")
    cat("Processing ProteoMaker dataset (COPF):", name, "\n")
    t0 <- Sys.time()

    data <- fread(file_path)
    intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

    # CCprofiler expects linear-scale intensities; prepared data are log2-transformed.
    # Build wide with peptide_id + protein_id for back-transform, then melt to long format.
    copf_wide <- data[,
      c("Peptidoform", "Accession", intensity_cols),
      with = FALSE
    ]
    setnames(
      copf_wide,
      c("Peptidoform", "Accession"),
      c("peptide_id", "protein_id")
    )
    for (col in intensity_cols) {
      set(
        copf_wide,
        j = col,
        value = ifelse(
          is.finite(copf_wide[[col]]),
          2^copf_wide[[col]],
          NA_real_
        )
      )
    }

    # Order fractions: already C1_R1, C1_R2 ... C5_R3 in natural order
    fraction_annot <- data.table(
      filename = intensity_cols,
      fraction_number = seq_along(intensity_cols)
    )

    # importPCPdata (proteoformLocationMapping branch) expects long format:
    # peptide_id, protein_id, filename, intensity
    copf_long <- data.table::melt(
      copf_wide,
      id.vars = c("peptide_id", "protein_id"),
      measure.vars = intensity_cols,
      variable.name = "filename",
      value.name = "intensity",
      variable.factor = FALSE
    )

    # protein_id is already in the long data (passed to importPCPdata); including it
    # in trace_annot would create protein_id.x/protein_id.y in annotateTraces merge.
    trace_annot <- unique(data[, .(peptide_id = Peptidoform)])

    res <- run_copf_pipeline(copf_long, fraction_annot, trace_annot, name)

    if (!is.null(res)) {
      fwrite(res, out_path)
      cat("  Saved to:", out_path, "\n")
    }

    cat("Finished:", name, "| Time:", format_time_diff(Sys.time() - t0), "\n")
    cat(paste0(rep("#", 80), collapse = ""), "\n")
  }
}

################################################################################
## ProteoMaker NA datasets — KNN imputed (k = 5)                              ##
################################################################################

if (run_proteomaker_na_imputed) {
  na_pm_datasets <- list(
    "LowNA"  = "./02-DCF/00-Data/02-Prepared/ProteoMakerLowNA_wide_processed.csv",
    "MedNA"  = "./02-DCF/00-Data/02-Prepared/ProteoMakerMedNA_wide_processed.csv",
    "HighNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerHighNA_wide_processed.csv"
  )

  for (name in names(na_pm_datasets)) {
    file_path <- na_pm_datasets[[name]]
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }

    out_path <- file.path(
      output_dir,
      paste0("COPF_ProteoMaker", name, "_imputed_results.csv")
    )
    if (file.exists(out_path)) {
      cat("Skipping", name, "(imputed) - result already exists\n")
      next
    }

    cat(paste0(rep("#", 80), collapse = ""), "\n")
    cat("Processing ProteoMaker NA dataset (COPF, imputed):", name, "\n")
    t0 <- Sys.time()

    data <- fread(file_path)
    intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

    # Impute on log2 scale, then back-transform to linear for CCprofiler
    imputed_vals <- VIM::kNN(data[, ..intensity_cols], k = 5, imp_var = FALSE)
    non_intensity_cols <- setdiff(names(data), intensity_cols)
    data_imp <- cbind(data[, ..non_intensity_cols], imputed_vals)

    copf_wide <- data_imp[,
      c("Peptidoform", "Accession", intensity_cols),
      with = FALSE
    ]
    setnames(copf_wide, c("Peptidoform", "Accession"), c("peptide_id", "protein_id"))
    for (col in intensity_cols) {
      set(copf_wide, j = col,
          value = ifelse(is.finite(copf_wide[[col]]), 2^copf_wide[[col]], NA_real_))
    }

    fraction_annot <- data.table(
      filename = intensity_cols,
      fraction_number = seq_along(intensity_cols)
    )

    copf_long <- data.table::melt(
      copf_wide,
      id.vars = c("peptide_id", "protein_id"),
      measure.vars = intensity_cols,
      variable.name = "filename",
      value.name = "intensity",
      variable.factor = FALSE
    )

    trace_annot <- unique(data_imp[, .(peptide_id = Peptidoform)])

    res <- run_copf_pipeline(
      copf_long, fraction_annot, trace_annot, paste0(name, "_imputed")
    )

    if (!is.null(res)) {
      fwrite(res, out_path)
      cat("  Saved to:", out_path, "\n")
    }

    cat("Finished:", name, "(imputed) | Time:", format_time_diff(Sys.time() - t0), "\n")
    cat(paste0(rep("#", 80), collapse = ""), "\n")
  }
}

################################################################################
## Protein complex datasets (WideRaw)                                          ##
################################################################################

cx_datasets <- list(
  "B-Wich" = "./02-DCF/00-Data/02-Prepared/B-Wich_WideNoLogNoNorm.csv",
  "CRD-mediated" = "./02-DCF/00-Data/02-Prepared/CRD-mediated_WideNoLogNoNorm.csv",
  "Multiaminoacyl-tRNA" = "./02-DCF/00-Data/02-Prepared/Multiaminoacyl-tRNA_WideNoLogNoNorm.csv",
  "26S-proteasome" = "./02-DCF/00-Data/02-Prepared/26S-proteasome_WideNoLogNoNorm.csv",
  "60S-cytosolic-large-ribosomal-subunit" = "./02-DCF/00-Data/02-Prepared/60S-cytosolic-large-ribosomal-subunit_WideNoLogNoNorm.csv",
  "BSSADCRC" = "./02-DCF/00-Data/02-Prepared/BSSADCRC_WideNoLogNoNorm.csv",
  "Dynactin-complex" = "./02-DCF/00-Data/02-Prepared/Dynactin-complex_WideNoLogNoNorm.csv",
  "Dynein-1-complex-variant-1" = "./02-DCF/00-Data/02-Prepared/Dynein-1-complex-variant-1_WideNoLogNoNorm.csv",
  "Eukaryotic-translation-initiation-factor-3-complex" = "./02-DCF/00-Data/02-Prepared/Eukaryotic-translation-initiation-factor-3-complex_WideNoLogNoNorm.csv",
  "Intraflagellar-transport-complex-B" = "./02-DCF/00-Data/02-Prepared/Intraflagellar-transport-complex-B_WideNoLogNoNorm.csv",
  "Laminin-213" = "./02-DCF/00-Data/02-Prepared/Laminin-213_WideNoLogNoNorm.csv",
  "Major-Spliceosomal-B" = "./02-DCF/00-Data/02-Prepared/Major-Spliceosomal-B_WideNoLogNoNorm.csv",
  "Nuclear-pore-complex" = "./02-DCF/00-Data/02-Prepared/Nuclear-pore-complex_WideNoLogNoNorm.csv"
)

if (run_complex_raw) {
  for (name in names(cx_datasets)) {
    file_path <- cx_datasets[[name]]
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }

    out_path <- file.path(
      output_dir,
      paste0("COPF_Complex_", name, "_results.csv")
    )
    if (file.exists(out_path)) {
      cat("Skipping", name, "- result already exists\n")
      next
    }

    cat(paste0(rep("#", 80), collapse = ""), "\n")
    cat("Processing protein complex dataset (COPF):", name, "\n")
    t0 <- Sys.time()

    data <- fread(file_path)
    intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)

    # Order fractions by condition (D-number) then replicate (B-number)
    condition <- sub("^B\\d+_(D\\d+)$", "\\1", intensity_cols)
    replicate <- sub("^(B\\d+)_D\\d+$", "\\1", intensity_cols)
    col_order <- order(
      as.numeric(gsub("D", "", condition)),
      as.numeric(gsub("B", "", replicate))
    )
    intensity_cols <- intensity_cols[col_order]

    fraction_annot <- data.table(
      filename = intensity_cols,
      fraction_number = seq_along(intensity_cols)
    )

    # importPCPdata (proteoformLocationMapping branch) expects long format:
    # peptide_id, protein_id, filename, intensity
    copf_wide <- data[,
      c("Peptide", "complex_id", intensity_cols),
      with = FALSE
    ]
    setnames(
      copf_wide,
      c("Peptide", "complex_id"),
      c("peptide_id", "protein_id")
    )
    copf_long <- data.table::melt(
      copf_wide,
      id.vars = c("peptide_id", "protein_id"),
      measure.vars = intensity_cols,
      variable.name = "filename",
      value.name = "intensity",
      variable.factor = FALSE
    )

    # protein_id is already in the long data; excluding it here avoids column conflict
    # in annotateTraces merge. For complex datasets, protein_id = complex_id
    # (one complex = one "protein"), so proteoforms found are complexoforms.
    trace_annot <- unique(data[, .(peptide_id = Peptide)])

    res <- run_copf_pipeline(copf_long, fraction_annot, trace_annot, name)

    if (!is.null(res)) {
      fwrite(res, out_path)
      cat("  Saved to:", out_path, "\n")
    }

    cat("Finished:", name, "| Time:", format_time_diff(Sys.time() - t0), "\n")
    cat(paste0(rep("#", 80), collapse = ""), "\n")
  }
}

################################################################################
## Protein complex datasets — KNN imputed (k = 5)                             ##
################################################################################

if (run_complex_imputed) {
  for (name in names(cx_datasets)) {
    file_path <- cx_datasets[[name]]
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }

    out_path <- file.path(
      output_dir,
      paste0("COPF_Complex_", name, "_imputed_results.csv")
    )
    if (file.exists(out_path)) {
      cat("Skipping", name, "(imputed) - result already exists\n")
      next
    }

    cat(paste0(rep("#", 80), collapse = ""), "\n")
    cat("Processing protein complex dataset (COPF, imputed):", name, "\n")
    t0 <- Sys.time()

    data <- fread(file_path)
    intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)

    imputed_vals <- VIM::kNN(data[, ..intensity_cols], k = 5, imp_var = FALSE)
    non_intensity_cols <- setdiff(names(data), intensity_cols)
    data_imp <- cbind(data[, ..non_intensity_cols], imputed_vals)

    condition <- sub("^B\\d+_(D\\d+)$", "\\1", intensity_cols)
    replicate <- sub("^(B\\d+)_D\\d+$", "\\1", intensity_cols)
    col_order <- order(
      as.numeric(gsub("D", "", condition)),
      as.numeric(gsub("B", "", replicate))
    )
    intensity_cols <- intensity_cols[col_order]

    fraction_annot <- data.table(
      filename = intensity_cols,
      fraction_number = seq_along(intensity_cols)
    )

    # importPCPdata (proteoformLocationMapping branch) expects long format:
    # peptide_id, protein_id, filename, intensity
    copf_wide <- data_imp[,
      c("Peptide", "complex_id", intensity_cols),
      with = FALSE
    ]
    setnames(
      copf_wide,
      c("Peptide", "complex_id"),
      c("peptide_id", "protein_id")
    )
    copf_long <- data.table::melt(
      copf_wide,
      id.vars = c("peptide_id", "protein_id"),
      measure.vars = intensity_cols,
      variable.name = "filename",
      value.name = "intensity",
      variable.factor = FALSE
    )

    trace_annot <- unique(data_imp[, .(peptide_id = Peptide)])

    res <- run_copf_pipeline(
      copf_long,
      fraction_annot,
      trace_annot,
      paste0(name, "_imputed")
    )

    if (!is.null(res)) {
      fwrite(res, out_path)
      cat("  Saved to:", out_path, "\n")
    }

    cat(
      "Finished:",
      name,
      "(imputed) | Time:",
      format_time_diff(Sys.time() - t0),
      "\n"
    )
    cat(paste0(rep("#", 80), collapse = ""), "\n")
  }
}
