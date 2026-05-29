options(warn = -1)

library(data.table, warn.conflicts = FALSE)

setwd(here::here())

source("ComplexoFinder/02-Identification/AssignComplexes.R")

start_time <- Sys.time()

output_dir <- "./Q-Benchmark/02-DCF/00-Data/03-Results/DCFBaseline"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Baseline normalization: for each peptide (row), subtract the median of its
# first available condition so all profiles start at zero.  Works in log-space,
# which is equivalent to dividing by the baseline in linear space.
baseline_normalize <- function(data, intensity_cols, cond_regex) {
  data <- copy(data)
  cond_from_cols <- sub(cond_regex, "\\1", intensity_cols, perl = TRUE)
  unique_conds <- unique(cond_from_cols)
  unique_conds <- unique_conds[order(as.numeric(gsub(
    "[^0-9]",
    "",
    unique_conds
  )))]

  for (i in seq_len(nrow(data))) {
    vals <- as.numeric(data[i, intensity_cols, with = FALSE])
    baseline_val <- NA_real_
    for (cond in unique_conds) {
      cond_cols <- intensity_cols[cond_from_cols == cond]
      cond_vals <- vals[match(cond_cols, intensity_cols)]
      if (any(!is.na(cond_vals))) {
        baseline_val <- median(cond_vals, na.rm = TRUE)
        break
      }
    }
    if (!is.na(baseline_val)) {
      for (col in intensity_cols) {
        if (!is.na(data[[col]][i])) {
          data[[col]][i] <- data[[col]][i] - baseline_val
        }
      }
    }
  }
  data
}

################################################################################
## ProteoMaker datasets                                                        ##
################################################################################

pm_datasets <- list(
  "LowNA" = "./Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerLowNA_wide_processed.csv",
  "LowNoNA" = "./Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerLowNoNA_wide_processed.csv",
  "MedNA" = "./Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerMedNA_wide_processed.csv",
  "MedNoNA" = "./Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerMedNoNA_wide_processed.csv",
  "HighNA" = "./Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerHighNA_wide_processed.csv",
  "HighNoNA" = "./Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerHighNoNA_wide_processed.csv"
)

PM_summary <- data.table(
  Dataset = names(pm_datasets),
  Num_Peptides = integer(length(pm_datasets)),
  Num_Discordant_Peptides = integer(length(pm_datasets)),
  Time_Taken = numeric(length(pm_datasets))
)

for (name in names(pm_datasets)) {
  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Processing ProteoMaker dataset (baseline):", name, "\n")

  dataset_start_time <- Sys.time()

  file_path <- pm_datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  out_path <- file.path(
    output_dir,
    paste0("DCFBaseline_ProteoMaker", name, "_results.csv")
  )
  if (file.exists(out_path)) {
    cat("Skipping", name, "- result already exists\n")
    next
  }

  data <- fread(file_path)
  PM_summary[Dataset == name, Num_Peptides := nrow(data)]

  intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

  data[,
    Peptide_ID := paste(
      Peptide,
      Position,
      Proteoform_ID,
      Peptidoform,
      sep = "_"
    )
  ]
  data <- data[, c("Peptide_ID", "Accession", intensity_cols), with = FALSE]

  all_res <- data.table()

  for (accession in unique(data$Accession)) {
    acc_data <- data[Accession == accession]
    if (nrow(acc_data) < 2) {
      next
    }

    acc_norm <- baseline_normalize(
      acc_data,
      intensity_cols,
      cond_regex = "^(C\\d+)_R\\d+$"
    )

    dcf_result <- tryCatch(
      discover_complexoforms(
        complex_data = acc_norm,
        intensity_cols = intensity_cols,
        group_col = "Accession",
        peptide_col = "Peptide_ID",
        alpha = 0.05,
        adjust_method = "BH",
        deep_split = 2,
        minClusterSize = 2,
        cond_regex = "^(C\\d+)_R\\d+$",
        min_total_non_na_frac = 0.4,
        min_conditions = 2,
        min_reps_per_condition = 2,
        canonical_label = "dCF0",
        singleton_label = "dCF-1",
        return_condition_tests = FALSE,
        verbose = FALSE
      ),
      error = function(e) {
        warning(paste(
          "Error processing accession",
          accession,
          ":",
          conditionMessage(e)
        ))
        NULL
      }
    )
    if (is.null(dcf_result)) {
      next
    }
    all_res <- rbind(all_res, dcf_result$peptide_results)
  }

  PM_summary[
    Dataset == name,
    Num_Discordant_Peptides := sum(all_res$discordant)
  ]

  fwrite(
    all_res,
    file.path(
      output_dir,
      paste0("DCFBaseline_ProteoMaker", name, "_results.csv")
    )
  )

  time_taken <- Sys.time() - dataset_start_time
  PM_summary[
    Dataset == name,
    Time_Taken := as.numeric(time_taken, units = "secs")
  ]
  cat("Finished:", name, "| Time:", round(time_taken, 2), "seconds\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")
}

fwrite(PM_summary, file.path(output_dir, "Summary_DCFBaseline_ProteoMaker.csv"))

################################################################################
## Protein complex datasets                                                    ##
################################################################################

cx_datasets <- list(
  "B-Wich" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/B-Wich_WideLogNorm.csv",
  "CRD-mediated" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/CRD-mediated_WideLogNorm.csv",
  "Multiaminoacyl-tRNA" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Multiaminoacyl-tRNA_WideLogNorm.csv",
  "Major-Spliceosomal-B" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Major-Spliceosomal-B_WideLogNorm.csv",
  "26S-proteasome" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/26S-proteasome_WideLogNorm.csv",
  "60S-cytosolic-large-ribosomal-subunit" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/60S-cytosolic-large-ribosomal-subunit_WideLogNorm.csv",
  "Nuclear-pore-complex" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Nuclear-pore-complex_WideLogNorm.csv",
  "Eukaryotic-translation-initiation-factor-3-complex" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Eukaryotic-translation-initiation-factor-3-complex_WideLogNorm.csv",
  "Dynein-1-complex-variant-1" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Dynein-1-complex-variant-1_WideLogNorm.csv",
  "Intraflagellar-transport-complex-B" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Intraflagellar-transport-complex-B_WideLogNorm.csv",
  "BSSADCRC" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/BSSADCRC_WideLogNorm.csv",
  "Dynactin-complex" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Dynactin-complex_WideLogNorm.csv",
  "Laminin-213" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Laminin-213_WideLogNorm.csv"
)

cx_summary <- data.table(
  Dataset = names(cx_datasets),
  Num_Peptides = integer(length(cx_datasets)),
  Num_Discordant_Peptides = integer(length(cx_datasets)),
  Num_dCF = integer(length(cx_datasets)),
  Time_Taken = numeric(length(cx_datasets))
)

for (name in names(cx_datasets)) {
  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Processing protein complex dataset (baseline):", name, "\n")

  dataset_start_time <- Sys.time()

  file_path <- cx_datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  out_path <- file.path(
    output_dir,
    paste0("DCFBaseline_Complex_", name, "_results.csv")
  )
  if (file.exists(out_path)) {
    cat("Skipping", name, "- result already exists\n")
    next
  }

  data <- fread(file_path)
  cx_summary[Dataset == name, Num_Peptides := nrow(data)]

  intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)

  data_norm <- baseline_normalize(
    data,
    intensity_cols,
    cond_regex = "^B\\d+_(D\\d+)$"
  )

  dcf_result <- tryCatch(
    discover_complexoforms(
      complex_data = data_norm,
      intensity_cols = intensity_cols,
      group_col = "complex_id",
      peptide_col = "Peptide",
      alpha = 0.05,
      adjust_method = "BH",
      deep_split = 2,
      minClusterSize = 2,
      cond_regex = "^B\\d+_(D\\d+)$",
      min_total_non_na_frac = 0.4,
      min_conditions = 2,
      min_reps_per_condition = 2,
      canonical_label = "dCF0",
      singleton_label = "dCF-1",
      return_condition_tests = FALSE,
      verbose = FALSE
    ),
    error = function(e) {
      warning(paste("Error processing complex", name, ":", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(dcf_result)) {
    next
  }

  res <- dcf_result$peptide_results
  unique_dCF <- unique(res$dCF)

  cx_summary[Dataset == name, Num_Discordant_Peptides := sum(res$discordant)]
  cx_summary[
    Dataset == name,
    Num_dCF := length(unique_dCF[!is.na(unique_dCF) & unique_dCF != "dCF-1"])
  ]

  fwrite(
    res,
    file.path(output_dir, paste0("DCFBaseline_Complex_", name, "_results.csv"))
  )

  time_taken <- Sys.time() - dataset_start_time
  cx_summary[
    Dataset == name,
    Time_Taken := as.numeric(time_taken, units = "secs")
  ]
  cat("Finished:", name, "| Time:", round(time_taken, 2), "seconds\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")
}

fwrite(cx_summary, file.path(output_dir, "Summary_DCFBaseline_Complex.csv"))

end_time <- Sys.time()
cat(
  "\nTotal time:",
  round(difftime(end_time, start_time, units = "secs"), 2),
  "seconds\n"
)
