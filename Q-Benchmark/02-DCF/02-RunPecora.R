options(warn = -1)

suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(PeCorA, warn.conflicts = FALSE))
suppressMessages(library(arrow, warn.conflicts = FALSE))
suppressMessages(library(VIM, warn.conflicts = FALSE))

# Set working directory
setwd(here::here("Q-Benchmark"))

output_dir <- "./02-DCF/00-Data/03-Results/PeCorA"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#####################################
## ProteoMaker datasets for PeCorA ##
#####################################

datasets <- list(
  "LowNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerLowNA_wide_processed.csv",
  "LowNoNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerLowNoNA_wide_processed.csv",
  "MedNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerMedNA_wide_processed.csv",
  "MedNoNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerMedNoNA_wide_processed.csv",
  "HighNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerHighNA_wide_processed.csv",
  "HighNoNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerHighNoNA_wide_processed.csv"
)

for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  output_path <- file.path(output_dir, paste0(name, "_PeCorA_results.csv"))
  if (file.exists(output_path)) {
    cat("Skipping", name, "- result already exists\n")
    next
  }

  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Processing ProteoMaker dataset:", name, "\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")
  data <- fread(file_path)

  # Convert to long format
  intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)
  long_data <- data.table::melt(
    data,
    id.vars = c("Peptidoform", "Accession"),
    measure.vars = intensity_cols,
    variable.name = "Condition_BioReplicate",
    value.name = "Normalized.Area"
  )

  # Extract Condition and BioReplicate from Condition_BioReplicate
  long_data[,
    c("Condition", "BioReplicate") := tstrsplit(
      Condition_BioReplicate,
      "_",
      fixed = TRUE
    )
  ]
  long_data[, Condition_BioReplicate := NULL]

  # NOTE:
  # The prepared data are already log2-transformed
  long_data[,
    Area_for_pecora := ifelse(
      is.finite(Normalized.Area),
      2^Normalized.Area,
      NA_real_
    )
  ]

  # Reorder columns
  setcolorder(
    long_data,
    c(
      "Peptidoform",
      "Accession",
      "Condition",
      "BioReplicate",
      "Normalized.Area",
      "Area_for_pecora"
    )
  )

  pecora_df <- long_data

  pecora_df[, Peptide.Modified.Sequence := Peptidoform]
  pecora_df[, Protein := Accession]

  scaled_peptides <- PeCorA::PeCorA_preprocessing(
    t = pecora_df,
    area_column_name = "Area_for_pecora",
    threshold_to_filter = min(pecora_df$Area_for_pecora, na.rm = TRUE),
    control_name = "C1"
  )

  proteins_testable <- scaled_peptides[,
    .(n_peptides = uniqueN(modpep_z)),
    by = Protein
  ][n_peptides > 1, .N]

  if (proteins_testable == 0) {
    warning(
      paste(
        "No proteins have at least 2 peptides after preprocessing in dataset",
        name,
        "; PeCorA cannot run."
      )
    )
    next
  }

  # Run PeCorA
  res <- PeCorA::PeCorA(
    scaled_peptides
  )

  fwrite(res, output_path)

  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Finished processing dataset:", name, "\n")
  cat(paste0("Saved PeCorA results to:", output_path), "\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")
}

##############################################################
## ProteoMaker NA datasets for PeCorA — kNN imputed (k = 5) ##
##############################################################

na_datasets <- list(
  "LowNA"  = "./02-DCF/00-Data/02-Prepared/ProteoMakerLowNA_wide_processed.csv",
  "MedNA"  = "./02-DCF/00-Data/02-Prepared/ProteoMakerMedNA_wide_processed.csv",
  "HighNA" = "./02-DCF/00-Data/02-Prepared/ProteoMakerHighNA_wide_processed.csv"
)

for (name in names(na_datasets)) {
  file_path <- na_datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  output_path <- file.path(output_dir, paste0(name, "_imputed_PeCorA_results.csv"))
  if (file.exists(output_path)) {
    cat("Skipping", name, "(imputed) - result already exists\n")
    next
  }

  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Processing ProteoMaker NA dataset with kNN imputation:", name, "\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")

  data <- fread(file_path)
  intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

  # Impute missing values on the log2 scale then back-transform for PeCorA
  imputed_vals <- VIM::kNN(data[, ..intensity_cols], k = 5, imp_var = FALSE)
  non_intensity_cols <- setdiff(names(data), intensity_cols)
  data_imp <- cbind(data[, ..non_intensity_cols], imputed_vals)

  long_data <- data.table::melt(
    data_imp,
    id.vars = c("Peptidoform", "Accession"),
    measure.vars = intensity_cols,
    variable.name = "Condition_BioReplicate",
    value.name = "Normalized.Area"
  )
  long_data[,
    c("Condition", "BioReplicate") := tstrsplit(
      Condition_BioReplicate, "_", fixed = TRUE
    )
  ]
  long_data[, Condition_BioReplicate := NULL]
  long_data[,
    Area_for_pecora := ifelse(
      is.finite(Normalized.Area), 2^Normalized.Area, NA_real_
    )
  ]
  setcolorder(
    long_data,
    c("Peptidoform", "Accession", "Condition", "BioReplicate",
      "Normalized.Area", "Area_for_pecora")
  )

  pecora_df <- long_data
  pecora_df[, Peptide.Modified.Sequence := Peptidoform]
  pecora_df[, Protein := Accession]

  scaled_peptides <- PeCorA::PeCorA_preprocessing(
    t = pecora_df,
    area_column_name = "Area_for_pecora",
    threshold_to_filter = min(pecora_df$Area_for_pecora, na.rm = TRUE),
    control_name = "C1"
  )

  proteins_testable <- scaled_peptides[,
    .(n_peptides = uniqueN(modpep_z)), by = Protein
  ][n_peptides > 1, .N]

  if (proteins_testable == 0) {
    warning(paste(
      "No proteins have >= 2 peptides after preprocessing in imputed dataset",
      name, "; PeCorA cannot run."
    ))
    next
  }

  res <- tryCatch(
    PeCorA::PeCorA(scaled_peptides),
    error = function(e) {
      warning(paste("PeCorA failed for imputed dataset", name, ":", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(res)) next

  fwrite(res, output_path)
  cat("Saved:", output_path, "\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")
}

#########################################
## Protein complex datasets for PeCorA ##
#########################################

datasets <- list(
  "B-Wich"               = "./02-DCF/00-Data/02-Prepared/B-Wich_WideNoLogNoNorm.csv",
  "CRD-mediated"         = "./02-DCF/00-Data/02-Prepared/CRD-mediated_WideNoLogNoNorm.csv",
  "Multiaminoacyl-tRNA"  = "./02-DCF/00-Data/02-Prepared/Multiaminoacyl-tRNA_WideNoLogNoNorm.csv",
  "26S-proteasome"       = "./02-DCF/00-Data/02-Prepared/26S-proteasome_WideNoLogNoNorm.csv",
  "60S-cytosolic-large-ribosomal-subunit" = "./02-DCF/00-Data/02-Prepared/60S-cytosolic-large-ribosomal-subunit_WideNoLogNoNorm.csv",
  "BSSADCRC"             = "./02-DCF/00-Data/02-Prepared/BSSADCRC_WideNoLogNoNorm.csv",
  "Dynactin-complex"     = "./02-DCF/00-Data/02-Prepared/Dynactin-complex_WideNoLogNoNorm.csv",
  "Dynein-1-complex-variant-1" = "./02-DCF/00-Data/02-Prepared/Dynein-1-complex-variant-1_WideNoLogNoNorm.csv",
  "Eukaryotic-translation-initiation-factor-3-complex" = "./02-DCF/00-Data/02-Prepared/Eukaryotic-translation-initiation-factor-3-complex_WideNoLogNoNorm.csv",
  "Intraflagellar-transport-complex-B" = "./02-DCF/00-Data/02-Prepared/Intraflagellar-transport-complex-B_WideNoLogNoNorm.csv",
  "Laminin-213"          = "./02-DCF/00-Data/02-Prepared/Laminin-213_WideNoLogNoNorm.csv",
  "Major-Spliceosomal-B" = "./02-DCF/00-Data/02-Prepared/Major-Spliceosomal-B_WideNoLogNoNorm.csv",
  "Nuclear-pore-complex" = "./02-DCF/00-Data/02-Prepared/Nuclear-pore-complex_WideNoLogNoNorm.csv"
)

for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  output_path <- file.path(output_dir, paste0(name, "_PeCorA_results.csv"))
  if (file.exists(output_path)) {
    cat("Skipping", name, "- result already exists\n")
    next
  }

  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Processing protein complex dataset:", name, "\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")

  data <- fread(file_path)

  # Convert to long format
  intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)
  long_data <- data.table::melt(
    data,
    id.vars = c("Peptide", "complex_id"),
    measure.vars = intensity_cols,
    variable.name = "BioReplicate_Condition",
    value.name = "Area_for_pecora"
  )

  # Extract BioReplicate and Condition from BioReplicate_Condition
  long_data[,
    c("BioReplicate", "Condition") := tstrsplit(
      BioReplicate_Condition,
      "_",
      fixed = TRUE
    )
  ]
  long_data[, BioReplicate_Condition := NULL]

  pecora_df <- long_data

  pecora_df[, Peptide.Modified.Sequence := Peptide]
  pecora_df[, Protein := complex_id]

  scaled_peptides <- PeCorA::PeCorA_preprocessing(
    t = pecora_df,
    area_column_name = "Area_for_pecora",
    threshold_to_filter = min(pecora_df$Area_for_pecora, na.rm = TRUE),
    control_name = "D17"
  )

  proteins_testable <- scaled_peptides[,
    .(n_peptides = uniqueN(modpep_z)),
    by = Protein
  ][n_peptides > 1, .N]

  if (proteins_testable == 0) {
    warning(
      paste(
        "No proteins have at least 2 peptides after preprocessing in dataset",
        name,
        "; PeCorA cannot run."
      )
    )
    next
  }

  # Run PeCorA
  res <- PeCorA::PeCorA(
    scaled_peptides
  )

  fwrite(res, output_path)
}

#################################################################
## Protein complex datasets for PeCorA - with imputation first ##
#################################################################

datasets <- list(
  "B-Wich"               = "./02-DCF/00-Data/02-Prepared/B-Wich_WideNoLogNoNorm.csv",
  "CRD-mediated"         = "./02-DCF/00-Data/02-Prepared/CRD-mediated_WideNoLogNoNorm.csv",
  "Multiaminoacyl-tRNA"  = "./02-DCF/00-Data/02-Prepared/Multiaminoacyl-tRNA_WideNoLogNoNorm.csv",
  "26S-proteasome"       = "./02-DCF/00-Data/02-Prepared/26S-proteasome_WideNoLogNoNorm.csv",
  "60S-cytosolic-large-ribosomal-subunit" = "./02-DCF/00-Data/02-Prepared/60S-cytosolic-large-ribosomal-subunit_WideNoLogNoNorm.csv",
  "BSSADCRC"             = "./02-DCF/00-Data/02-Prepared/BSSADCRC_WideNoLogNoNorm.csv",
  "Dynactin-complex"     = "./02-DCF/00-Data/02-Prepared/Dynactin-complex_WideNoLogNoNorm.csv",
  "Dynein-1-complex-variant-1" = "./02-DCF/00-Data/02-Prepared/Dynein-1-complex-variant-1_WideNoLogNoNorm.csv",
  "Eukaryotic-translation-initiation-factor-3-complex" = "./02-DCF/00-Data/02-Prepared/Eukaryotic-translation-initiation-factor-3-complex_WideNoLogNoNorm.csv",
  "Intraflagellar-transport-complex-B" = "./02-DCF/00-Data/02-Prepared/Intraflagellar-transport-complex-B_WideNoLogNoNorm.csv",
  "Laminin-213"          = "./02-DCF/00-Data/02-Prepared/Laminin-213_WideNoLogNoNorm.csv",
  "Major-Spliceosomal-B" = "./02-DCF/00-Data/02-Prepared/Major-Spliceosomal-B_WideNoLogNoNorm.csv",
  "Nuclear-pore-complex" = "./02-DCF/00-Data/02-Prepared/Nuclear-pore-complex_WideNoLogNoNorm.csv"
)

imputed_datasets <- list()

# Impute missing values with knn = 5 before running PeCorA
for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  imputed_out <- file.path(output_dir, paste0(name, "_imputed_PeCorA_results.csv"))
  if (file.exists(imputed_out)) {
    cat("Skipping imputation for", name, "- result already exists\n")
    next
  }

  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Processing protein complex dataset with imputation:", name, "\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")

  data <- fread(file_path)
  dataset_intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)

  # Impute missing values using knn = 5
  imputed_data <- VIM::kNN(
    data[, ..dataset_intensity_cols],
    k = 5,
    imp_var = FALSE
  )

  # Combine imputed data with non-intensity columns
  non_intensity_cols <- setdiff(names(data), dataset_intensity_cols)
  imputed_data <- cbind(data[, ..non_intensity_cols], imputed_data)

  imputed_datasets[[name]] <- imputed_data
}

# Now run PeCorA on the imputed datasets
for (name in names(imputed_datasets)) {
  data <- imputed_datasets[[name]]

  imputed_out <- file.path(output_dir, paste0(name, "_imputed_PeCorA_results.csv"))
  if (file.exists(imputed_out)) {
    cat("Skipping", name, "- imputed result already exists\n")
    next
  }

  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Running PeCorA on imputed protein complex dataset:", name, "\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")

  # Convert to long format
  intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)
  long_data <- data.table::melt(
    data,
    id.vars = c("Peptide", "complex_id"),
    measure.vars = intensity_cols,
    variable.name = "BioReplicate_Condition",
    value.name = "Area_for_pecora"
  )

  # Extract BioReplicate and Condition from BioReplicate_Condition
  long_data[,
    c("BioReplicate", "Condition") := tstrsplit(
      BioReplicate_Condition,
      "_",
      fixed = TRUE
    )
  ]
  long_data[, BioReplicate_Condition := NULL]

  pecora_df <- long_data
  pecora_df[, Peptide.Modified.Sequence := Peptide]
  pecora_df[, Protein := complex_id]

  scaled_peptides <- PeCorA::PeCorA_preprocessing(
    t = pecora_df,
    area_column_name = "Area_for_pecora",
    threshold_to_filter = min(pecora_df$Area_for_pecora, na.rm = TRUE),
    control_name = "D17"
  )

  proteins_testable <- scaled_peptides[,
    .(n_peptides = uniqueN(modpep_z)),
    by = Protein
  ][n_peptides > 1, .N]

  if (proteins_testable == 0) {
    warning(
      paste(
        "No proteins have at least 2 peptides after preprocessing in imputed dataset",
        name,
        "; PeCorA cannot run."
      )
    )
    next
  }

  res <- tryCatch(
    PeCorA::PeCorA(scaled_peptides),
    error = function(e) {
      warning(
        paste(
          "PeCorA failed for imputed dataset",
          name,
          "with error:",
          conditionMessage(e)
        )
      )
      NULL
    }
  )

  if (is.null(res)) {
    next
  }

  fwrite(res, imputed_out)

  cat("Saved imputed PeCorA results to:", imputed_out, "\n")
}
