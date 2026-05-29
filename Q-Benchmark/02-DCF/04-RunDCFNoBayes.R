options(warn = -1)

# Libraries
suppressMessages(library(data.table, warn.conflicts = FALSE))

setwd("/Users/malthe/Desktop/0MastersProject/9")

# Source functions
source("quantification_complex/02-Identification/RunPoolComplex.R")

start_time <- Sys.time()

output_dir <- "Q-Benchmark/02-DCF/00-Data/03-Results/DCFNoBayes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

##################################
## ProteoMaker datasets for RPC ##
##################################
datasets <- list(
  "LowNA" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerLowNA_wide_processed.csv",
  "LowNoNA" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerLowNoNA_wide_processed.csv",
  "MedNA" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerMedNA_wide_processed.csv",
  "MedNoNA" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerMedNoNA_wide_processed.csv",
  "HighNA" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerHighNA_wide_processed.csv",
  "HighNoNA" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/ProteoMakerHighNoNA_wide_processed.csv"
)

# Create data table to store #peptides, #discordant peptides, and time
PM_summary <- data.table(
  Dataset = names(datasets),
  Num_Peptides = integer(length(datasets)),
  Num_Discordant_Peptides = integer(length(datasets)),
  Time_Taken = numeric(length(datasets))
)

for (name in names(datasets)) {
  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Processing ProteoMaker dataset:", name, "\n")

  # Record start time for this dataset
  dataset_start_time <- Sys.time()

  # Read data
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  data <- fread(file_path)

  # Record number of peptides in summary table
  PM_summary[Dataset == name, Num_Peptides := nrow(data)]

  # Get intensity columns
  intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

  # Combine non-intensity columns into a single "Peptidoform" column
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

  # Find unique Accessions in data
  unique_accessions <- unique(data$Accession)

  all_res <- data.table()

  # Run RPC for each accession and combine results
  for (accession in unique_accessions) {
    acc_data <- data[Accession == accession]
    if (nrow(acc_data) < 2) {
      next
    }

    PM_DCFNoBayes_result <- tryCatch(
      {
        discover_complexoforms_no_bayes(
          complex_data = acc_data,
          intensity_cols = intensity_cols,
          group_col = "Accession",
          peptide_col = "Peptide_ID",
          alpha = 0.05,
          deep_split = 2,
          adjust_method = "BH",
          cond_regex = "^C(\\d+)_R(\\d+)$",
          min_total_non_na_frac = 0.4,
          min_conditions = 2,
          min_reps_per_condition = 2,
          canonical_label = "dCF0",
          singleton_label = "dCF-1",
          return_condition_tests = FALSE,
          verbose = FALSE
        )
      }
    )
    if (is.null(PM_DCFNoBayes_result)) {
      next
    }
    all_res <- rbind(all_res, PM_DCFNoBayes_result$peptide_results)
  }
  # Count number of discordant peptides (those with adjusted p-value < 0.01)
  num_discordant <- sum(all_res$discordant)
  PM_summary[Dataset == name, Num_Discordant_Peptides := num_discordant]

  # Write results for this dataset
  output_path <- paste0(
    output_dir,
    "/DCFNoBayes_ProteoMaker",
    name,
    "_results.csv"
  )
  fwrite(all_res, output_path)

  # Record time taken for DCFNoBayes for this dataset
  dataset_end_time <- Sys.time()
  time_taken <- dataset_end_time - dataset_start_time

  PM_summary[
    Dataset == name,
    Time_Taken := as.numeric(time_taken, units = "secs")
  ]

  cat(
    "Finished processing dataset:",
    name,
    "\n",
    "Time taken for RPC:",
    round(time_taken, 2),
    "seconds\n"
  )
  cat(paste0(rep("#", 80), collapse = ""), "\n")
}

# Write summary table to CSV
summary_output_path <- paste0(output_dir, "/Summary_DCFNoBayes_ProteoMaker.csv")
fwrite(PM_summary, summary_output_path)

#############################################
## Protein complex datasets for DCFNoBayes ##
#############################################

datasets <- list(
  "CRD-mediated" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/CRD-mediated_WideLogNorm.csv",
  "Multiaminoacyl-tRNA" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Multiaminoacyl-tRNA_WideLogNorm.csv",
  "B-Wich" = "Q-Benchmark/02-DCF/00-Data/02-Prepared/B-Wich_WideLogNorm.csv",
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

complex_summary <- data.table(
  Dataset = names(datasets),
  Num_Peptides = integer(length(datasets)),
  Num_Discordant_Peptides = integer(length(datasets)),
  Num_dCF = integer(length(datasets)),
  Time_Taken = numeric(length(datasets))
)

for (name in names(datasets)) {
  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Processing protein complex dataset:", name, "\n")

  # Record start time for this dataset
  dataset_start_time <- Sys.time()

  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  data <- fread(file_path)

  num_peptides <- nrow(data)
  complex_summary[Dataset == name, Num_Peptides := num_peptides]

  intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)

  PC_DCFNoBayes_result <- tryCatch(
    {
      discover_complexoforms_no_bayes(
        complex_data = data,
        intensity_cols = intensity_cols,
        group_col = "complex_id",
        peptide_col = "Peptide",
        alpha = 0.05,
        adjust_method = "BH",
        cond_regex = "^B(\\d+)_D(\\d+)$",
        min_total_non_na_frac = 0.4,
        min_conditions = 2,
        min_reps_per_condition = 2,
        canonical_label = "dCF0",
        singleton_label = "dCF-1",
        return_condition_tests = FALSE,
        verbose = FALSE
      )
    }
  )
  if (is.null(PC_DCFNoBayes_result)) {
    next
  }
  res <- PC_DCFNoBayes_result$peptide_results

  num_discordant <- sum(res$discordant)
  unique_dCF <- unique(res$dCF)
  num_dCF <- length(unique_dCF[!is.na(unique_dCF) & unique_dCF != "dCF-1"])
  complex_summary[Dataset == name, Num_Discordant_Peptides := num_discordant]
  complex_summary[Dataset == name, Num_dCF := num_dCF]

  output_path <- paste0(
    output_dir,
    "/DCFNoBayes_Complex_",
    name,
    "_results.csv"
  )
  fwrite(res, output_path)

  dataset_end_time <- Sys.time()
  time_taken <- dataset_end_time - dataset_start_time
  complex_summary[
    Dataset == name,
    Time_Taken := as.numeric(time_taken, units = "secs")
  ]
  cat(
    "Finished processing dataset:",
    name,
    "\n",
    "Time taken for RPC:",
    round(time_taken, 2),
    "seconds\n"
  )
  cat(paste0(rep("#", 80), collapse = ""), "\n")
}

# Write summary table to CSV
complex_summary_output_path <- paste0(
  output_dir,
  "/Summary_DCFNoBayes_Complex.csv"
)
fwrite(complex_summary, complex_summary_output_path)
