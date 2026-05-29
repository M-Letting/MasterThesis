# Libraries
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Set working directory to project root
setwd(here::here())

# Source functions
source("ComplexoFinder/03-Constraints/Constraints.R")

EBI_LT <- readRDS("ComplexoFinder/00-LookupTables/ebi_cp_lookup.rds")

# Set working directory to Q-Benchmark
setwd(paste0(getwd(), "/Q-Benchmark"))

# Set output directory for results
output_dir <- "03-Clustering/00-Data/02-CannotLinkMatrices"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cannot_link_variant <- getOption("qbenchmark.cannot_link_variant", "corr")
if (!cannot_link_variant %in% c("corr", "rmsd", "ccc")) {
  stop("Unsupported cannot-link variant: ", cannot_link_variant)
}

create_cannot_link_fun <- switch(
  cannot_link_variant,
  corr = create_cannot_link_Pearson,
  rmsd = create_cannot_link_ZRMSD,
  ccc = create_cannot_link_CCC
)

################################################################################
## Generate Cannot-Link matrices for all datasets ##############################
################################################################################

start_time <- Sys.time()

get_complex_rpc_datasets <- function(
  method_prefix,
  rpc_dir = "./02-DCF/00-Data/03-Results/DCF"
) {
  pattern <- paste0("^", method_prefix, "_Complex_.*_results\\.csv$")
  files <- sort(list.files(rpc_dir, pattern = pattern, full.names = TRUE))

  if (length(files) == 0L) {
    stop(paste0(
      "No complex datasets found for prefix '",
      method_prefix,
      "' in ",
      rpc_dir
    ))
  }

  dataset_names <- sub(
    "_results\\.csv$",
    "",
    sub("^.*_Complex_", "", basename(files))
  )
  setNames(files, dataset_names)
}

##########################
## ProteoMaker datasets ##
##########################

cat(paste0(rep("#", 80), collapse = ""), "\n")
cat(
  "Generating Cannot-Link matrices for ProteoMaker datasets (variant:",
  cannot_link_variant,
  ")\n"
)

PM_start_time <- Sys.time()

# Define dataset paths
datasets <- list(
  "LowNA" = "02-DCF/00-Data/02-Prepared/ProteoMakerLowNA_wide_processed.csv",
  "LowNoNA" = "02-DCF/00-Data/02-Prepared/ProteoMakerLowNoNA_wide_processed.csv",
  "MedNA" = "02-DCF/00-Data/02-Prepared/ProteoMakerMedNA_wide_processed.csv",
  "MedNoNA" = "02-DCF/00-Data/02-Prepared/ProteoMakerMedNoNA_wide_processed.csv",
  "HighNA" = "02-DCF/00-Data/02-Prepared/ProteoMakerHighNA_wide_processed.csv",
  "HighNoNA" = "02-DCF/00-Data/02-Prepared/ProteoMakerHighNoNA_wide_processed.csv"
)

# Get accessions for top 25 proteins with most peptides
high_no_na_data <- fread(datasets[["HighNoNA"]])
top_proteins_high_no_na <- high_no_na_data[,
  .N,
  by = Accession
][order(-N)][1:25, Accession]

# Set multiple cannot-link thresholds to test
cannot_link_thresholds <- c(0.4, 0.5, 0.6, 0.7, 0.8)

# Initialize summary table for logging results
PM_summary_table <- data.table(
  CannotLinkMethod = character(),
  Dataset = character(),
  Threshold = numeric(),
  Protein = character(),
  nPeptides = integer(),
  nRestrictions = integer(),
  TimeTaken = numeric()
)

# Loop through thresholds and datasets to create Cannot-Link matrices
for (th in cannot_link_thresholds) {
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat("Processing with cannot-link threshold:", th, "\n")

  for (name in names(datasets)) {
    cat(paste0(rep("-", 40), collapse = ""), "\n")
    cat("Processing dataset:", name, "\n")

    # Read data
    file_path <- datasets[[name]]
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }

    data <- fread(file_path)

    # Get intensity columns for this dataset (constant across proteins)
    intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

    # Subset to peptides from top 10 proteins in HighNoNA dataset
    for (protein in top_proteins_high_no_na) {
      protein_start_time <- Sys.time()

      data_protein <- data[Accession == protein]
      if (nrow(data_protein) == 0) {
        next
      }

      # Create Cannot-Link matrix for this protein
      CLM <- create_cannot_link_fun(
        data = data_protein,
        complex_name = NULL,
        id_cols = c(
          "Peptide",
          "Accession",
          "Position",
          "Proteoform_ID",
          "Peptidoform"
        ),
        identifier_info = c(
          "Peptide",
          "Accession",
          "Position",
          "Proteoform_ID",
          "Peptidoform"
        ),
        intensity_cols = intensity_cols,
        n_rep = 3,
        n_cond = 10,
        min_shared_cond = 2,
        cannot_link_th = th
      )

      # Save Cannot-Link matrix
      th_no_dot <- gsub("\\.", "", th) # Remove dots from threshold for file naming
      output_file <- file.path(
        output_dir,
        paste0(
          "/PM-",
          name,
          "_",
          protein,
          "_CLM",
          th_no_dot,
          "_",
          cannot_link_variant,
          ".csv"
        )
      )
      write.csv(CLM, output_file, row.names = TRUE)

      # Log results in summary table
      n_peptides <- nrow(data_protein)
      n_restrictions <- sum(CLM) / 2 # Each pair is counted twice in the matrix
      time_taken <- round(
        as.numeric(difftime(Sys.time(), protein_start_time, units = "secs")),
        2
      )

      PM_summary_table <- rbind(
        PM_summary_table,
        data.table(
          CannotLinkMethod = cannot_link_variant,
          Dataset = name,
          Threshold = th,
          Protein = protein,
          nPeptides = n_peptides,
          nRestrictions = n_restrictions,
          TimeTaken = time_taken
        )
      )
    }
  }
}

# Save summary table
summary_output_file <- file.path(output_dir, "/00-PM_CannotLinkSummary.csv")
fwrite(PM_summary_table, summary_output_file)

summary_output_file <- file.path(
  output_dir,
  paste0("/00-PM_CannotLinkSummary_", cannot_link_variant, ".csv")
)
fwrite(PM_summary_table, summary_output_file)

# Log total time taken for ProteoMaker datasets
PM_end_time <- Sys.time()
PM_time_taken <- round(
  as.numeric(difftime(PM_end_time, PM_start_time, units = "mins")),
  2
)
cat(
  "Finished generating Cannot-Link matrices for ProteoMaker datasets in",
  PM_time_taken,
  "minutes\n"
)
cat(paste0(rep("#", 80), collapse = ""), "\n")

##############################
## Protein complex datasets ##
##############################

cat(
  "Generating Cannot-Link matrices for Protein complex datasets (variant:",
  cannot_link_variant,
  ")\n"
)

PC_start_time <- Sys.time()

# Define dataset paths from expanded DCF roster
datasets <- get_complex_rpc_datasets("DCF")

# Set multiple cannot-link thresholds to test
cannot_link_thresholds <- c(0.4, 0.5, 0.6, 0.7, 0.8)

# Initialize summary table for logging results
PC_summary_table <- data.table(
  CannotLinkMethod = character(),
  Dataset = character(),
  Threshold = numeric(),
  Protein = character(),
  nPeptides = integer(),
  nRestrictions = integer(),
  TimeTaken = numeric()
)


# Loop through thresholds and datasets to create Cannot-Link matrices
for (th in cannot_link_thresholds) {
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat("Processing with cannot-link threshold:", th, "\n")

  for (name in names(datasets)) {
    cat(paste0(rep("-", 40), collapse = ""), "\n")
    cat("Processing dataset:", name, "\n")

    complex_start_time <- Sys.time()

    # Read data
    file_path <- datasets[[name]]
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }

    data <- fread(file_path)

    # Get intensity columns for this dataset (constant across proteins)
    intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)

    # Create Cannot-Link matrix for this dataset
    CLM <- create_cannot_link_fun(
      data = data,
      complex_name = NULL,
      id_cols = c("complex_id", "Peptide"),
      identifier_info = c(
        "complex_id",
        "Gene",
        "Accession",
        "Position",
        "Datatype"
      ),
      intensity_cols = intensity_cols,
      n_rep = 3,
      n_cond = 28,
      min_shared_cond = NULL,
      cannot_link_th = th
    )

    # Save Cannot-Link matrix
    th_no_dot <- gsub("\\.", "", th) # Remove dots from threshold for file naming
    output_file <- file.path(
      output_dir,
      paste0(
        "/PC-",
        name,
        "_CLM",
        th_no_dot,
        "_",
        cannot_link_variant,
        ".csv"
      )
    )
    write.csv(CLM, output_file, row.names = TRUE)

    # Log results in summary table
    n_peptides <- nrow(data)
    n_restrictions <- sum(CLM) / 2 # Each pair is counted twice in the matrix
    time_taken <- round(
      as.numeric(difftime(Sys.time(), complex_start_time, units = "secs")),
      2
    )

    PC_summary_table <- rbind(
      PC_summary_table,
      data.table(
        CannotLinkMethod = cannot_link_variant,
        Dataset = name,
        Threshold = th,
        Protein = "All",
        nPeptides = n_peptides,
        nRestrictions = n_restrictions,
        TimeTaken = time_taken
      )
    )
  }
}

# Save summary table
summary_output_file <- file.path(output_dir, "/00-PC_CannotLinkSummary.csv")
fwrite(PC_summary_table, summary_output_file)

summary_output_file <- file.path(
  output_dir,
  paste0("/00-PC_CannotLinkSummary_", cannot_link_variant, ".csv")
)
fwrite(PC_summary_table, summary_output_file)

# Log total time taken for Protein complex datasets
PC_end_time <- Sys.time()
PC_time_taken <- round(
  as.numeric(difftime(PC_end_time, PC_start_time, units = "mins")),
  2
)
cat(
  "Finished generating Cannot-Link matrices for Protein complex datasets in",
  PC_time_taken,
  "minutes\n"
)
cat(paste0(rep("#", 80), collapse = ""), "\n")

# Log total time taken for entire process
end_time <- Sys.time()
total_time_taken <- round(
  as.numeric(difftime(end_time, start_time, units = "mins")),
  2
)
cat(
  "Finished generating Cannot-Link matrices for all datasets in",
  total_time_taken,
  "minutes\n"
)
cat(paste0(rep("#", 80), collapse = ""), "\n")
