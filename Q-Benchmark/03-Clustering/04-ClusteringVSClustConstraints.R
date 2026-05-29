# Libraries
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Set working directory to project root for sourcing functions
setwd(here::here())

# Source functions
source("ComplexoFinder/03-Constraints/Constraints.R")
source("ComplexoFinder/04-Clustering/VsClustWrappers.R")
source("ComplexoFinder/04-Clustering/VSClust.R")

# Set working directory to Q-Benchmark
setwd(here::here("Q-Benchmark"))

# Set output directory for results
output_dir <- "03-Clustering/00-Data/04-VSClustConstraints"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cannot_link_variant <- getOption("qbenchmark.cannot_link_variant", "corr")
if (!cannot_link_variant %in% c("corr", "rmsd", "ccc")) {
  stop("Unsupported cannot-link variant: ", cannot_link_variant)
}

include_nona <- isTRUE(getOption("qbenchmark.include_nona", TRUE))

################################################################################
## Run VSClust on datasets from 02-DCF with VSClust constraints ################
################################################################################

start_time <- Sys.time()

##########################
## ProteoMaker datasets ##
##########################

PM_start_time <- Sys.time()

PM_summary <- data.table(
  CannotLinkMethod = character(),
  Method = character(),
  Threshold = numeric(),
  Dataset = character(),
  Protein = character(),
  nPeptides = integer(),
  TimeTaken = numeric()
)

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

constraints_thresholds <- c("04", "05", "06", "07", "08")

###############
# RPC results #
###############

# PM_RPC_start_time <- Sys.time()

# cat(paste("#", 80), "\n")
# cat(
#   "Running VSClust with VSClust constraints on ProteoMaker datasets (RPC results, variant: ",
#   cannot_link_variant,
#   ")\n"
# )

# # Define dataset paths
# datasets <- list(
#   "LowNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerLowNA_results.csv",
#   "LowNoNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerLowNoNA_results.csv",
#   "MedNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerMedNA_results.csv",
#   "MedNoNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerMedNoNA_results.csv",
#   "HighNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerHighNA_results.csv",
#   "HighNoNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerHighNoNA_results.csv"
# )

# # Get accessions for top 25 proteins with most peptides
# high_no_na_data <- fread(datasets[["HighNoNA"]])
# top_proteins_high_no_na <- high_no_na_data[,
#   .N,
#   by = Accession
# ][order(-N)][1:25, Accession]

# # Define different constraint thresholds
# constraints_thresholds <- c("04", "05", "06", "07", "08")

# PM_RPC <- list()

# # Run VSClust with DCF constraints for each protein, dataset, and threshold
# for (th in constraints_thresholds) {
#   cat(paste(rep("=", 60), collapse = ""), "\n")
#   cat(paste0(
#     "Running VSClust with DCF constraints (threshold ",
#     th,
#     ") on ProteoMaker datasets\n"
#   ))

#   # Loop through datasets
#   for (name in names(datasets)) {
#     cat(paste(rep("-", 40), collapse = ""), "\n")
#     cat(paste0("Processing dataset: ", name, "\n"))

#     # Load dataset
#     filepath <- datasets[[name]]
#     if (!file.exists(filepath)) {
#       cat("  Dataset file not found:", filepath, "\n")
#       next
#     }
#     data <- fread(filepath)

#     # Get intensity columns
#     intensity_cols <- grep("^C\\d+_R\\d+$", colnames(data), value = TRUE)

#     # Loop through proteins
#     for (protein in top_proteins_high_no_na) {
#       cat(paste(rep("-", 20), collapse = ""), "\n")
#       cat(paste0("Processing protein: ", protein, "\n"))

#       clustering_start_time <- Sys.time()

#       # Subset data for current protein
#       protein_data <- data[Accession == protein]

#       if (nrow(protein_data) == 0) {
#         cat("No data found for protein:", protein, "!\n")
#         next
#       }

#       # Get cannot-link matrix for current protein, dataset, and threshold
#       CLM_file <- paste0(
#         "03-Clustering/00-Data/02-CannotLinkMatrices/PM-",
#         name,
#         "_",
#         protein,
#         "_CLM",
#         th,
#         "_",
#         cannot_link_variant,
#         ".csv"
#       )
#       if (!file.exists(CLM_file)) {
#         cat("  Cannot-link matrix file not found:", CLM_file, "\n")
#         next
#       }
#       CLM <- as.matrix(read.csv(CLM_file, row.names = 1))

#       # Align CLM identifiers to the Identifier scheme used by
#       # dCF_to_vsclust_constraints() with id_cols = c("Accession", "Peptide_ID"):
#       # - protein_data Identifier:   Accession_Peptide_ID
#       # - CLM rownames (from create_cannot_link): Peptide_Accession_<rest>
#       #   where Peptide_ID == Peptide_<rest>
#       expected_ids <- protein_data[, paste(Accession, Peptide_ID, sep = "_")]

#       old_ids <- rownames(CLM)
#       peptide_id_from_clm <- sub(
#         paste0("^([^_]+)_", protein, "_"),
#         "\\1_",
#         old_ids,
#         perl = TRUE
#       )
#       new_ids <- paste(protein, peptide_id_from_clm, sep = "_")

#       if (!is.null(colnames(CLM))) {
#         colnames(CLM) <- new_ids
#       }
#       rownames(CLM) <- new_ids

#       common_ids <- intersect(expected_ids, rownames(CLM))
#       if (length(common_ids) == 0) {
#         cat(
#           "  Identifier mismatch: no overlap between protein_data and CLM for protein ",
#           protein,
#           "\n",
#           sep = ""
#         )
#         next
#       }

#       CLM <- CLM[common_ids, common_ids, drop = FALSE]

#       # Convert to clustering constraints
#       constraint_result <- vsclust_to_restrictions(
#         data = protein_data,
#         cannotlink_matrix = CLM,
#         n_clusters = find_nclust(protein_data),
#         n_rep = 3,
#         n_cond = 10,
#         id_cols = c("Accession", "Peptide_ID"),
#         value_cols = intensity_cols,
#         verbose = FALSE
#       )
#       constraint_matrix <- constraint_result$restriction_matrix

#       # Find number of clusters
#       nClust <- ncol(constraint_matrix)
#       if (nClust == 0) {
#         cat("  No clusters available after applying constraints. Skipping.\n")
#         next
#       }

#       # Prepare protein_data for VSClust
#       if (!"Identifier" %in% names(protein_data)) {
#         protein_data[, Identifier := paste(Accession, Peptide_ID, sep = "_")]
#       }
#       vsclust_data <- protein_data[,
#         c("Identifier", intensity_cols),
#         with = FALSE
#       ]

#       # Run VSClust with constraints
#       vsclust_result <- vsclust_on_complex(
#         data = vsclust_data,
#         id_cols = "Identifier",
#         intensity_cols = intensity_cols,
#         n_rep = 3,
#         n_cond = 10,
#         n_clusters = nClust,
#         grouped_replicates = TRUE,
#         restriction_matrix = constraint_matrix,
#         sds_amplification_factor = 1
#       )

#       # Store results in list
#       if (is.null(PM_RPC[[th]])) {
#         PM_RPC[[th]] <- list()
#       }
#       if (is.null(PM_RPC[[th]][[name]])) {
#         PM_RPC[[th]][[name]] <- list()
#       }
#       PM_RPC[[th]][[name]][[protein]] <- list(
#         vsclust_result = vsclust_result,
#         constraint_matrix = constraint_matrix
#       )

#       # Log summary info
#       nPeptides <- nrow(protein_data)
#       time_taken <- round(
#         as.numeric(difftime(Sys.time(), clustering_start_time, units = "secs")),
#         2
#       )
#       PM_summary <- rbind(
#         PM_summary,
#         data.table(
#           CannotLinkMethod = cannot_link_variant,
#           Method = "RPC DCF",
#           Threshold = th,
#           Dataset = name,
#           Protein = protein,
#           nPeptides = nPeptides,
#           TimeTaken = time_taken
#         )
#       )
#     }
#   }
# }

# # Save clustering results and constraints
# saveRDS(
#   PM_RPC,
#   file.path(
#     output_dir,
#     paste0(
#       "/01-PM_RPC_VSClust_ClusteringResults_",
#       cannot_link_variant,
#       ".rds"
#     )
#   )
# )

# # Log total time taken for ProteoMaker RPC datasets
# total_time <- round(
#   as.numeric(difftime(Sys.time(), PM_start_time, units = "mins")),
#   2
# )
# cat(paste(
#   "Total time taken for ProteoMaker RPC datasets:",
#   total_time,
#   "minutes\n"
# ))
# cat(paste(rep("#", 80), collapse = ""), "\n")

##################
# DCF results #
##################

PM_DCF_start_time <- Sys.time()

# Define dataset paths
if (include_nona) {
  datasets <- list(
    "LowNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerLowNA_results.csv",
    "LowNoNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerLowNoNA_results.csv",
    "MedNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerMedNA_results.csv",
    "MedNoNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerMedNoNA_results.csv",
    "HighNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerHighNA_results.csv",
    "HighNoNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerHighNoNA_results.csv"
  )
} else {
  datasets <- list(
    "LowNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerLowNA_results.csv",
    "MedNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerMedNA_results.csv",
    "HighNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerHighNA_results.csv"
  )
}

PM_DCF <- list()

# Run VSClust with DCF constraints for each protein, dataset, and threshold
for (th in constraints_thresholds) {
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat(paste0(
    "Running VSClust with DCF constraints (threshold ",
    th,
    ") on ProteoMaker datasets\n"
  ))

  # Loop through datasets
  for (name in names(datasets)) {
    cat(paste(rep("-", 40), collapse = ""), "\n")
    cat(paste0("Processing dataset: ", name, "\n"))

    # Load dataset
    filepath <- datasets[[name]]
    if (!file.exists(filepath)) {
      cat("  Dataset file not found:", filepath, "\n")
      next
    }
    data <- fread(filepath)

    # Get intensity columns
    intensity_cols <- grep("^C\\d+_R\\d+$", colnames(data), value = TRUE)

    # Loop through proteins
    for (protein in top_proteins_high_no_na) {
      cat(paste(rep("-", 20), collapse = ""), "\n")
      cat(paste0("Processing protein: ", protein, "\n"))

      clustering_start_time <- Sys.time()

      # Subset data for current protein
      protein_data <- data[Accession == protein]

      if (nrow(protein_data) == 0) {
        cat("No data found for protein:", protein, "!\n")
        next
      }

      # Get cannot-link matrix for current protein, dataset, and threshold
      CLM_file <- paste0(
        "03-Clustering/00-Data/02-CannotLinkMatrices/PM-",
        name,
        "_",
        protein,
        "_CLM",
        th,
        "_",
        cannot_link_variant,
        ".csv"
      )
      if (!file.exists(CLM_file)) {
        cat("  Cannot-link matrix file not found:", CLM_file, "\n")
        next
      }
      CLM <- as.matrix(read.csv(CLM_file, row.names = 1))

      # Align CLM identifiers to the Identifier scheme used by
      # dCF_to_vsclust_constraints() with id_cols = c("Accession", "Peptide_ID"):
      # - protein_data Identifier:   Accession_Peptide_ID
      # - CLM rownames (from create_cannot_link): Peptide_Accession_<rest>
      #   where Peptide_ID == Peptide_<rest>
      expected_ids <- protein_data[, paste(Accession, Peptide_ID, sep = "_")]

      old_ids <- rownames(CLM)
      peptide_id_from_clm <- sub(
        paste0("^([^_]+)_", protein, "_"),
        "\\1_",
        old_ids,
        perl = TRUE
      )
      new_ids <- paste(protein, peptide_id_from_clm, sep = "_")

      if (!is.null(colnames(CLM))) {
        colnames(CLM) <- new_ids
      }
      rownames(CLM) <- new_ids

      common_ids <- intersect(expected_ids, rownames(CLM))
      if (length(common_ids) == 0) {
        cat(
          "  Identifier mismatch: no overlap between protein_data and CLM for protein ",
          protein,
          "\n",
          sep = ""
        )
        next
      }

      CLM <- CLM[common_ids, common_ids, drop = FALSE]

      # Convert to clustering constraints
      constraint_result <- vsclust_to_restrictions(
        data = protein_data,
        cannotlink_matrix = CLM,
        n_clusters = find_nclust(protein_data),
        n_rep = 3,
        n_cond = 10,
        id_cols = c("Accession", "Peptide_ID"),
        value_cols = intensity_cols,
        verbose = FALSE
      )
      constraint_matrix <- constraint_result$restriction_matrix

      # Find number of clusters
      nClust <- ncol(constraint_matrix)
      if (nClust == 0) {
        cat("  No clusters available after applying constraints. Skipping.\n")
        next
      }

      # Prepare protein_data for VSClust
      if (!"Identifier" %in% names(protein_data)) {
        protein_data[, Identifier := paste(Accession, Peptide_ID, sep = "_")]
      }
      vsclust_data <- protein_data[,
        c("Identifier", intensity_cols),
        with = FALSE
      ]

      # Run VSClust with constraints
      vsclust_result <- vsclust_on_complex(
        data = vsclust_data,
        id_cols = "Identifier",
        intensity_cols = intensity_cols,
        n_rep = 3,
        n_cond = 10,
        n_clusters = nClust,
        grouped_replicates = TRUE,
        restriction_matrix = constraint_matrix,
        sds_amplification_factor = 1
      )

      # Store results in list
      if (is.null(PM_DCF[[th]])) {
        PM_DCF[[th]] <- list()
      }
      if (is.null(PM_DCF[[th]][[name]])) {
        PM_DCF[[th]][[name]] <- list()
      }
      PM_DCF[[th]][[name]][[protein]] <- list(
        vsclust_result = vsclust_result,
        constraint_matrix = constraint_matrix
      )

      # Log summary info
      nPeptides <- nrow(protein_data)
      time_taken <- round(
        as.numeric(difftime(Sys.time(), clustering_start_time, units = "secs")),
        2
      )
      PM_summary <- rbind(
        PM_summary,
        data.table(
          CannotLinkMethod = cannot_link_variant,
          Method = "RPC2-2 DCF",
          Threshold = th,
          Dataset = name,
          Protein = protein,
          nPeptides = nPeptides,
          TimeTaken = time_taken
        )
      )
    }
  }
}

# Save ProteoMaker RPC2-2 DCF results
saveRDS(
  PM_DCF,
  file.path(
    output_dir,
    paste0(
      "01-PM_DCF_VSClust_ClusteringResults_",
      cannot_link_variant,
      ".rds"
    )
  )
)

# Save summary table
summary_output_file <- file.path(
  output_dir,
  paste0("/01-PM_DCF_VSClust_Summary_", cannot_link_variant, ".csv")
)
fwrite(PM_summary, summary_output_file)

# Log total time taken for ProteoMaker RPC2-2 datasets
total_time <- round(
  as.numeric(difftime(Sys.time(), PM_DCF_start_time, units = "mins")),
  2
)
cat(paste(
  "Total time taken for ProteoMaker RPC2-2 datasets:",
  total_time,
  "minutes\n"
))
cat(paste(rep("#", 80), collapse = ""), "\n")

# Log total time taken for all ProteoMaker datasets
total_time <- round(
  as.numeric(difftime(Sys.time(), PM_start_time, units = "mins")),
  2
)
cat(paste(
  "Total time taken for all ProteoMaker datasets:",
  total_time,
  "minutes\n"
))
cat(paste(rep("#", 80), collapse = ""), "\n")

##############################
## Protein complex datasets ##
##############################

cat(paste("#", 80), "\n")
cat("Running VSClust with DCF constraints on Protein complex datasets\n")
cat(paste("#", 80), "\n")

PC_start_time <- Sys.time()

PC_summary <- data.table(
  CannotLinkMethod = character(),
  Method = character(),
  Threshold = numeric(),
  Dataset = character(),
  Protein = character(),
  nPeptides = integer(),
  TimeTaken = numeric()
)

###############
# RPC results #
###############

# cat(
#   "Running VSClust with VSClust constraints on Protein complex datasets (RPC results, variant: ",
#   cannot_link_variant,
#   ")\n"
# )

# # Define dataset paths
# datasets <- get_complex_rpc_datasets("RPC")

# PC_RPC <- list()

# # Run VSClust with DCF constraints for each dataset and threshold
# for (th in constraints_thresholds) {
#   cat(paste(rep("=", 60), collapse = ""), "\n")
#   cat(paste0(
#     "Running VSClust with DCF constraints (threshold ",
#     th,
#     ") on Protein complex datasets\n"
#   ))

#   # Loop through datasets
#   for (name in names(datasets)) {
#     cat(paste(rep("-", 40), collapse = ""), "\n")
#     cat(paste0("Processing dataset: ", name, "\n"))

#     complex_start_time <- Sys.time()

#     # Load dataset
#     filepath <- datasets[[name]]
#     if (!file.exists(filepath)) {
#       cat("  Dataset file not found:", filepath, "\n")
#       next
#     }
#     data <- fread(filepath)

#     # Get intensity columns
#     intensity_cols <- grep("^B\\d+_D\\d+$", colnames(data), value = TRUE)

#     # Get cannot-link matrix for current dataset and threshold
#     CLM_file <- paste0(
#       "03-Clustering/00-Data/02-CannotLinkMatrices/PC-",
#       name,
#       "_CLM",
#       th,
#       "_",
#       cannot_link_variant,
#       ".csv"
#     )
#     if (!file.exists(CLM_file)) {
#       cat("  Cannot-link matrix file not found:", CLM_file, "\n")
#       next
#     }
#     CLM_dt <- fread(CLM_file)
#     if (ncol(CLM_dt) < 2) {
#       cat("  Cannot-link matrix has invalid format:", CLM_file, "\n")
#       next
#     }

#     clm_row_ids <- CLM_dt[[1]]
#     CLM <- as.matrix(CLM_dt[, -1, with = FALSE])
#     storage.mode(CLM) <- "logical"
#     rownames(CLM) <- clm_row_ids
#     colnames(CLM) <- names(CLM_dt)[-1]

#     # Keep only square/overlapping identifier set
#     clm_ids <- intersect(rownames(CLM), colnames(CLM))
#     if (length(clm_ids) == 0) {
#       cat("  Cannot-link matrix has no overlapping row/col identifiers.\n")
#       next
#     }
#     CLM <- CLM[clm_ids, clm_ids, drop = FALSE]

#     # Align CLM IDs with IDs in current dataset
#     expected_ids <- data[, paste(complex_id, Peptide, sep = "_")]
#     common_ids <- intersect(expected_ids, rownames(CLM))
#     if (length(common_ids) == 0) {
#       cat(
#         "  Identifier mismatch: no overlap between dataset and CLM for ",
#         name,
#         "\n",
#         sep = ""
#       )
#       next
#     }
#     CLM <- CLM[common_ids, common_ids, drop = FALSE]

#     # Convert cannot-link matrix to clustering constraints
#     constraint_result <- vsclust_to_restrictions(
#       data = data,
#       cannotlink_matrix = CLM,
#       n_clusters = find_nclust(data),
#       n_rep = 3,
#       n_cond = 28,
#       id_cols = c("complex_id", "Peptide"),
#       value_cols = intensity_cols,
#       verbose = FALSE
#     )
#     constraint_matrix <- constraint_result$restriction_matrix

#     nClust <- ncol(constraint_matrix)
#     if (nClust == 0) {
#       cat("  No clusters available after applying constraints. Skipping.\n")
#       next
#     }

#     # Prepare data for VSClust
#     data[, Identifier := paste(complex_id, Peptide, sep = "_")]
#     vsclust_data <- data[, c("Identifier", intensity_cols), with = FALSE]

#     # Run VSClust with constraints
#     vsclust_result <- vsclust_on_complex(
#       data = vsclust_data,
#       id_cols = "Identifier",
#       intensity_cols = intensity_cols,
#       n_rep = 3,
#       n_cond = 28,
#       n_clusters = nClust,
#       grouped_replicates = TRUE,
#       restriction_matrix = constraint_matrix,
#       sds_amplification_factor = 1
#     )

#     # Store results in list
#     if (is.null(PC_RPC[[th]])) {
#       PC_RPC[[th]] <- list()
#     }
#     PC_RPC[[th]][[name]] <- list(
#       vsclust_result = vsclust_result,
#       constraint_matrix = constraint_matrix
#     )

#     # Log summary info
#     nPeptides <- nrow(data)
#     time_taken <- round(
#       as.numeric(difftime(Sys.time(), complex_start_time, units = "secs")),
#       2
#     )
#     PC_summary <- rbind(
#       PC_summary,
#       data.table(
#         CannotLinkMethod = cannot_link_variant,
#         Method = "RPC DCF",
#         Threshold = th,
#         Dataset = name,
#         Protein = "All",
#         nPeptides = nPeptides,
#         TimeTaken = time_taken
#       )
#     )
#   }
# }

# # Save complex RPC DCF results and summary
# saveRDS(
#   PC_RPC,
#   file.path(
#     output_dir,
#     paste0(
#       "01-PC_RPC_VSClust_ClusteringResults_",
#       cannot_link_variant,
#       ".rds"
#     )
#   )
# )

# # Log total time taken for Protein complex RPC datasets
# PC_RPC_time_taken <- round(
#   as.numeric(difftime(Sys.time(), PC_start_time, units = "mins")),
#   2
# )
# cat(
#   "Total time taken for Protein complex RPC datasets:",
#   PC_RPC_time_taken,
#   "minutes\n"
# )

###############
# DCF results #
###############

cat(
  "Running VSClust with VSClust constraints on Protein complex datasets (DCF results, variant: ",
  cannot_link_variant,
  ")\n"
)

DCF_start_time <- Sys.time()

# Define dataset paths
datasets <- get_complex_rpc_datasets("DCF")

PC_DCF <- list()

constraints_thresholds <- c("04", "05", "06", "07", "08")

# Run VSClust with DCF constraints for each dataset and threshold
for (th in constraints_thresholds) {
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat(paste0(
    "Running VSClust with DCF constraints (threshold ",
    th,
    ") on Protein complex datasets\n"
  ))

  # Loop through datasets
  for (name in names(datasets)) {
    cat(paste(rep("-", 40), collapse = ""), "\n")
    cat(paste0("Processing dataset: ", name, "\n"))

    complex_start_time <- Sys.time()

    # Load dataset
    filepath <- datasets[[name]]
    if (!file.exists(filepath)) {
      cat("  Dataset file not found:", filepath, "\n")
      next
    }
    data <- fread(filepath)

    # Get intensity columns
    intensity_cols <- grep("^B\\d+_D\\d+$", colnames(data), value = TRUE)

    # Get cannot-link matrix for current dataset and threshold
    CLM_file <- paste0(
      "03-Clustering/00-Data/02-CannotLinkMatrices/PC-",
      name,
      "_CLM",
      th,
      "_",
      cannot_link_variant,
      ".csv"
    )
    if (!file.exists(CLM_file)) {
      cat("  Cannot-link matrix file not found:", CLM_file, "\n")
      next
    }
    CLM_dt <- fread(CLM_file)
    if (ncol(CLM_dt) < 2) {
      cat("  Cannot-link matrix has invalid format:", CLM_file, "\n")
      next
    }

    clm_row_ids <- CLM_dt[[1]]
    CLM <- as.matrix(CLM_dt[, -1, with = FALSE])
    storage.mode(CLM) <- "logical"
    rownames(CLM) <- clm_row_ids
    colnames(CLM) <- names(CLM_dt)[-1]

    # Keep only square/overlapping identifier set
    clm_ids <- intersect(rownames(CLM), colnames(CLM))
    if (length(clm_ids) == 0) {
      cat("  Cannot-link matrix has no overlapping row/col identifiers.\n")
      next
    }
    CLM <- CLM[clm_ids, clm_ids, drop = FALSE]

    # Align CLM IDs with IDs in current dataset
    expected_ids <- data[, paste(complex_id, Peptide, sep = "_")]
    common_ids <- intersect(expected_ids, rownames(CLM))
    if (length(common_ids) == 0) {
      cat(
        "  Identifier mismatch: no overlap between dataset and CLM for ",
        name,
        "\n",
        sep = ""
      )
      next
    }
    CLM <- CLM[common_ids, common_ids, drop = FALSE]

    # Convert cannot-link matrix to clustering constraints
    constraint_result <- vsclust_to_restrictions(
      data = data,
      cannotlink_matrix = CLM,
      n_clusters = find_nclust(data),
      n_rep = 3,
      n_cond = 28,
      id_cols = c("complex_id", "Peptide"),
      value_cols = intensity_cols,
      verbose = FALSE
    )
    constraint_matrix <- constraint_result$restriction_matrix

    nClust <- ncol(constraint_matrix)
    if (nClust == 0) {
      cat("  No clusters available after applying constraints. Skipping.\n")
      next
    }

    # Prepare data for VSClust
    data[, Identifier := paste(complex_id, Peptide, sep = "_")]
    vsclust_data <- data[, c("Identifier", intensity_cols), with = FALSE]

    # Run VSClust with constraints
    vsclust_result <- vsclust_on_complex(
      data = vsclust_data,
      id_cols = "Identifier",
      intensity_cols = intensity_cols,
      n_rep = 3,
      n_cond = 28,
      n_clusters = nClust,
      grouped_replicates = TRUE,
      restriction_matrix = constraint_matrix,
      sds_amplification_factor = 1
    )

    # Store results in list
    if (is.null(PC_DCF[[th]])) {
      PC_DCF[[th]] <- list()
    }
    PC_DCF[[th]][[name]] <- list(
      vsclust_result = vsclust_result,
      constraint_matrix = constraint_matrix
    )

    # Log summary info
    nPeptides <- nrow(data)
    time_taken <- round(
      as.numeric(difftime(Sys.time(), complex_start_time, units = "secs")),
      2
    )
    PC_summary <- rbind(
      PC_summary,
      data.table(
        CannotLinkMethod = cannot_link_variant,
        Method = "DCF VSClust",
        Threshold = th,
        Dataset = name,
        Protein = "All",
        nPeptides = nPeptides,
        TimeTaken = time_taken
      )
    )
  }
}

# Save complex DCF results and summary
saveRDS(
  PC_DCF,
  file.path(
    output_dir,
    paste0(
      "01-PC_DCF_VSClust_ClusteringResults_",
      cannot_link_variant,
      ".rds"
    )
  )
)

# Log total time taken for Protein complex DCF datasets
PC_DCF_time_taken <- round(
  as.numeric(difftime(Sys.time(), DCF_start_time, units = "mins")),
  2
)
cat(
  "Total time taken for Protein complex DCF datasets:",
  PC_DCF_time_taken,
  "minutes\n"
)
cat(paste(rep("#", 80), collapse = ""), "\n")

# Save summary table
fwrite(
  PC_summary,
  file.path(
    output_dir,
    paste0("01-PC_DCF_VSClust_Summary_", cannot_link_variant, ".csv")
  )
)

# Log total time taken for all Protein complex datasets
total_time <- round(
  as.numeric(difftime(Sys.time(), PC_start_time, units = "mins")),
  2
)
cat(
  "Total time taken for all Protein complex datasets:",
  total_time,
  "minutes\n"
)
cat(paste(rep("#", 80), collapse = ""), "\n")

# Log total time taken for all datasets
total_time <- round(
  as.numeric(difftime(Sys.time(), start_time, units = "mins")),
  2
)
cat(
  "Total time taken for all datasets:",
  total_time,
  "minutes\n"
)
cat(paste(rep("#", 80), collapse = ""), "\n")
