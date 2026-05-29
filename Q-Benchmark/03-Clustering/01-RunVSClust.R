# Libraries
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Set working directory to project root for sourcing functions
setwd(here::here())

# Source functions
source("ComplexoFinder/03-Constraints/Constraints.R")
source("ComplexoFinder/04-Clustering/VsClustWrappers.R")
source("ComplexoFinder/04-Clustering/VSClust.R")

# Set working directory to Q-Benchmark
setwd(paste0(getwd(), "/Q-Benchmark"))

# Set output directory for results
output_dir <- "03-Clustering/00-Data/01-BlindVSClust"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

################################################################################
## Run VSClust on datasets from 02-DCF #########################################
################################################################################

start_time <- Sys.time()

summary <- data.table(
  Dataset = character(),
  Protein = character(),
  Method = character(),
  Nclust = integer(),
  NDiscordant = integer(),
  Time_Taken_Seconds = numeric()
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

######################################
## ProteoMaker datasets for VSClust ##
######################################

PM_start_time <- Sys.time()

###############
# DCF results #
###############

DCF_start_time <- Sys.time()

# Define paths to datasets for VSClust
datasets <- list(
  "LowNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerLowNA_results.csv",
  "LowNoNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerLowNoNA_results.csv",
  "MedNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerMedNA_results.csv",
  "MedNoNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerMedNoNA_results.csv",
  "HighNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerHighNA_results.csv",
  "HighNoNA" = "./02-DCF/00-Data/03-Results/DCF/DCF_ProteoMakerHighNoNA_results.csv"
)

# Get accessions for top 25 proteins with most peptides
high_no_na_data <- fread(datasets[["HighNoNA"]])
top_proteins_high_no_na <- high_no_na_data[,
  .N,
  by = Accession
][order(-N)][1:25, Accession]

PM_DCF <- list()

# Run VSClust on each dataset
for (name in names(datasets)) {
  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Running VSClust on ProteoMaker dataset:", name, "\n")

  dataset_start_time <- Sys.time()

  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  dataset <- fread(file_path)

  # Run VSClust on top 10 proteins with most peptides in HighNoNA dataset
  for (acc in top_proteins_high_no_na) {
    protein_start_time <- Sys.time()

    if (!acc %in% dataset$Accession) {
      warning(paste("Accession", acc, "not found in dataset", name))
      next
    }
    cat(paste0(rep("-", 40), collapse = ""), "\n")
    cat("Running VSClust on protein:", acc, "\n")

    # Subset to accession of interest
    data <- subset(dataset, Accession == acc)

    # Get intensity columns
    intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

    # Combine non-intensity columns into a single "Identifier" column
    data[, Identifier := paste(Accession, Peptide_ID, sep = "_")]
    data$Accession <- NULL
    data$Peptide_ID <- NULL

    # Get number of clusters based on DCF results
    nclust <- find_nclust(data)
    NDiscordant <- sum(data$discordant)
    data <- data[, c("Identifier", intensity_cols), with = FALSE]

    # Run VSClust
    result <- vsclust_on_complex(
      data = data,
      id_cols = "Identifier",
      intensity_cols = intensity_cols,
      n_rep = 3,
      n_cond = 10,
      n_clusters = nclust,
      grouped_replicates = TRUE,
      restriction_matrix = NULL,
      sds_amplification_factor = 1
    )

    # Store result in list
    PM_DCF[[name]][[acc]] <- result

    # Record time taken for this protein
    protein_end_time <- Sys.time()
    time_taken <- round(
      as.numeric(difftime(
        protein_end_time,
        protein_start_time,
        units = "secs"
      )),
      2
    )
    summary <- rbind(
      summary,
      data.table(
        Dataset = name,
        Protein = acc,
        Method = "DCF",
        Nclust = nclust,
        NDiscordant = NDiscordant,
        Time_Taken_Seconds = time_taken
      )
    )

    # fmt: skip
    cat(
      "Time taken for protein", acc, "in dataset", name, ":",
      time_taken, "seconds\n"
    )
  }

  # Record time taken for this dataset
  dataset_end_time <- Sys.time()
  time_taken <- round(
    as.numeric(difftime(
      dataset_end_time,
      dataset_start_time,
      units = "secs"
    )),
    2
  )
  cat("Time taken for dataset", name, ":", time_taken, "seconds\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")
}

# Save results
saveRDS(PM_DCF, file = paste0(output_dir, "/01-PM_DCF.rds"))

# Record total time taken for DCF datasets
DCF_end_time <- Sys.time()
DCF_time_taken <- round(
  as.numeric(difftime(
    DCF_end_time,
    DCF_start_time,
    units = "secs"
  )),
  2
)
cat("Total time taken for DCF datasets:", DCF_time_taken, "seconds\n")

###############
# RPC results #
###############

# RPC_start_time <- Sys.time()

# # Define paths to datasets for VSClust
# datasets <- list(
#   "LowNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerLowNA_results.csv",
#   "LowNoNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerLowNoNA_results.csv",
#   "MedNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerMedNA_results.csv",
#   "MedNoNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerMedNoNA_results.csv",
#   "HighNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerHighNA_results.csv",
#   "HighNoNA" = "./02-DCF/00-Data/03-Results/RPC/RPC_ProteoMakerHighNoNA_results.csv"
# )

# PM_RPC <- list()

# # Run VSClust on each dataset
# for (name in names(datasets)) {
#   cat(paste0(rep("#", 80), collapse = ""), "\n")
#   cat("Running VSClust on ProteoMaker dataset:", name, "\n")

#   dataset_start_time <- Sys.time()

#   file_path <- datasets[[name]]
#   if (!file.exists(file_path)) {
#     warning(paste("File not found:", file_path))
#     next
#   }

#   dataset <- fread(file_path)

#   # Run VSClust on top 10 proteins with most peptides in HighNoNA dataset
#   for (acc in top_proteins_high_no_na) {
#     protein_start_time <- Sys.time()

#     if (!acc %in% dataset$Accession) {
#       warning(paste("Accession", acc, "not found in dataset", name))
#       next
#     }
#     cat(paste0(rep("-", 40), collapse = ""), "\n")
#     cat("Running VSClust on protein:", acc, "\n")

#     # Subset to accession of interest
#     data <- subset(dataset, Accession == acc)

#     # Get intensity columns
#     intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

#     # Combine non-intensity columns into a single "Identifier" column
#     data[, Identifier := paste(Accession, Peptide_ID, sep = "_")]
#     data$Accession <- NULL
#     data$Peptide_ID <- NULL

#     # Get number of clusters based on RPC results
#     nclust <- find_nclust(data)
#     NDiscordant <- sum(data$discordant)
#     data <- data[, c("Identifier", intensity_cols), with = FALSE]

#     # Run VSClust
#     result <- vsclust_on_complex(
#       data = data,
#       id_cols = "Identifier",
#       intensity_cols = intensity_cols,
#       n_rep = 3,
#       n_cond = 10,
#       n_clusters = nclust,
#       grouped_replicates = TRUE,
#       restriction_matrix = NULL,
#       sds_amplification_factor = 1
#     )

#     # Store result in list
#     PM_RPC[[name]][[acc]] <- result

#     # Record time taken for this protein
#     protein_end_time <- Sys.time()
#     time_taken <- round(
#       as.numeric(difftime(
#         protein_end_time,
#         protein_start_time,
#         units = "secs"
#       )),
#       2
#     )
#     summary <- rbind(
#       summary,
#       data.table(
#         Dataset = name,
#         Protein = acc,
#         Method = "RPC",
#         Nclust = nclust,
#         NDiscordant = NDiscordant,
#         Time_Taken_Seconds = time_taken
#       )
#     )

#     # fmt: skip
#     cat(
#       "Time taken for protein", acc, "in dataset", name, ":",
#       time_taken, "seconds\n"
#     )
#   }

#   # Record time taken for this dataset
#   dataset_end_time <- Sys.time()
#   time_taken <- round(
#     as.numeric(difftime(
#       dataset_end_time,
#       dataset_start_time,
#       units = "secs"
#     )),
#     2
#   )
#   cat("Time taken for dataset", name, ":", time_taken, "seconds\n")
#   cat(paste0(rep("#", 80), collapse = ""), "\n")
# }

# # Save results
# saveRDS(PM_RPC, file = paste0(output_dir, "/01-PM_RPC.rds"))

# # Record total time taken for RPC datasets
# RPC_end_time <- Sys.time()
# RPC_time_taken <- round(
#   as.numeric(difftime(
#     RPC_end_time,
#     RPC_start_time,
#     units = "secs"
#   )),
#   2
# )
# cat("Total time taken for RPC datasets:", RPC_time_taken, "seconds\n\n")

# # Record total time taken for all ProteoMaker datasets
# PM_end_time <- Sys.time()
# total_time_taken <- round(
#   as.numeric(difftime(
#     PM_end_time,
#     PM_start_time,
#     units = "secs"
#   )),
#   2
# )
# cat(
#   "Total time taken for all ProteoMaker datasets:",
#   total_time_taken,
#   "seconds\n"
# )

##########################################
## Protein complex datasets for VSClust ##
##########################################

complex_start_time <- Sys.time()

###################
# RPC2-2 datasets #
###################

DCF_start_time <- Sys.time()

datasets <- get_complex_rpc_datasets("DCF")

complex_DCF <- list()

# Run VSClust on each dataset
for (name in names(datasets)) {
  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("Running VSClust on protein complex dataset:", name, "\n")

  dataset_start_time <- Sys.time()

  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  dataset <- fread(file_path)

  # Get intensity columns
  intensity_cols <- grep("^B\\d+_D\\d+$", names(dataset), value = TRUE)

  # Combine non-intensity columns into a single "Identifier" column
  dataset[, Identifier := paste(complex_id, Peptide, sep = "_")]
  dataset$complex_id <- NULL
  dataset$Peptide <- NULL

  # Get number of clusters based on RPC results
  nclust <- find_nclust(dataset)
  NDiscordant <- sum(dataset$discordant)
  dataset <- dataset[, c("Identifier", intensity_cols), with = FALSE]

  # Run VSClust
  result <- vsclust_on_complex(
    data = dataset,
    id_cols = "Identifier",
    intensity_cols = intensity_cols,
    n_rep = 3,
    n_cond = 28,
    n_clusters = nclust,
    grouped_replicates = TRUE,
    restriction_matrix = NULL,
    sds_amplification_factor = 1
  )

  # Store result in list
  complex_DCF[[name]] <- result

  # Record time taken for this dataset
  dataset_end_time <- Sys.time()
  time_taken <- round(
    as.numeric(difftime(
      dataset_end_time,
      dataset_start_time,
      units = "secs"
    )),
    2
  )
  cat("Time taken for dataset", name, ":", time_taken, "seconds\n")
  cat(paste0(rep("#", 80), collapse = ""), "\n")

  summary <- rbind(
    summary,
    data.table(
      Dataset = name,
      Protein = NA,
      Method = "DCF",
      Nclust = nclust,
      NDiscordant = NDiscordant,
      Time_Taken_Seconds = time_taken
    )
  )
}

# Save results
saveRDS(complex_DCF, file = paste0(output_dir, "/01-Complex_DCF.rds"))

# Record total time taken for DCF complex datasets
DCF_end_time <- Sys.time()
DCF_time_taken <- round(
  as.numeric(difftime(
    DCF_end_time,
    DCF_start_time,
    units = "secs"
  )),
  2
)
cat(
  "Total time taken for DCF complex datasets:",
  DCF_time_taken,
  "seconds\n"
)

################
# RPC datasets #
################

# RPC_start_time <- Sys.time()

# datasets <- get_complex_rpc_datasets("RPC")

# complex_RPC <- list()

# # Run VSClust on each dataset
# for (name in names(datasets)) {
#   cat(paste0(rep("#", 80), collapse = ""), "\n")
#   cat("Running VSClust on protein complex dataset:", name, "\n")

#   dataset_start_time <- Sys.time()

#   file_path <- datasets[[name]]
#   if (!file.exists(file_path)) {
#     warning(paste("File not found:", file_path))
#     next
#   }

#   dataset <- fread(file_path)

#   # Get intensity columns
#   intensity_cols <- grep("^B\\d+_D\\d+$", names(dataset), value = TRUE)

#   # Combine non-intensity columns into a single "Identifier" column
#   dataset[, Identifier := paste(complex_id, Peptide, sep = "_")]
#   dataset$complex_id <- NULL
#   dataset$Peptide <- NULL

#   # Get number of clusters based on RPC results
#   nclust <- find_nclust(dataset)
#   NDiscordant <- sum(dataset$discordant)
#   dataset <- dataset[, c("Identifier", intensity_cols), with = FALSE]

#   # Run VSClust
#   result <- vsclust_on_complex(
#     data = dataset,
#     id_cols = "Identifier",
#     intensity_cols = intensity_cols,
#     n_rep = 3,
#     n_cond = 28,
#     n_clusters = nclust,
#     grouped_replicates = TRUE,
#     restriction_matrix = NULL,
#     sds_amplification_factor = 1
#   )

#   # Store result in list
#   complex_RPC[[name]] <- result

#   # Record time taken for this dataset
#   dataset_end_time <- Sys.time()
#   time_taken <- round(
#     as.numeric(difftime(
#       dataset_end_time,
#       dataset_start_time,
#       units = "secs"
#     )),
#     2
#   )
#   cat("Time taken for dataset", name, ":", time_taken, "seconds\n")
#   cat(paste0(rep("#", 80), collapse = ""), "\n")

#   summary <- rbind(
#     summary,
#     data.table(
#       Dataset = name,
#       Protein = NA,
#       Method = "RPC",
#       Nclust = nclust,
#       NDiscordant = NDiscordant,
#       Time_Taken_Seconds = time_taken
#     )
#   )
# }

# # Save results
# saveRDS(complex_RPC, file = paste0(output_dir, "/01-Complex_RPC.rds"))

# # Record total time taken for RPC complex datasets
# RPC_end_time <- Sys.time()
# RPC_time_taken <- round(
#   as.numeric(difftime(
#     RPC_end_time,
#     RPC_start_time,
#     units = "secs"
#   )),
#   2
# )
# cat("Total time taken for RPC complex datasets:", RPC_time_taken, "seconds\n\n")

# Write sumamry table to CSV
fwrite(summary, file = paste0(output_dir, "/01-VSClust_Summary.csv"))

# Record total time taken for all VSClust analyses
end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
total_time_minutes <- floor(total_time / 60)
total_time_seconds <- round(total_time %% 60)
cat(
  "Total time taken for VSClust analysis:",
  total_time_minutes,
  "min",
  total_time_seconds,
  "seconds\n"
)
