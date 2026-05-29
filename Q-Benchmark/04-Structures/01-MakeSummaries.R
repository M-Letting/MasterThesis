# Libraries
library(data.table)
library(arrow)

# Set working directory to project root for sourcing
setwd(here::here())

# Source functions for creating CF summaries
source("ComplexoFinder/04-Clustering/ClusterSummary.R")

# Set output directory for summaries
output_dir <- "Q-Benchmark/04-Structures/00-Data/Summaries"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

################################################################################
## CF summaries ################################################################
################################################################################

# Define input
CF_input <- "Q-Benchmark/03-Clustering/00-Data/04-VSClustConstraints/01-PC_DCF_VSClust_ClusteringResults_corr.rds"

# Load data
CF_data <- readRDS(CF_input)

# Select the relevant subset of data
CF_data <- CF_data$`05`

# Subset to relevant complexes
CF_data <- CF_data[c(
  "26S-proteasome",
  "BSSADCRC",
  "Dynein-1-complex-variant-1",
  "Major-Spliceosomal-B",
  "Eukaryotic-translation-initiation-factor-3-complex"
)]

# Rename BSSADCRC to Brain-SWI-SNF for consistency with structure mapping
names(CF_data)[names(CF_data) == "BSSADCRC"] <- "Brain-SWI-SNF"

# Create summaries
CF_summaries <- list()

for (i in seq_along(CF_data)) {
  name <- names(CF_data)[i]
  membership_matrix <- CF_data[[i]]$vsclust_result$ClustOut$Bestcl$membership
  summary <- make_cluster_summary(
    membership_matrix,
    id_cols = c("Complex ID", "Gene name", "Accession", "Position", "Datatype")
  )

  # Remove Datatype = "NM" if present
  summary <- summary[Datatype != "NM"]

  CF_summaries[[name]] <- summary
}

# Save CF summaries
for (name in names(CF_summaries)) {
  summary <- CF_summaries[[name]]
  output_path <- file.path(
    output_dir,
    paste0(name, "_ComplexoFinderSummary.csv")
  )
  fwrite(summary, output_path, sep = ",", row.names = FALSE)
}

################################################################################
## PF summaries ################################################################
################################################################################

# Define input
# PF_input <- list(
#   "26S-proteasome" = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_26S-proteasome_result.feather",
#   "Brain-SWI-SNF" = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_BSSADCRC_result.feather",
#   "Dynein-1-complex-variant-1" = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Dynein-1-complex-variant-1_result.feather",
#   "Eukaryotic-translation-initiation-factor-3-complex" = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Eukaryotic-translation-initiation-factor-3-complex_result.feather",
#   "Major-Spliceosomal-B" = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Major-Spliceosomal-B_result.feather"
# )

# # Load data
# PF_data <- lapply(PF_input, read_feather)
# PF_data <- lapply(PF_data, as.data.table)

# # Subset to relevant columns
# PF_data <- lapply(PF_data, function(dt) dt[, .(peptide_id, ClusterID)])

# # Keep only 1 row per peptide_id and ClusterID combination
# PF_data <- lapply(PF_data, function(dt) {
#   dt[!duplicated(paste(peptide_id, ClusterID))]
# })

# # Create summaries
# PF_summaries <- list()

# for (name in names(PF_data)) {
#   # Split peptide_id into "Complex ID", "Gene name", "Accession", "Position", "Datatype"
#   dt <- PF_data[[name]]
#   dt[,
#     c(
#       "Gene name",
#       "Accession",
#       "Position",
#       "Datatype"
#     ) := tstrsplit(
#       peptide_id,
#       "_",
#       fixed = TRUE
#     )
#   ]
#   dt[, peptide_id := NULL]

#   # Remove Datatype = "NM" if present
#   dt <- dt[Datatype != "NM"]

#   PF_summaries[[name]] <- dt
# }

# # Save PF summaries
# for (name in names(PF_summaries)) {
#   summary <- PF_summaries[[name]]
#   output_path <- file.path(output_dir, paste0(name, "_ProteoForgeSummary.csv"))
#   fwrite(summary, output_path, sep = ",", row.names = FALSE)
# }
