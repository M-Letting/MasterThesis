setwd(here::here())

cat(paste(rep("=", 80), collapse = ""), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
overall_start_time <- Sys.time()

# Run all clustering comparison scripts

# Run blind VSClust to get baseline results without constraints
source("Q-Benchmark/03-Clustering/01-RunVSClust.R", echo = FALSE)
setwd(here::here())

# Options for cannot-link constraint variants
cannot_link_variants <- c("corr", "ccc") #rmsd not run

# Option to include proteomaker datasets with no NA values in the comparison
include_nona <- FALSE # FALSE to exclude proteomaker datasets with no NA values

options(qbenchmark.include_nona = include_nona)

# Run clustering with different cannot-link constraint variants
for (variant in cannot_link_variants) {
  # Create cannot-link matrices based on the current variant
  options(qbenchmark.cannot_link_variant = variant)
  cat("Running cannot-link variant:", variant, "\n")
  source("Q-Benchmark/03-Clustering/02-CannotLink.R", echo = FALSE)

  setwd(here::here())

  # Run clustering with DCF-based constraints
  source(
    "Q-Benchmark/03-Clustering/03-ClusteringDCFConstraints.R",
    echo = FALSE
  )

  setwd(here::here())

  # Run clustering with VSClust constraints
  source(
    "Q-Benchmark/03-Clustering/04-ClusteringVSClustConstraints.R",
    echo = FALSE
  )

  setwd(here::here())
}
options(qbenchmark.cannot_link_variant = NULL)

setwd(here::here())

# Compare clustering results using different definitions of cluster membership
membership_threshold_modes <- c("fixed", "one_over_k")

for (mode in membership_threshold_modes) {
  options(qbenchmark.membership_threshold_mode = mode)
  options(qbenchmark.membership_threshold_value = 0.5)
  cat(
    "Running comparison with membership threshold mode:",
    mode,
    "\n"
  )
  source("Q-Benchmark/03-Clustering/05-CompareClustering.R", echo = FALSE)
  setwd(here::here())

  # Plot comparison results
  source("Q-Benchmark/03-Clustering/06-PlotComparison.R", echo = FALSE)

  setwd(here::here())
}

source("Q-Benchmark/03-Clustering/07-FindBestMethod.R", echo = FALSE)

options(qbenchmark.membership_threshold_mode = NULL)
options(qbenchmark.membership_threshold_value = NULL)
options(qbenchmark.include_nona = NULL)

overall_end_time <- Sys.time()
overall_time <- round(
  difftime(overall_end_time, overall_start_time, units = "hours"),
  2
)

cat("All clustering comparison scripts completed.\n")
cat("Total time taken for all comparisons:", overall_time, "hours\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
