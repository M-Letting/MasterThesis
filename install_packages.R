# Install all R packages used across scripts in this repository.
# Run this script once before executing any other scripts.

# --- CRAN packages ---
cran_packages <- c(
  "arrow",
  "cowplot",
  "data.table",
  "dendextend",
  "dplyr",
  "dynamicTreeCut",
  "ggbeeswarm",
  "ggplot2",
  "gt",
  "matrixStats",
  "patchwork",
  "readxl",
  "scico",
  "tidyr"
)

install.packages(cran_packages, repos = "https://cloud.r-project.org")

# --- Bioconductor packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c(
  "ComplexHeatmap",  # also requires circlize
  "circlize",
  "limma",
  "pcaMethods"
)

BiocManager::install(bioc_packages)

# --- GitHub packages ---
# vsclust: fuzzy c-means clustering for omics data
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_github("veitveit/vsclust")

# COPF
remotes::install_github("CCprofiler/CCprofiler", ref = "proteoformLocationMapping")