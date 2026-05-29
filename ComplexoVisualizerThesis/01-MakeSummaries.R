# Libraries
library(data.table)

# Set working directory to project root
setwd(here::here())

# Output directory for summaries
output_dir <- "ComplexoVisualiserThesis/00-Data/Summaries"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

################################################################################
## Read cluster summaries directly from ComplexoFinder results #################
################################################################################

# Mapping: summary file name -> ComplexoFinder result folder name
complexes <- list(
  "26S-proteasome" =
    "ComplexoFinder/10-Results/26S_proteasome_complex/cluster_summary.csv",
  "Brain-SWI-SNF" = paste0(
    "ComplexoFinder/10-Results/",
    "Brain-specific_SWI_SNF_ATP-dependent_chromatin_remodeling_complex,_ARID1A-SMARCA2_variant",
    "/cluster_summary.csv"
  ),
  "Dynein-1-complex-variant-1" =
    "ComplexoFinder/10-Results/Dynein-1_complex,_variant_1/cluster_summary.csv",
  "Eukaryotic-translation-initiation-factor-3-complex" =
    "ComplexoFinder/10-Results/Eukaryotic_translation_initiation_factor_3_complex/cluster_summary.csv"
)

for (name in names(complexes)) {
  input_path <- complexes[[name]]

  if (!file.exists(input_path)) {
    warning("cluster_summary.csv not found for ", name, ": ", input_path)
    next
  }

  summary <- fread(input_path)

  # Drop Modification column if present (not used by downstream scripts)
  if ("Modification" %in% names(summary)) {
    summary[, Modification := NULL]
  }

  # Remove NM peptides (protein-level data)
  summary <- summary[!Datatype %in% c("NM", "NMpeptide")]

  output_path <- file.path(output_dir, paste0(name, "_ComplexoFinderSummary.csv"))
  fwrite(summary, output_path)
  cat("Saved:", output_path, "(", nrow(summary), "rows )\n")
}
