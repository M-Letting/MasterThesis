# End-to-end runner for the ComplexoVisualiserThesis pipeline.
#
# All clustering results are taken directly from ComplexoFinder/10-Results/,
# ensuring summaries match exactly what RunComplexoFinder.R produced.
#
# Prerequisites:
#   - ComplexoFinder has been run (ComplexoFinder/10-Results/ populated)
#   - ChimeraX installed at the path in 02-VisualiseStructures.R
#
# Steps:
#   1  Copy cluster summaries from ComplexoFinder results
#   2  Assemble structures and compute pairwise PTM distances
#   3  Plot per-cluster pairwise distance heatmaps
#   4  Plot combined cross-cluster distance heatmaps

setwd(here::here())

run_step <- function(n, label, script) {
  message("\n", strrep("=", 72), "\n", "Step ", n, ": ", label, "\n", strrep("=", 72))
  t0 <- proc.time()["elapsed"]
  source(script, local = FALSE)
  message("Step ", n, " completed in ", round(proc.time()["elapsed"] - t0, 1), "s")
}

run_step(1, "Make summaries from ComplexoFinder results",
  "ComplexoVisualiserThesis/01-MakeSummaries.R")

run_step(2, "Assemble structures and compute pairwise distances",
  "ComplexoVisualiserThesis/02-VisualiseStructures.R")

run_step(3, "Plot per-cluster distance heatmaps",
  "ComplexoVisualiserThesis/03-PlotDistances.R")

run_step(4, "Plot combined cross-cluster distance heatmaps",
  "ComplexoVisualiserThesis/04-PlotCombinedHeatmap.R")

message("\n", strrep("=", 72))
message("ComplexoVisualiserThesis pipeline complete.")
message("Figures written to: ComplexoVisualiserThesis/00-Plots/")
message(strrep("=", 72), "\n")
