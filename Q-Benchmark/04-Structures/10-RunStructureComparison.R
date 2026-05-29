# End-to-end runner for the structure benchmark pipeline.
#
# Prerequisites:
#   - ChimeraX installed at the path defined below (required by step 2)
#   - Upstream clustering results present (inputs to step 1)
#
# Steps:
#   1  Make cluster summary CSVs from CF and PF results
#   2  Run ChimeraX to assemble structures, write CXC scripts and distance files
#   3  Plot per-cluster pairwise distance heatmaps
#   4  Plot combined cross-cluster distance heatmaps (threshold and highest modes)

setwd(here::here())

# ChimeraX path – required by step 2. Adjust if your installation differs.
chimeraX_path <- "/Applications/ChimeraX-1.11.1.app/Contents/bin/ChimeraX"

step <- function(n, label) {
  message(
    "\n", strrep("=", 72), "\n",
    "Step ", n, ": ", label, "\n",
    strrep("=", 72)
  )
}

run_step <- function(n, label, script) {
  step(n, label)
  t0 <- proc.time()["elapsed"]
  source(script, local = FALSE)
  elapsed <- round(proc.time()["elapsed"] - t0, 1)
  message("Step ", n, " completed in ", elapsed, "s")
}

run_step(1, "Make cluster summary CSVs",
  "Q-Benchmark/04-Structures/01-MakeSummaries.R")

run_step(2, "Assemble structures and compute pairwise distances",
  "Q-Benchmark/04-Structures/02-BenchmarkStructures.R")

run_step(3, "Plot per-cluster distance heatmaps",
  "Q-Benchmark/04-Structures/03-PlotDistances.R")

run_step(4, "Plot combined cross-cluster distance heatmaps",
  "Q-Benchmark/04-Structures/04-PlotCombinedHeatmap.R")

message("\n", strrep("=", 72))
message("Structure benchmark complete.")
message("Figures written to: Q-Benchmark/04-Structures/00-Plots/")
message(strrep("=", 72), "\n")
