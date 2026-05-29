# Libraries
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(scico)
library(cowplot)

# Set working directory
setwd(paste0(here::here(), "/Q-Benchmark/03-Clustering"))

# ==============================================================================
# Target method — edit these two lines to switch methods
# ==============================================================================
target_method <- "VSClustDCF_Corr"
target_threshold <- 0.6 # cannot-link threshold
membership_threshold_mode <- "one_over_k" # "fixed" or "one_over_k"
axis_text <- "VSC Cor, 0.6" # for plot labels
# ==============================================================================

plot_output_dir <- "00-Plots"
compare_dir <- "00-Data/05-CompareClustering"

if (!dir.exists(plot_output_dir)) {
  dir.create(plot_output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Helper: load a comparison CSV, keep only baseline and target, add MethodLabel
# ------------------------------------------------------------------------------
load_metric <- function(
  filename,
  y_var,
  log10_transform = FALSE,
  exclude_nona = FALSE
) {
  path <- file.path(
    compare_dir,
    paste0(sub("\\.csv$", "", filename), "_", membership_threshold_mode, ".csv")
  )

  dt <- fread(path)
  dt <- dt[!grepl("RPC", Method)]
  if (exclude_nona) {
    dt <- dt[!grepl("NoNA", Dataset)]
  }
  dt[, Dataset := gsub(".*_", "", Dataset)]

  dt <- dt[
    (Method == "BlindDCF" & is.na(ConstraintThreshold)) |
      (Method == target_method & ConstraintThreshold == target_threshold)
  ]

  if (nrow(dt) == 0L) {
    stop(
      "No data found for method '",
      target_method,
      "' with threshold ",
      target_threshold,
      " in ",
      path
    )
  }

  dt[,
    MethodLabel := fifelse(Method == "BlindDCF", "Baseline", axis_text)
  ]
  dt[, (y_var) := as.numeric(get(y_var))]

  if (log10_transform) {
    dt[get(y_var) <= 0, (y_var) := NA_real_]
    dt[, (y_var) := log10(get(y_var))]
  }

  dt
}

# ------------------------------------------------------------------------------
# Helper: beeswarm plot — one point per observation, colored by Dataset,
# with mean (black diamond) and median (red diamond) overlaid
# ------------------------------------------------------------------------------
make_beeswarm <- function(dt, y_var, y_label, title, legend_position = "none") {
  dt <- dt[is.finite(get(y_var))]

  ggplot(dt, aes(x = MethodLabel, y = get(y_var), fill = Dataset)) +
    geom_beeswarm(
      size = 3,
      shape = 21,
      alpha = 0.7,
      method = "compactswarm",
      corral = "gutter",
      corral.width = 0.7
    ) +
    stat_summary(
      aes(
        x = MethodLabel,
        y = get(y_var),
        group = MethodLabel,
        color = "Median"
      ),
      fun = median,
      geom = "point",
      inherit.aes = FALSE,
      size = 2.5,
      shape = 18
    ) +
    stat_summary(
      aes(
        x = MethodLabel,
        y = get(y_var),
        group = MethodLabel,
        color = "Mean"
      ),
      fun = mean,
      geom = "point",
      inherit.aes = FALSE,
      size = 2.5,
      shape = 18
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = legend_position,
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank()
    ) +
    labs(
      title = title,
      y = y_label,
      fill = "Dataset",
      color = "Summary"
    ) +
    scale_fill_scico_d(palette = "batlow", na.value = "grey80") +
    scale_color_manual(values = c("Mean" = "black", "Median" = "red"))
}

# ------------------------------------------------------------------------------
# ProteoMaker panels
# ------------------------------------------------------------------------------

PM_Fraction <- make_beeswarm(
  load_metric(
    "ProteoMakerFractionEvaluation.csv",
    "Top_Fraction",
    exclude_nona = TRUE
  ),
  y_var = "Top_Fraction",
  y_label = "Fraction",
  title = "PM: Fraction of Proteoform ID 1 in Top Cluster"
)

PM_Purity <- make_beeswarm(
  load_metric(
    "ProteoMakerClusterPurityEvaluation.csv",
    "Weighted_Mean_Cluster_Purity",
    exclude_nona = TRUE
  ),
  y_var = "Weighted_Mean_Cluster_Purity",
  y_label = "WMCP",
  title = "PM: Weighted Mean Cluster Purity"
)

PM_TotalPotential <- make_beeswarm(
  load_metric(
    "ProteoMakerTotalPotentialProteoforms.csv",
    "Sum_Potential_Proteoforms",
    log10_transform = TRUE,
    exclude_nona = TRUE
  ),
  y_var = "Sum_Potential_Proteoforms",
  y_label = "Log10(Total)",
  title = "PM: Total Potential Proteoforms (log10)"
)

pm_sse_dt <- load_metric(
  "ProteoMakerWeightedSSEEvaluation.csv",
  "Total_Weighted_SSE",
  exclude_nona = TRUE
)
pm_sse_dt <- unique(pm_sse_dt, by = c("Dataset", "Accession", "MethodLabel"))
PM_TotalSSE <- make_beeswarm(
  pm_sse_dt,
  y_var = "Total_Weighted_SSE",
  y_label = "Total Weighted SSE",
  title = "PM: Total Weighted SSE"
)

# ------------------------------------------------------------------------------
# Complex panels
# ------------------------------------------------------------------------------

Complex_Fraction <- make_beeswarm(
  load_metric("ComplexFractionEvaluation.csv", "Top_Fraction"),
  y_var = "Top_Fraction",
  y_label = "Fraction",
  title = "Complex: Fraction of NMPeptide in Top Cluster"
)

cx_sse_dt <- load_metric(
  "ComplexWeightedSSEEvaluation.csv",
  "Total_Weighted_SSE"
)
cx_sse_dt <- unique(cx_sse_dt, by = c("Dataset", "Complex", "MethodLabel"))
Complex_TotalSSE <- make_beeswarm(
  cx_sse_dt,
  y_var = "Total_Weighted_SSE",
  y_label = "Total Weighted SSE",
  title = "Complex: Total Weighted SSE"
)

# ------------------------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------------------------

combined <- plot_grid(
  PM_Fraction,
  PM_Purity,
  PM_TotalPotential,
  PM_TotalSSE,
  Complex_Fraction,
  Complex_TotalSSE,
  ncol = 3,
  labels = "AUTO"
)

legend <- get_legend(
  PM_Fraction +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 1))
)

final_plot <- plot_grid(
  combined,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

method_safe <- gsub("[^A-Za-z0-9]", "_", target_method)
ggsave(
  filename = file.path(
    plot_output_dir,
    sprintf(
      "SingleMethod_%s_%.2f_%s.pdf",
      method_safe,
      target_threshold,
      membership_threshold_mode
    )
  ),
  plot = final_plot,
  width = 14,
  height = 8,
  dpi = 300
)
