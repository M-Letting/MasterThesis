library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(scico)

setwd(here::here())

# ── Constants ──────────────────────────────────────────────────────────────────

TARGET_METHOD <- "VSClustDCF_Corr"
TARGET_THRESHOLD <- 0.5
DATA_DIR <- "Q-Benchmark/03-Clustering/00-Data/05-CompareClustering"
OUT_DIR <- "00-ThesisFigures/Figures/Clustering"

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ── Palette ────────────────────────────────────────────────────────────────────

na_levels <- c("LowNA", "MedNA", "HighNA")
na_labels <- c("Low NA", "Med NA", "High NA")
na_pal <- scico(3, palette = "batlow")
na_colors <- setNames(na_pal, na_labels)

x_levels <- c("Baseline\n(Fixed)", "Fixed", "Baseline\n(1/k)", "1/k")

# ── Loaders ────────────────────────────────────────────────────────────────────
#
# Each file is loaded for both threshold modes. Baseline rows (BlindDCF) and
# constrained rows (TARGET_METHOD) are kept. NALevel is extracted from Dataset
# via gsub(".*_", ""), so both "LowNA" and "04_LowNA" become "LowNA".
# Constrained rows are averaged across simulation runs within (NALevel, Accession).
# XGroup places each row at one of the four x-positions; CompGroup identifies
# which comparison pair a row belongs to ("Fixed" or "1/k").

load_pm <- function(stem, value_col, dedup_by = NULL, log_base = NULL) {
  rbindlist(lapply(c("fixed", "one_over_k"), function(mode) {
    dt <- fread(file.path(DATA_DIR, paste0(stem, "_", mode, ".csv")))

    dt[, NALevel := gsub(".*_", "", Dataset)]
    dt <- dt[NALevel %in% na_levels]
    dt <- dt[
      (Method == "BlindDCF" & is.na(ConstraintThreshold)) |
        (Method == TARGET_METHOD & ConstraintThreshold == TARGET_THRESHOLD)
    ]

    dt[, v := as.numeric(get(value_col))]
    if (!is.null(dedup_by)) {
      dt <- unique(dt, by = dedup_by)
    }
    if (!is.null(log_base)) {
      dt[v <= 0, v := NA_real_]
      dt[, v := log(v, base = log_base)]
    }

    agg <- dt[,
      .(v = mean(v, na.rm = TRUE)),
      by = .(NALevel, Accession, Method, ConstraintThreshold)
    ]

    ml <- if (mode == "fixed") "Fixed" else "1/k"
    agg[, CompGroup := ml]
    agg[,
      XGroup := fifelse(
        Method == "BlindDCF",
        paste0("Baseline\n(", ml, ")"),
        ml
      )
    ]
    agg[, NALevel := factor(NALevel, levels = na_levels, labels = na_labels)]
    agg
  }))
}

load_cx <- function(stem, value_col, dedup_by = NULL, log10_tf = FALSE) {
  rbindlist(lapply(c("fixed", "one_over_k"), function(mode) {
    dt <- fread(file.path(DATA_DIR, paste0(stem, "_", mode, ".csv")))

    dt <- dt[
      (Method == "BlindDCF" & is.na(ConstraintThreshold)) |
        (Method == TARGET_METHOD & ConstraintThreshold == TARGET_THRESHOLD)
    ]

    dt[, v := as.numeric(get(value_col))]
    if (!is.null(dedup_by)) {
      dt <- unique(dt, by = dedup_by)
    }
    if (log10_tf) {
      dt[v <= 0, v := NA_real_]
      dt[, v := log10(v)]
    }

    agg <- dt[,
      .(v = mean(v, na.rm = TRUE)),
      by = .(Complex, Method, ConstraintThreshold)
    ]

    ml <- if (mode == "fixed") "Fixed" else "1/k"
    agg[, CompGroup := ml]
    agg[,
      XGroup := fifelse(
        Method == "BlindDCF",
        paste0("Baseline\n(", ml, ")"),
        ml
      )
    ]
    agg
  }))
}

# ── Load data ──────────────────────────────────────────────────────────────────

pm_assigned <- load_pm("ProteoMakerAssignedEvaluation", "Assigned_Features")
pm_fraction <- load_pm("ProteoMakerFractionEvaluation", "Top_Fraction")
pm_purity <- load_pm(
  "ProteoMakerClusterPurityEvaluation",
  "Weighted_Mean_Cluster_Purity"
)
pm_sumuniq <- load_pm(
  "ProteoMakerSumUniqueProteoforms",
  "Sum_Unique_Proteoforms"
)
pm_totpot <- load_pm(
  "ProteoMakerTotalPotentialProteoforms",
  "Sum_Potential_Proteoforms",
  log_base = 2
)
pm_sse <- load_pm(
  "ProteoMakerWeightedSSEEvaluation",
  "Total_Weighted_SSE",
  dedup_by = c("Dataset", "Accession", "Method", "ConstraintThreshold")
)

cx_fraction <- load_cx("ComplexFractionEvaluation", "Top_Fraction")
cx_assigned <- load_cx("ComplexAssignedEvaluation", "Assigned_Features")
cx_position <- load_cx(
  "ComplexPositionPotentialProteoforms",
  "Sum_Position_Proteoforms"
)
cx_sse <- load_cx(
  "ComplexWeightedSSEEvaluation",
  "Total_Weighted_SSE",
  dedup_by = c("Complex", "Method", "ConstraintThreshold")
)

# Remove Dynein-1 complex from CX panels
cx_fraction <- cx_fraction[Complex != "Dynein-1-complex-variant-1"]
cx_assigned <- cx_assigned[Complex != "Dynein-1-complex-variant-1"]
cx_position <- cx_position[Complex != "Dynein-1-complex-variant-1"]
cx_sse <- cx_sse[Complex != "Dynein-1-complex-variant-1"]

# ── Shared theme ───────────────────────────────────────────────────────────────

panel_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.position = "none"
  )

pair_divider <- geom_vline(
  xintercept = 2.5,
  linetype = "dotted",
  color = "grey50",
  linewidth = 0.5
)

# Shared stat_summary layers — one mean line and one median line per CompGroup
mean_median_layers <- function() {
  list(
    stat_summary(
      aes(x = factor(XGroup, levels = x_levels), y = v, group = CompGroup),
      fun = mean,
      geom = "line",
      color = "black",
      linewidth = 1.1,
      inherit.aes = FALSE
    ),
    stat_summary(
      aes(x = factor(XGroup, levels = x_levels), y = v, group = CompGroup),
      fun = mean,
      geom = "point",
      color = "black",
      shape = 18,
      size = 5,
      inherit.aes = FALSE
    ),
    stat_summary(
      aes(x = factor(XGroup, levels = x_levels), y = v, group = CompGroup),
      fun = median,
      geom = "line",
      color = "red3",
      linewidth = 1.1,
      linetype = "dashed",
      inherit.aes = FALSE
    ),
    stat_summary(
      aes(x = factor(XGroup, levels = x_levels), y = v, group = CompGroup),
      fun = median,
      geom = "point",
      color = "red3",
      shape = 18,
      size = 5,
      inherit.aes = FALSE
    )
  )
}

# ── PM panel: beeswarm coloured by NA level ────────────────────────────────────

make_pm_panel <- function(dt, y_label, title) {
  ggplot(
    dt,
    aes(x = factor(XGroup, levels = x_levels), y = v, fill = NALevel)
  ) +
    pair_divider +
    geom_beeswarm(
      size = 1.8,
      shape = 21,
      alpha = 0.6,
      stroke = 0.15,
      method = "compactswarm",
      corral = "gutter",
      corral.width = 0.55
    ) +
    mean_median_layers() +
    scale_fill_manual(values = na_colors, name = "ProteoMaker Dataset") +
    scale_x_discrete(limits = x_levels, expand = expansion(add = 0.5)) +
    labs(title = title, y = y_label) +
    panel_theme
}

# ── CX panel: beeswarm coloured by complex (batlow gradient) ──────────────────

make_cx_panel <- function(dt, y_label, title) {
  dt <- copy(dt)
  cx_order <- sort(unique(dt$Complex))
  dt[, ComplexIdx := as.numeric(factor(Complex, levels = cx_order))]

  ggplot(
    dt,
    aes(x = factor(XGroup, levels = x_levels), y = v, fill = ComplexIdx)
  ) +
    pair_divider +
    geom_beeswarm(
      size = 2.5,
      shape = 21,
      alpha = 0.7,
      stroke = 0.2,
      method = "compactswarm",
      corral = "gutter",
      corral.width = 0.55
    ) +
    mean_median_layers() +
    scale_fill_scico(palette = "batlow", guide = "none") +
    scale_x_discrete(limits = x_levels, expand = expansion(add = 0.5)) +
    labs(title = title, y = y_label) +
    panel_theme
}

# ── Build panels ───────────────────────────────────────────────────────────────

p_pm_assign <- make_pm_panel(
  pm_assigned,
  "Count",
  "Number of Assigned Features\n"
)
p_pm_frac <- make_pm_panel(
  pm_fraction,
  "Fraction",
  "Fraction of Top Proteoform\nin Top Cluster"
)
p_pm_purity <- make_pm_panel(
  pm_purity,
  "Weighted Mean Cluster Purity",
  "Cluster Purity\n"
)
p_pm_sumuniq <- make_pm_panel(
  pm_sumuniq,
  "Count (lower is better)",
  "Sum of Unique Proteoform IDs\n"
)
p_pm_totpot <- make_pm_panel(
  pm_totpot,
  expression(log[2] * "(Count) (lower is better)"),
  "Total Potential Proteoforms\n"
)
p_pm_sse <- make_pm_panel(
  pm_sse,
  "Total Weighted SSE (lower is better)",
  "Weighted SSE\n"
)

p_cx_assign <- make_cx_panel(
  cx_assigned,
  "Count",
  "Number of Assigned Features"
)
p_cx_frac <- make_cx_panel(
  cx_fraction,
  "Fraction",
  "Fraction of NM-Peptide in Top Cluster"
)
p_cx_pos <- make_cx_panel(
  cx_position,
  "Count (lower is better)",
  "Position Based Potential Proteoforms"
)
p_cx_sse <- make_cx_panel(
  cx_sse,
  "Total Weighted SSE (lower is better)",
  "Weighted SSE"
)

# ── Legends ────────────────────────────────────────────────────────────────────

# NA-level legend (between rows)
legend_na <- get_legend(
  make_pm_panel(pm_fraction, "", "") +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        title = "ProteoMaker Dataset",
        override.aes = list(size = 3.5, shape = 21, alpha = 1, stroke = 0.5)
      )
    )
)

# Complex colour-bar legend (bottom)
legend_cx <- get_legend(
  make_cx_panel(cx_assigned, "", "") +
    scale_fill_scico(
      palette = "batlow",
      name = "Protein Complex (A-Z)",
      guide = guide_colorbar(
        barwidth = 8,
        barheight = 0.5,
        title.position = "top",
        title.hjust = 0.5,
        label = FALSE,
        ticks = FALSE
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9)
    )
)

# Summary (mean / median) legend (bottom)
summary_dummy <- data.frame(
  x = c(1, 2, 1, 2),
  y = c(1, 1, 2, 2),
  Summary = c("Mean", "Mean", "Median", "Median")
)
legend_summary <- get_legend(
  ggplot(
    summary_dummy,
    aes(x = x, y = y, color = Summary, linetype = Summary)
  ) +
    geom_line() +
    geom_point(shape = 18, size = 4) +
    scale_color_manual(
      values = c("Mean" = "black", "Median" = "red3"),
      name = "Summary"
    ) +
    scale_linetype_manual(
      values = c("Mean" = "solid", "Median" = "dashed"),
      name = "Summary"
    ) +
    guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.6, "cm")
    )
)

# ── Vertical row labels ────────────────────────────────────────────────────────

label_pm <- ggdraw() +
  draw_label("ProteoMaker", angle = 90, fontface = "bold", size = 12)
label_cx <- ggdraw() +
  draw_label("Protein Complex", angle = 90, fontface = "bold", size = 12)

# ── Assemble ───────────────────────────────────────────────────────────────────

pm_panels <- plot_grid(
  p_pm_assign,
  p_pm_frac,
  p_pm_purity,
  p_pm_sumuniq,
  p_pm_totpot,
  p_pm_sse,
  nrow = 1,
  label_size = 12
)

cx_panels <- plot_grid(
  p_cx_assign,
  p_cx_frac,
  p_cx_pos,
  p_cx_sse,
  nrow = 1,
  label_size = 12
)

# Row A (PM) and B (CX): label sits on the combined row via the outer plot_grid
pm_row <- plot_grid(label_pm, pm_panels, ncol = 2, rel_widths = c(0.03, 1))
cx_row <- plot_grid(label_cx, cx_panels, ncol = 2, rel_widths = c(0.03, 1))

bottom_legends <- plot_grid(
  legend_cx,
  legend_summary,
  nrow = 1,
  rel_widths = c(1.4, 0.8)
)

fig_title <- ggdraw() +
  draw_label(
    sprintf(
      "VSClust With vs Without DCF Constraints (%s | threshold %.2f)",
      TARGET_METHOD,
      TARGET_THRESHOLD
    ),
    fontface = "bold",
    size = 12,
    hjust = 0.5
  )

final_fig <- plot_grid(
  pm_row,
  legend_na,
  cx_row,
  bottom_legends,
  ncol = 1,
  rel_heights = c(0.8, 0.07, 1, 0.09),
  labels = c("A", "", "B", ""),
  label_size = 14
)

print(final_fig)

# ── Save ───────────────────────────────────────────────────────────────────────

ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints.pdf"),
  final_fig,
  width = 16,
  height = 8
)
ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints.png"),
  final_fig,
  width = 16,
  height = 8,
  dpi = 300
)

cat("Saved:", file.path(OUT_DIR, "ClusteringWithConstraints.pdf/.png\n"))

# ── Tall version: PM 2×3, CX 2×2 ──────────────────────────────────────────────

pm_panels_tall <- plot_grid(
  p_pm_assign,
  p_pm_frac,
  p_pm_purity,
  p_pm_sumuniq,
  p_pm_totpot,
  p_pm_sse,
  ncol = 3,
  nrow = 2
)

cx_panels_tall <- plot_grid(
  p_cx_assign,
  p_cx_frac,
  p_cx_pos,
  p_cx_sse,
  ncol = 2,
  nrow = 2
)

pm_row_tall <- plot_grid(
  label_pm,
  pm_panels_tall,
  ncol = 2,
  rel_widths = c(0.03, 1)
)
cx_row_tall <- plot_grid(
  label_cx,
  cx_panels_tall,
  ncol = 2,
  rel_widths = c(0.03, 1)
)

final_fig_tall <- plot_grid(
  pm_row_tall,
  legend_na,
  cx_row_tall,
  bottom_legends,
  ncol = 1,
  rel_heights = c(1, 0.06, 0.72, 0.08),
  labels = c("A", "", "B", ""),
  label_size = 14
)

ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_tall.pdf"),
  final_fig_tall,
  width = 10,
  height = 12
)
ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_tall.png"),
  final_fig_tall,
  width = 10,
  height = 12,
  dpi = 300
)

cat("Saved:", file.path(OUT_DIR, "ClusteringWithConstraints_tall.pdf/.png\n"))

# ── Summary table of means ─────────────────────────────────────────────────────

summarise_metric <- function(
  dt,
  metric_name,
  data_type,
  lower_is_better = FALSE
) {
  s <- dt[, .(Mean = mean(v, na.rm = TRUE)), by = .(CompGroup, Method)]
  s[, Method := fifelse(Method == "BlindDCF", "Baseline", "Constrained")]
  s[, Metric := metric_name]
  s[, DataType := data_type]
  wide <- dcast(s, DataType + Metric + CompGroup ~ Method, value.var = "Mean")
  setnames(wide, "CompGroup", "Mode")
  wide[, Delta := round(Constrained - Baseline, 4)]
  wide[, Baseline := round(Baseline, 4)]
  wide[, Constrained := round(Constrained, 4)]
  wide[, Lower_is_better := lower_is_better]
  wide
}

means_table <- rbindlist(list(
  summarise_metric(pm_assigned, "Assigned Features", "ProteoMaker"),
  summarise_metric(pm_fraction, "Fraction in Top Cluster", "ProteoMaker"),
  summarise_metric(pm_purity, "Cluster Purity", "ProteoMaker"),
  summarise_metric(
    pm_sumuniq,
    "Sum Unique Proteoforms",
    "ProteoMaker",
    lower_is_better = TRUE
  ),
  summarise_metric(
    pm_totpot,
    "Total Potential Proteoforms (log2)",
    "ProteoMaker",
    lower_is_better = TRUE
  ),
  summarise_metric(
    pm_sse,
    "Weighted SSE",
    "ProteoMaker",
    lower_is_better = TRUE
  ),
  summarise_metric(cx_assigned, "Assigned Features", "Protein Complex"),
  summarise_metric(
    cx_fraction,
    "Fraction NM-Peptide in Top Cluster",
    "Protein Complex"
  ),
  summarise_metric(
    cx_position,
    "Position Based Proteoforms",
    "Protein Complex",
    lower_is_better = TRUE
  ),
  summarise_metric(
    cx_sse,
    "Weighted SSE",
    "Protein Complex",
    lower_is_better = TRUE
  )
))

setcolorder(
  means_table,
  c(
    "DataType",
    "Metric",
    "Mode",
    "Baseline",
    "Constrained",
    "Delta",
    "Lower_is_better"
  )
)
setorder(means_table, DataType, Metric, Mode)

print(means_table, digits = 4)

fwrite(means_table, file.path(OUT_DIR, "ClusteringMeans.csv"))
cat("Saved:", file.path(OUT_DIR, "ClusteringMeans.csv\n"))

# ── Summary-only figures (mean + median lines, no beeswarm) ────────────────────
#
# PM: one line per NA level (color), solid = mean / dashed = median.
# CX: single black mean line + single red dashed median line across all complexes.

make_pm_summary_panel <- function(dt, y_label, title) {
  ggplot(
    dt,
    aes(
      x = factor(XGroup, levels = x_levels),
      y = v,
      group = CompGroup
    )
  ) +
    pair_divider +
    stat_summary(fun = mean, geom = "line", color = "black", linewidth = 1.1) +
    stat_summary(
      fun = mean,
      geom = "point",
      color = "black",
      shape = 18,
      size = 4.5
    ) +
    stat_summary(
      fun = median,
      geom = "line",
      color = "red3",
      linewidth = 1.1,
      linetype = "dashed"
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      color = "red3",
      shape = 17,
      size = 3.5
    ) +
    scale_x_discrete(limits = x_levels, expand = expansion(add = 0.5)) +
    labs(title = title, y = y_label) +
    panel_theme
}

make_cx_summary_panel <- function(dt, y_label, title) {
  ggplot(
    dt,
    aes(
      x = factor(XGroup, levels = x_levels),
      y = v,
      group = CompGroup
    )
  ) +
    pair_divider +
    stat_summary(fun = mean, geom = "line", color = "black", linewidth = 1.1) +
    stat_summary(
      fun = mean,
      geom = "point",
      color = "black",
      shape = 18,
      size = 4.5
    ) +
    stat_summary(
      fun = median,
      geom = "line",
      color = "red3",
      linewidth = 1.1,
      linetype = "dashed"
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      color = "red3",
      shape = 17,
      size = 3.5
    ) +
    scale_x_discrete(limits = x_levels, expand = expansion(add = 0.5)) +
    labs(title = title, y = y_label) +
    panel_theme
}

# Build summary panels (reuse same titles / y-labels as full figure)
sp_pm_assign <- make_pm_summary_panel(
  pm_assigned,
  "Count",
  "Number of Assigned Features\n"
)
sp_pm_frac <- make_pm_summary_panel(
  pm_fraction,
  "Fraction",
  "Fraction of Top Proteoform\nin Top Cluster"
)
sp_pm_purity <- make_pm_summary_panel(
  pm_purity,
  "Weighted Mean Cluster Purity",
  "Cluster Purity\n"
)
sp_pm_sumuniq <- make_pm_summary_panel(
  pm_sumuniq,
  "Count (lower is better)",
  "Sum of Unique Proteoform IDs\n"
)
sp_pm_totpot <- make_pm_summary_panel(
  pm_totpot,
  expression(log[2] * "(Count, lower is better)"),
  "Total Potential Proteoforms\n"
)
sp_pm_sse <- make_pm_summary_panel(
  pm_sse,
  "Total Weighted SSE (lower is better)",
  "Weighted SSE\n"
)

sp_cx_assign <- make_cx_summary_panel(
  cx_assigned,
  "Count",
  "Number of Assigned Features"
)
sp_cx_frac <- make_cx_summary_panel(
  cx_fraction,
  "Fraction",
  "Fraction of NM-Peptide in Top Cluster"
)
sp_cx_pos <- make_cx_summary_panel(
  cx_position,
  "Count (lower is better)",
  "Position Based Potential Proteoforms"
)
sp_cx_sse <- make_cx_summary_panel(
  cx_sse,
  "Total Weighted SSE (lower is better)",
  "Weighted SSE"
)

# NA-level colour legend (uses colour, not fill, so extract fresh)
legend_na_sum <- get_legend(
  make_pm_summary_panel(pm_fraction, "", "") +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm")
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        override.aes = list(shape = 15, size = 4, alpha = 1)
      )
    )
)

# Mean / median line legend
legend_sum_lines <- get_legend(
  ggplot(
    data.frame(
      x = c(1, 2, 1, 2),
      y = c(1, 1, 2, 2),
      Summary = c("Mean", "Mean", "Median", "Median")
    ),
    aes(x = x, y = y, color = Summary, linetype = Summary, shape = Summary)
  ) +
    geom_line() +
    geom_point(size = 4) +
    scale_color_manual(
      values = c("Mean" = "black", "Median" = "red3"),
      name = "Summary"
    ) +
    scale_linetype_manual(
      values = c("Mean" = "solid", "Median" = "dashed"),
      name = "Summary"
    ) +
    scale_shape_manual(
      values = c("Mean" = 18, "Median" = 17),
      name = "Summary"
    ) +
    guides(
      color = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1),
      shape = guide_legend(nrow = 1)
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.6, "cm")
    )
)

# ── Wide summary figure ────────────────────────────────────────────────────────

sp_pm_row <- plot_grid(
  label_pm,
  plot_grid(
    sp_pm_assign,
    sp_pm_frac,
    sp_pm_purity,
    sp_pm_sumuniq,
    sp_pm_totpot,
    sp_pm_sse,
    nrow = 1
  ),
  ncol = 2,
  rel_widths = c(0.03, 1)
)
sp_cx_row <- plot_grid(
  label_cx,
  plot_grid(sp_cx_assign, sp_cx_frac, sp_cx_pos, sp_cx_sse, nrow = 1),
  ncol = 2,
  rel_widths = c(0.03, 1)
)

final_sum <- plot_grid(
  sp_pm_row,
  sp_cx_row,
  legend_sum_lines,
  ncol = 1,
  rel_heights = c(0.8, 1, 0.09),
  labels = c("A", "B", ""),
  label_size = 14
)

ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_summary.pdf"),
  final_sum,
  width = 16,
  height = 8
)
ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_summary.png"),
  final_sum,
  width = 16,
  height = 8,
  dpi = 300
)

# ── Tall summary figure ────────────────────────────────────────────────────────

sp_pm_row_tall <- plot_grid(
  label_pm,
  plot_grid(
    sp_pm_assign,
    sp_pm_frac,
    sp_pm_purity,
    sp_pm_sumuniq,
    sp_pm_totpot,
    sp_pm_sse,
    nrow = 1
  ),
  ncol = 2,
  rel_widths = c(0.03, 1)
)
sp_cx_row_tall <- plot_grid(
  label_cx,
  plot_grid(sp_cx_assign, sp_cx_frac, sp_cx_pos, sp_cx_sse, nrow = 1),
  ncol = 2,
  rel_widths = c(0.03, 1)
)

final_sum_tall <- plot_grid(
  sp_pm_row_tall,
  sp_cx_row_tall,
  legend_sum_lines,
  ncol = 1,
  rel_heights = c(1, 0.72, 0.08),
  labels = c("A", "B", ""),
  label_size = 14
)

ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_summary_tall.pdf"),
  final_sum_tall,
  width = 10,
  height = 14
)
ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_summary_tall.png"),
  final_sum_tall,
  width = 10,
  height = 14,
  dpi = 300
)

cat(
  "Saved: ClusteringWithConstraints_summary(.pdf/.png) and _summary_tall(.pdf/.png)\n"
)

# ── Dynamite (bar + error) figures ─────────────────────────────────────────────
#
# Two variants:
#   _dynamite_mean   : bar = mean, error = ±1 SEM
#   _dynamite_median : bar = median, error = IQR (Q25–Q75)
#
# PM panels: dodged bars by NA level (same colour palette as beeswarm figures).
# CX panels: single bar per x-group, filled by Baseline vs Constrained.

# Custom stat helper: median + IQR
median_iqr <- function(x) {
  x <- x[!is.na(x)]
  data.frame(
    y = median(x),
    ymin = quantile(x, 0.25),
    ymax = quantile(x, 0.75)
  )
}

cx_method_colors <- c("Baseline" = "#7BAFD4", "Constrained" = "#E07B54")

prep_cx_method <- function(dt) {
  dt <- copy(dt)
  dt[, MethodType := fifelse(Method == "BlindDCF", "Baseline", "Constrained")]
  dt
}

# ── Mean dynamite panel builders ───────────────────────────────────────────────

make_pm_dynamite_mean <- function(dt, y_label, title) {
  ggplot(
    dt,
    aes(x = factor(XGroup, levels = x_levels), y = v, fill = NALevel)
  ) +
    pair_divider +
    stat_summary(
      fun = mean,
      geom = "bar",
      position = position_dodge(width = 0.8),
      width = 0.7,
      alpha = 0.85,
      color = "black",
      linewidth = 0.3
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      position = position_dodge(width = 0.8),
      width = 0.25,
      linewidth = 0.6
    ) +
    scale_fill_manual(values = na_colors, name = "ProteoMaker Dataset") +
    scale_x_discrete(limits = x_levels, expand = expansion(add = 0.5)) +
    labs(title = title, y = y_label) +
    panel_theme
}

make_cx_dynamite_mean <- function(dt, y_label, title) {
  dt <- prep_cx_method(dt)
  ggplot(
    dt,
    aes(x = factor(XGroup, levels = x_levels), y = v, fill = MethodType)
  ) +
    pair_divider +
    stat_summary(
      fun = mean,
      geom = "bar",
      width = 0.65,
      alpha = 0.85,
      color = "black",
      linewidth = 0.3
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.25,
      linewidth = 0.6
    ) +
    scale_fill_manual(
      values = cx_method_colors,
      name = "Method",
      guide = "none"
    ) +
    scale_x_discrete(limits = x_levels, expand = expansion(add = 0.5)) +
    labs(title = title, y = y_label) +
    panel_theme
}

# ── Median dynamite panel builders ─────────────────────────────────────────────

make_pm_dynamite_median <- function(dt, y_label, title) {
  ggplot(
    dt,
    aes(x = factor(XGroup, levels = x_levels), y = v, fill = NALevel)
  ) +
    pair_divider +
    stat_summary(
      fun = median,
      geom = "bar",
      position = position_dodge(width = 0.8),
      width = 0.7,
      alpha = 0.85,
      color = "black",
      linewidth = 0.3
    ) +
    stat_summary(
      fun.data = median_iqr,
      geom = "errorbar",
      position = position_dodge(width = 0.8),
      width = 0.25,
      linewidth = 0.6
    ) +
    scale_fill_manual(values = na_colors, name = "ProteoMaker Dataset") +
    scale_x_discrete(limits = x_levels, expand = expansion(add = 0.5)) +
    labs(title = title, y = y_label) +
    panel_theme
}

make_cx_dynamite_median <- function(dt, y_label, title) {
  dt <- prep_cx_method(dt)
  ggplot(
    dt,
    aes(x = factor(XGroup, levels = x_levels), y = v, fill = MethodType)
  ) +
    pair_divider +
    stat_summary(
      fun = median,
      geom = "bar",
      width = 0.65,
      alpha = 0.85,
      color = "black",
      linewidth = 0.3
    ) +
    stat_summary(
      fun.data = median_iqr,
      geom = "errorbar",
      width = 0.25,
      linewidth = 0.6
    ) +
    scale_fill_manual(
      values = cx_method_colors,
      name = "Method",
      guide = "none"
    ) +
    scale_x_discrete(limits = x_levels, expand = expansion(add = 0.5)) +
    labs(title = title, y = y_label) +
    panel_theme
}

# ── Build mean dynamite panels ─────────────────────────────────────────────────

dp_mean_pm_assign <- make_pm_dynamite_mean(
  pm_assigned,
  "Count",
  "Number of Assigned Features\n"
)
dp_mean_pm_frac <- make_pm_dynamite_mean(
  pm_fraction,
  "Fraction",
  "Fraction of Top Proteoform\nin Top Cluster"
)
dp_mean_pm_purity <- make_pm_dynamite_mean(
  pm_purity,
  "Weighted Mean Cluster Purity",
  "Cluster Purity\n"
)
dp_mean_pm_sumuniq <- make_pm_dynamite_mean(
  pm_sumuniq,
  "Count (lower is better)",
  "Sum of Unique Proteoform IDs\n"
)
dp_mean_pm_totpot <- make_pm_dynamite_mean(
  pm_totpot,
  expression(log[2] * "(Count) (lower is better)"),
  "Total Potential Proteoforms\n"
)
dp_mean_pm_sse <- make_pm_dynamite_mean(
  pm_sse,
  "Total Weighted SSE (lower is better)",
  "Weighted SSE\n"
)

dp_mean_cx_assign <- make_cx_dynamite_mean(
  cx_assigned,
  "Count",
  "Number of Assigned Features"
)
dp_mean_cx_frac <- make_cx_dynamite_mean(
  cx_fraction,
  "Fraction",
  "Fraction of NM-Peptide in Top Cluster"
)
dp_mean_cx_pos <- make_cx_dynamite_mean(
  cx_position,
  "Count (lower is better)",
  "Position Based Potential Proteoforms"
)
dp_mean_cx_sse <- make_cx_dynamite_mean(
  cx_sse,
  "Total Weighted SSE (lower is better)",
  "Weighted SSE"
)

# ── Build median dynamite panels ───────────────────────────────────────────────

dp_med_pm_assign <- make_pm_dynamite_median(
  pm_assigned,
  "Count",
  "Number of Assigned Features\n"
)
dp_med_pm_frac <- make_pm_dynamite_median(
  pm_fraction,
  "Fraction",
  "Fraction of Top Proteoform\nin Top Cluster"
)
dp_med_pm_purity <- make_pm_dynamite_median(
  pm_purity,
  "Weighted Mean Cluster Purity",
  "Cluster Purity\n"
)
dp_med_pm_sumuniq <- make_pm_dynamite_median(
  pm_sumuniq,
  "Count (lower is better)",
  "Sum of Unique Proteoform IDs\n"
)
dp_med_pm_totpot <- make_pm_dynamite_median(
  pm_totpot,
  expression(log[2] * "(Count) (lower is better)"),
  "Total Potential Proteoforms\n"
)
dp_med_pm_sse <- make_pm_dynamite_median(
  pm_sse,
  "Total Weighted SSE (lower is better)",
  "Weighted SSE\n"
)

dp_med_cx_assign <- make_cx_dynamite_median(
  cx_assigned,
  "Count",
  "Number of Assigned Features"
)
dp_med_cx_frac <- make_cx_dynamite_median(
  cx_fraction,
  "Fraction",
  "Fraction of NM-Peptide in Top Cluster"
)
dp_med_cx_pos <- make_cx_dynamite_median(
  cx_position,
  "Count (lower is better)",
  "Position Based Potential Proteoforms"
)
dp_med_cx_sse <- make_cx_dynamite_median(
  cx_sse,
  "Total Weighted SSE (lower is better)",
  "Weighted SSE"
)

# ── Legends for dynamite figures ───────────────────────────────────────────────

# NA-level legend (PM, shared between mean and median variants)
legend_na_dp <- get_legend(
  make_pm_dynamite_mean(pm_fraction, "", "") +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        title = "ProteoMaker Dataset",
        override.aes = list(alpha = 1)
      )
    )
)

# CX method fill legend (Baseline vs Constrained)
legend_cx_method <- get_legend(
  ggplot(
    data.frame(
      x = c(1, 2),
      y = c(1, 2),
      MethodType = c("Baseline", "Constrained")
    ),
    aes(x = x, y = y, fill = MethodType)
  ) +
    geom_col() +
    scale_fill_manual(values = cx_method_colors, name = "Protein Complex") +
    guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    )
)

# Error-bar description legend (one for SEM, one for IQR)
make_errorbar_legend <- function(label) {
  get_legend(
    ggplot(
      data.frame(x = 1, y = 1, g = label),
      aes(x = x, y = y, fill = g)
    ) +
      geom_col(alpha = 0.85) +
      geom_errorbar(aes(ymin = 0.7, ymax = 1.3), width = 0.25) +
      scale_fill_manual(
        values = setNames("grey50", label),
        name = "Error bars"
      ) +
      guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm")
      )
  )
}

legend_sem <- make_errorbar_legend("± SEM")
legend_iqr <- make_errorbar_legend("IQR (Q25–Q75)")

# ── Assemble mean dynamite figure (wide) ───────────────────────────────────────

assemble_dynamite <- function(
  pm_panels,
  cx_panels,
  legend_na,
  legend_cx,
  legend_err
) {
  pm_row <- plot_grid(
    label_pm,
    plot_grid(plotlist = pm_panels, nrow = 1),
    ncol = 2,
    rel_widths = c(0.03, 1)
  )
  cx_row <- plot_grid(
    label_cx,
    plot_grid(plotlist = cx_panels, nrow = 1),
    ncol = 2,
    rel_widths = c(0.03, 1)
  )
  bottom <- plot_grid(
    legend_cx,
    legend_err,
    nrow = 1,
    rel_widths = c(1, 0.9)
  )
  plot_grid(
    pm_row,
    legend_na,
    cx_row,
    bottom,
    ncol = 1,
    rel_heights = c(0.8, 0.07, 1, 0.09),
    labels = c("A", "", "B", ""),
    label_size = 14
  )
}

assemble_dynamite_tall <- function(
  pm_panels,
  cx_panels,
  legend_na,
  legend_cx,
  legend_err
) {
  pm_row <- plot_grid(
    label_pm,
    plot_grid(plotlist = pm_panels, ncol = 3, nrow = 2),
    ncol = 2,
    rel_widths = c(0.03, 1)
  )
  cx_row <- plot_grid(
    label_cx,
    plot_grid(plotlist = cx_panels, ncol = 2, nrow = 2),
    ncol = 2,
    rel_widths = c(0.03, 1)
  )
  bottom <- plot_grid(
    legend_cx,
    legend_err,
    nrow = 1,
    rel_widths = c(1, 0.9)
  )
  plot_grid(
    pm_row,
    legend_na,
    cx_row,
    bottom,
    ncol = 1,
    rel_heights = c(1, 0.06, 0.72, 0.08),
    labels = c("A", "", "B", ""),
    label_size = 14
  )
}

# ── Save mean dynamite ─────────────────────────────────────────────────────────

fig_dp_mean <- assemble_dynamite(
  pm_panels = list(
    dp_mean_pm_assign,
    dp_mean_pm_frac,
    dp_mean_pm_purity,
    dp_mean_pm_sumuniq,
    dp_mean_pm_totpot,
    dp_mean_pm_sse
  ),
  cx_panels = list(
    dp_mean_cx_assign,
    dp_mean_cx_frac,
    dp_mean_cx_pos,
    dp_mean_cx_sse
  ),
  legend_na = legend_na_dp,
  legend_cx = legend_cx_method,
  legend_err = legend_sem
)

ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_dynamite_mean.pdf"),
  fig_dp_mean,
  width = 16,
  height = 8
)
ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_dynamite_mean.png"),
  fig_dp_mean,
  width = 16,
  height = 8,
  dpi = 300
)

fig_dp_mean_tall <- assemble_dynamite_tall(
  pm_panels = list(
    dp_mean_pm_assign,
    dp_mean_pm_frac,
    dp_mean_pm_purity,
    dp_mean_pm_sumuniq,
    dp_mean_pm_totpot,
    dp_mean_pm_sse
  ),
  cx_panels = list(
    dp_mean_cx_assign,
    dp_mean_cx_frac,
    dp_mean_cx_pos,
    dp_mean_cx_sse
  ),
  legend_na = legend_na_dp,
  legend_cx = legend_cx_method,
  legend_err = legend_sem
)

ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_dynamite_mean_tall.pdf"),
  fig_dp_mean_tall,
  width = 10,
  height = 12
)
ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_dynamite_mean_tall.png"),
  fig_dp_mean_tall,
  width = 10,
  height = 12,
  dpi = 300
)

cat(
  "Saved: ClusteringWithConstraints_dynamite_mean(.pdf/.png) and _tall(.pdf/.png)\n"
)

# ── Save median dynamite ───────────────────────────────────────────────────────

fig_dp_med <- assemble_dynamite(
  pm_panels = list(
    dp_med_pm_assign,
    dp_med_pm_frac,
    dp_med_pm_purity,
    dp_med_pm_sumuniq,
    dp_med_pm_totpot,
    dp_med_pm_sse
  ),
  cx_panels = list(
    dp_med_cx_assign,
    dp_med_cx_frac,
    dp_med_cx_pos,
    dp_med_cx_sse
  ),
  legend_na = legend_na_dp,
  legend_cx = legend_cx_method,
  legend_err = legend_iqr
)

ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_dynamite_median.pdf"),
  fig_dp_med,
  width = 16,
  height = 8
)
ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_dynamite_median.png"),
  fig_dp_med,
  width = 16,
  height = 8,
  dpi = 300
)

fig_dp_med_tall <- assemble_dynamite_tall(
  pm_panels = list(
    dp_med_pm_assign,
    dp_med_pm_frac,
    dp_med_pm_purity,
    dp_med_pm_sumuniq,
    dp_med_pm_totpot,
    dp_med_pm_sse
  ),
  cx_panels = list(
    dp_med_cx_assign,
    dp_med_cx_frac,
    dp_med_cx_pos,
    dp_med_cx_sse
  ),
  legend_na = legend_na_dp,
  legend_cx = legend_cx_method,
  legend_err = legend_iqr
)

ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_dynamite_median_tall.pdf"),
  fig_dp_med_tall,
  width = 10,
  height = 12
)
ggsave(
  file.path(OUT_DIR, "ClusteringWithConstraints_dynamite_median_tall.png"),
  fig_dp_med_tall,
  width = 10,
  height = 12,
  dpi = 300
)

cat(
  "Saved: ClusteringWithConstraints_dynamite_median(.pdf/.png) and _tall(.pdf/.png)\n"
)
