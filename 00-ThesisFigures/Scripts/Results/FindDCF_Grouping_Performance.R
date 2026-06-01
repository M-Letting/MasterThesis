library(data.table)
library(ggplot2)
library(cowplot)
library(scico)

setwd(here::here())

# ── Data ──────────────────────────────────────────────────────────────────────

data <- fread(
  "Q-Benchmark/01-FindDCF/data/results/peptide_grouping_performance_data.csv"
)

data[,
  F1 := ifelse(
    Precision + Recall > 0,
    2 * Precision * Recall / (Precision + Recall),
    NA_real_
  )
]

# For each method × perturbation find the row with the highest F1 score
optimal <- data[
  !is.nan(F1) & !is.na(F1) & is.finite(F1),
  .SD[which.max(F1)],
  by = .(method, perturbation)
]

# ── Factor ordering ───────────────────────────────────────────────────────────

pert_levels <- c(
  "2 Peptides",
  "Random (2 to %50) Peptides",
  "%50 Peptides"
)
pert_labels <- c("Two", "Random", "50%")

method_levels <- c(
  "ComplexoFinder",
  "ComplexoFinderBaseline",
  "ProteoForge",
  "COPF"
)

optimal[,
  method := fcase(
    method == "DCF"         , "ComplexoFinder"         ,
    method == "DCFBaseline" , "ComplexoFinderBaseline" ,
    default = method
  )
]
optimal[,
  perturbation := factor(
    perturbation,
    levels = pert_levels,
    labels = pert_labels
  )
]
optimal[, method := factor(method, levels = method_levels)]

# ── Colours ───────────────────────────────────────────────────────────────────

all_method_levels <- c(
  "ComplexoFinder",
  "ComplexoFinderBaseline",
  "ProteoForge",
  "PeCorA",
  "COPF"
)
pal <- scico(5, palette = "batlow", direction = 1)
method_colors <- setNames(pal, all_method_levels)[method_levels]
method_faces <- ifelse(
  method_levels %in% c("ComplexoFinder", "ComplexoFinderBaseline"),
  "bold",
  "plain"
)

# ── Bar theme ─────────────────────────────────────────────────────────────────

bar_theme <- theme_classic(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10),
    legend.position = "none"
  )

# ── TPR bar chart ─────────────────────────────────────────────────────────────

p_tpr <- ggplot(optimal, aes(x = perturbation, y = TPR, fill = method)) +
  geom_col(
    position = position_dodge(width = 0.85),
    width = 0.75,
    alpha = 0.9,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.2f", TPR), color = method),
    position = position_dodge(width = 0.85),
    vjust = -0.35,
    size = 3.0,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1.2),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_color_manual(values = method_colors) +
  labs(title = "True Positive Rate", y = "TPR") +
  bar_theme

# ── FPR bar chart ─────────────────────────────────────────────────────────────

p_fpr <- ggplot(optimal, aes(x = perturbation, y = FPR, fill = method)) +
  geom_col(
    position = position_dodge(width = 0.85),
    width = 0.75,
    alpha = 0.9,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.2f", FPR), color = method),
    position = position_dodge(width = 0.85),
    vjust = -0.35,
    size = 3.0,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1.2),
    breaks = seq(0, 1, 0.25),
    labels = scales::number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_color_manual(values = method_colors) +
  labs(title = "False Positive Rate", y = "FPR") +
  bar_theme

# ── Shared legend (horizontal, bottom) ───────────────────────────────────────

legend <- get_legend(
  p_tpr +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        override.aes = list(face = method_faces)
      )
    )
)

# ── Threshold annotation table ────────────────────────────────────────────────

thresh_annot <- {
  thresh <- optimal[, .(method, perturbation, val = threshold)]
  thresh[, method := factor(method, levels = rev(method_levels))]

  ggplot(
    thresh,
    aes(
      x = perturbation,
      y = method,
      label = formatC(val, format = "e", digits = 1)
    )
  ) +
    geom_text(aes(color = method), fontface = "bold", size = 3.0) +
    scale_color_manual(values = method_colors) +
    scale_x_discrete(expand = expansion(add = 0.55)) +
    labs(title = "Optimal Thresholds") +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 9,
        color = "grey30",
        margin = margin(b = 4)
      ),
      axis.text.y = element_text(
        hjust = 1,
        size = 8,
        margin = margin(r = 4),
        color = method_colors[rev(method_levels)],
        face = rev(method_faces)
      ),
      legend.position = "none",
      plot.margin = margin(t = 4, r = 0, b = 4, l = 0)
    )
}

# ── Assemble with cowplot ─────────────────────────────────────────────────────

fig_title <- ggdraw() +
  draw_label(
    "Discordant Peptide Grouping",
    fontface = "bold",
    size = 14,
    hjust = 0.5
  )

combined <- plot_grid(
  fig_title,
  p_tpr,
  p_fpr,
  legend,
  thresh_annot,
  ncol = 1,
  rel_heights = c(0.1, 1, 1, 0.12, 0.4),
  labels = c("", "A", "B", "", "C"),
  label_size = 13
)

print(combined)

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave(
  "00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Performance.pdf",
  combined,
  width = 10,
  height = 6
)
ggsave(
  "00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Performance.png",
  combined,
  width = 10,
  height = 6,
  dpi = 300
)

cat(
  "Saved: 00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Performance.pdf/.png\n"
)

# ── Scatter: TPR vs FPR, one panel per perturbation level ────────────────────

scatter_data <- copy(data)
scatter_data[,
  method := fcase(
    method == "DCF"         , "ComplexoFinder"         ,
    method == "DCFBaseline" , "ComplexoFinderBaseline" ,
    default = method
  )
]
scatter_data[,
  perturbation := factor(
    perturbation,
    levels = pert_levels,
    labels = pert_labels
  )
]
scatter_data[, method := factor(method, levels = method_levels)]
scatter_data <- scatter_data[
  !is.nan(TPR) &
    !is.na(TPR) &
    !is.nan(FPR) &
    !is.na(FPR) &
    is.finite(TPR) &
    is.finite(FPR) &
    !is.na(perturbation)
]

# Highlight points at p = 0.001, 0.01, 0.05
pval_highlights <- c(0.001, 0.01, 0.05)

scatter_highlights <- rbindlist(lapply(pval_highlights, function(thr) {
  scatter_data[,
    .SD[which.min(abs(threshold - thr))],
    by = .(method, perturbation)
  ][, thr_label := paste0("p = ", thr)]
}))
scatter_highlights[,
  thr_label := factor(thr_label, levels = paste0("p = ", pval_highlights))
]

# Panel builder
make_scatter_panel <- function(lvl, title) {
  ggplot(
    scatter_data[perturbation == lvl],
    aes(x = FPR, y = TPR, colour = method)
  ) +
    geom_point(size = 1.5, alpha = 0.5, shape = 16) +
    geom_point(
      data = scatter_highlights[perturbation == lvl],
      aes(x = FPR, y = TPR, shape = thr_label, fill = method),
      size = 3.5,
      color = "black",
      stroke = 0.6
    ) +
    scale_colour_manual(values = method_colors, guide = "none") +
    scale_fill_manual(values = method_colors, guide = "none") +
    scale_shape_manual(
      values = c("p = 0.001" = 21, "p = 0.01" = 22, "p = 0.05" = 23),
      guide = "none"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = title, x = "False Positive Rate", y = "True Positive Rate") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
    )
}

scatter_panels <- mapply(
  make_scatter_panel,
  lvl = pert_labels,
  title = c("Two Peptides", "Random Peptides", "50% Peptides"),
  SIMPLIFY = FALSE
)

# ── Legends ───────────────────────────────────────────────────────────────────

legend_dummy <- data.table(
  method = factor(method_levels, levels = method_levels),
  x = seq_along(method_levels),
  y = seq_along(method_levels)
)
scatter_legend_method <- get_legend(
  ggplot(legend_dummy, aes(x = x, y = y, colour = method)) +
    geom_point() +
    scale_colour_manual(values = method_colors, name = "Method") +
    guides(
      colour = guide_legend(
        nrow = 1,
        override.aes = list(
          shape = 22,
          size = 5,
          fill = unname(method_colors[method_levels]),
          colour = "black"
        )
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm")
    )
)

scatter_legend_shape <- get_legend(
  ggplot(scatter_highlights, aes(x = FPR, y = TPR, shape = thr_label)) +
    geom_point(size = 3.5, fill = "grey60", colour = "black", stroke = 0.6) +
    scale_shape_manual(
      values = c("p = 0.001" = 21, "p = 0.01" = 22, "p = 0.05" = 23),
      name = "Threshold"
    ) +
    guides(shape = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm")
    )
)

# ── Assemble 1×3 grid ─────────────────────────────────────────────────────────

scatter_grid <- plot_grid(
  plotlist = scatter_panels,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_size = 13
)

scatter_legends_row <- plot_grid(
  scatter_legend_method,
  scatter_legend_shape,
  nrow = 1,
  rel_widths = c(1.5, 0.8)
)

p_scatter_fig <- plot_grid(
  scatter_grid,
  scatter_legends_row,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

print(p_scatter_fig)

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave(
  "00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Scatter.pdf",
  p_scatter_fig,
  width = 12,
  height = 4
)
ggsave(
  "00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Scatter.png",
  p_scatter_fig,
  width = 12,
  height = 4,
  dpi = 300
)

cat(
  "Saved: 00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Scatter.pdf/.png\n"
)

# ── Extended figure: scatter + combined F1 bar ────────────────────────────────

f1_at_001 <- scatter_highlights[thr_label == "p = 0.01"]

p_f1_bar <- ggplot(f1_at_001, aes(x = perturbation, y = F1, fill = method)) +
  geom_col(
    position = position_dodge(width = 0.85),
    width = 0.75,
    alpha = 0.9,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.2f", F1), color = method),
    position = position_dodge(width = 0.85),
    vjust = -0.35,
    size = 3.0,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1.2),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_color_manual(values = method_colors) +
  labs(title = "F1 Score (p = 0.01)", y = "F1") +
  bar_theme

p_scatter_extended <- plot_grid(
  scatter_grid,
  p_f1_bar,
  scatter_legends_row,
  ncol = 1,
  labels = c("", "D", ""),
  label_size = 13,
  rel_heights = c(1, 0.55, 0.08)
)

print(p_scatter_extended)

ggsave(
  "00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Scatter_Extended.pdf",
  p_scatter_extended,
  width = 12,
  height = 6
)
ggsave(
  "00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Scatter_Extended.png",
  p_scatter_extended,
  width = 12,
  height = 6,
  dpi = 300
)

cat(
  "Saved: 00-ThesisFigures/Figures/FindDCF/FindDCF_Grouping_Scatter_Extended.pdf/.png\n"
)
