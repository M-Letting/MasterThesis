library(data.table)
library(ggplot2)
library(cowplot)
library(scico)

setwd(here::here())

# ── Data ──────────────────────────────────────────────────────────────────────

data <- fread(
  "Q-Benchmark/02-DCF/00-Data/04-PerformanceMetrics/Combined_NA_Metrics.csv"
)

# Remove rows with "NoNA" in Dataset
data <- data[!grepl("NoNA", Dataset)]

# For each method × NA level find the row with the highest F1 score
optimal <- data[
  !is.nan(F1) & !is.na(F1) & is.finite(F1),
  .SD[which.max(F1)],
  by = .(method, NALevel)
]

# ── Factor ordering ───────────────────────────────────────────────────────────

na_levels <- c("Low", "Med", "High")
na_labels <- c("Low", "Medium", "High")

method_levels <- c(
  "ComplexoFinder",
  "ComplexoFinderBaseline",
  "ProteoForge",
  "PeCorA",
  "COPF"
)

optimal[, NALevel := factor(NALevel, levels = na_levels, labels = na_labels)]
optimal[, method := factor(method, levels = method_levels)]

# ── Colours ───────────────────────────────────────────────────────────────────

pal <- scico(5, palette = "batlow", direction = 1)
method_colors <- setNames(pal, method_levels)
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

p_tpr <- ggplot(optimal, aes(x = NALevel, y = TPR, fill = method)) +
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

# ── FDR bar chart ─────────────────────────────────────────────────────────────

p_fdr <- ggplot(optimal, aes(x = NALevel, y = FDR, fill = method)) +
  geom_col(
    position = position_dodge(width = 0.85),
    width = 0.75,
    alpha = 0.9,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.2f", FDR), color = method),
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
  labs(title = "False Discovery Rate", y = "FDR") +
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
    guides(fill = guide_legend(nrow = 1))
)

# ── Threshold annotation table ────────────────────────────────────────────────

thresh_annot <- {
  thresh <- optimal[, .(method, NALevel, val = Threshold)]
  thresh[, method := factor(method, levels = rev(method_levels))]

  ggplot(
    thresh,
    aes(x = NALevel, y = method, label = formatC(val, format = "e", digits = 1))
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
    "Discordant Peptide Identification — ProteoMaker Datasets",
    fontface = "bold",
    size = 14,
    hjust = 0.5
  )

combined <- plot_grid(
  fig_title,
  p_tpr,
  p_fdr,
  legend,
  thresh_annot,
  ncol = 1,
  rel_heights = c(0.07, 1, 1, 0.12, 0.45),
  labels = c("", "A", "B", "", "C"),
  label_size = 13
)

print(combined)

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave(
  "00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_Performance.pdf",
  combined,
  width = 10,
  height = 10
)
ggsave(
  "00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_Performance.png",
  combined,
  width = 10,
  height = 10,
  dpi = 300
)

cat(
  "Saved: 00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_Performance.pdf/.png\n"
)

# ── Scatter: 2×3 combined grid (NA row A + NoNA row B) ───────────────────────

all_scatter <- fread(
  "Q-Benchmark/02-DCF/00-Data/04-PerformanceMetrics/Combined_NA_Metrics.csv"
)

# NA scatter data
na_scatter <- all_scatter[!grepl("NoNA", Dataset)]
na_scatter[,
  method := fcase(
    method == "DCF"         , "ComplexoFinder"         ,
    method == "DCFBaseline" , "ComplexoFinderBaseline" ,
    default = method
  )
]
na_scatter[, NALevel := factor(NALevel, levels = na_levels, labels = na_labels)]
na_scatter[, method := factor(method, levels = method_levels)]
na_scatter <- na_scatter[
  !is.nan(TPR) &
    !is.na(TPR) &
    !is.nan(FDR) &
    !is.na(FDR) &
    is.finite(TPR) &
    is.finite(FDR) &
    !is.na(NALevel)
]

# NoNA scatter data
nona_scatter <- all_scatter[grepl("NoNA", Dataset)]
nona_scatter[,
  NALevel := fcase(
    NALevel == "LowNo"  , "Low"  ,
    NALevel == "MedNo"  , "Med"  ,
    NALevel == "HighNo" , "High"
  )
]
nona_scatter[,
  NALevel := factor(NALevel, levels = na_levels, labels = na_labels)
]
nona_scatter[, method := factor(method, levels = method_levels)]
nona_scatter <- nona_scatter[
  !is.nan(TPR) &
    !is.na(TPR) &
    !is.nan(FDR) &
    !is.na(FDR) &
    is.finite(TPR) &
    is.finite(FDR) &
    !is.na(NALevel)
]

# Highlight points at p = 0.001, 0.01, 0.05
pval_highlights <- c(0.001, 0.01, 0.05)

make_highlights <- function(sdata) {
  out <- rbindlist(lapply(pval_highlights, function(thr) {
    sdata[,
      .SD[which.min(abs(Threshold - thr))],
      by = .(method, NALevel)
    ][, thr_label := paste0("p = ", thr)]
  }))
  out[,
    thr_label := factor(thr_label, levels = paste0("p = ", pval_highlights))
  ]
  out
}

na_highlights <- make_highlights(na_scatter)
nona_highlights <- make_highlights(nona_scatter)

# Panel builder
make_panel <- function(sdata, highlights, lvl, title) {
  ggplot(sdata[NALevel == lvl], aes(x = TPR, y = 1 - FDR, colour = method)) +
    geom_point(size = 1.5, alpha = 0.5, shape = 16) +
    geom_point(
      data = highlights[NALevel == lvl],
      aes(x = TPR, y = 1 - FDR, shape = thr_label, fill = method),
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
    labs(title = title, x = "TPR", y = "1 - FDR") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
    )
}

lvl_ordered <- c("Low", "Medium", "High")
na_titles <- c("Low — 50% NA", "Medium — 50% NA", "High — 50% NA")
nona_titles <- c("Low — No NA", "Medium — No NA", "High — No NA")

na_panels <- mapply(
  make_panel,
  lvl = lvl_ordered,
  title = na_titles,
  MoreArgs = list(sdata = na_scatter, highlights = na_highlights),
  SIMPLIFY = FALSE
)
nona_panels <- mapply(
  make_panel,
  lvl = lvl_ordered,
  title = nona_titles,
  MoreArgs = list(sdata = nona_scatter, highlights = nona_highlights),
  SIMPLIFY = FALSE
)

# ── Legends ───────────────────────────────────────────────────────────────────

# Method legend — dummy dataset ensures COPF (no NA scatter points) appears
legend_dummy <- data.table(
  method = factor(method_levels, levels = method_levels),
  x = seq_along(method_levels),
  y = seq_along(method_levels)
)
legend_method <- get_legend(
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

legend_shape <- get_legend(
  ggplot(na_highlights, aes(x = TPR, y = 1 - FDR, shape = thr_label)) +
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

# ── Shared legends row ────────────────────────────────────────────────────────

legends_row <- plot_grid(
  legend_method,
  legend_shape,
  nrow = 1,
  rel_widths = c(1.5, 0.8)
)

# ── F1 bar builder ────────────────────────────────────────────────────────────

make_f1_bar <- function(highlights) {
  f1_dat <- highlights[thr_label == "p = 0.01"]
  ggplot(f1_dat, aes(x = NALevel, y = F1, fill = method)) +
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
}

# ── Figure 1: 50% NA ──────────────────────────────────────────────────────────

row_na <- plot_grid(
  plotlist = na_panels,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_size = 13
)

p_na_extended <- plot_grid(
  row_na,
  make_f1_bar(na_highlights),
  legends_row,
  ncol = 1,
  labels = c("", "D", ""),
  label_size = 13,
  rel_heights = c(0.8, 0.55, 0.08)
)

print(p_na_extended)

ggsave(
  "00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_NA_Extended.pdf",
  p_na_extended,
  width = 12,
  height = 7
)
ggsave(
  "00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_NA_Extended.png",
  p_na_extended,
  width = 12,
  height = 7,
  dpi = 300
)

cat(
  "Saved: 00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_NA_Extended.pdf/.png\n"
)

# ── Figure 2: No NA ───────────────────────────────────────────────────────────

row_nona <- plot_grid(
  plotlist = nona_panels,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_size = 13
)

p_nona_extended <- plot_grid(
  row_nona,
  make_f1_bar(nona_highlights),
  legends_row,
  ncol = 1,
  labels = c("", "D", ""),
  label_size = 13,
  rel_heights = c(0.8, 0.55, 0.08)
)

print(p_nona_extended)

ggsave(
  "00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_NoNA_Extended.pdf",
  p_nona_extended,
  width = 12,
  height = 7
)
ggsave(
  "00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_NoNA_Extended.png",
  p_nona_extended,
  width = 12,
  height = 7,
  dpi = 300
)

cat(
  "Saved: 00-ThesisFigures/Figures/DCF_PM/DCF_ProteoMaker_NoNA_Extended.pdf/.png\n"
)
