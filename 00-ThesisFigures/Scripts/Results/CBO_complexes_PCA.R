library(data.table)
library(ggplot2)
library(cowplot)
library(scico)
library(arrow)
library(pcaMethods)

setwd(here::here())

output_dir <- "00-ThesisFigures/Figures/CBO_PCA"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ── Complex definitions ───────────────────────────────────────────────────────

complexes <- list(
  "B-WICH" = list(
    title = "B-WICH Chromatin Remodeling Complex",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/B-Wich_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_B-Wich_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_B-Wich_result.feather"
  ),
  "CRD" = list(
    title = "CRD-mediated mRNA stability complex",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/CRD-mediated_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_CRD-mediated_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_CRD-mediated_result.feather"
  ),
  "eIF3" = list(
    title = "Eukaryotic translation initiation factor 3",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Eukaryotic-translation-initiation-factor-3-complex_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_Eukaryotic-translation-initiation-factor-3-complex_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Eukaryotic-translation-initiation-factor-3-complex_result.feather"
  )
)

# ── Shared theme ──────────────────────────────────────────────────────────────

pca_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    aspect.ratio = 1
  )

# ── Helper: run ppca and return scores + var_expl ────────────────────────────

run_pca <- function(mat) {
  for (col in names(mat)) {
    mat[[col]][is.nan(mat[[col]])] <- NA
  }
  pca <- pcaMethods::pca(
    mat,
    method = "ppca",
    nPcs = 2,
    scale = "none",
    center = TRUE
  )
  list(scores = as.data.frame(pca@scores), var_expl = pca@R2)
}

# Maps raw cluster IDs to CF0, CF1, ... (numeric sort when possible, else lex).
# Anything already marked "CF-1" (singletons/unassigned) is kept as-is.
# Uses base:: qualifiers throughout to guard against BiocGenerics masking.
relabel_to_cf <- function(ids) {
  keep <- ids != "CF-1"
  unique_ids <- base::unique(ids[keep])
  num_vals <- base::suppressWarnings(base::as.numeric(unique_ids))
  sorted_ids <- if (!base::anyNA(num_vals)) {
    unique_ids[base::order(num_vals)]
  } else {
    base::sort(unique_ids)
  }
  cf_labels <- base::paste0("CF", base::seq_along(sorted_ids) - 1L)
  idx <- base::match(ids, sorted_ids)
  result <- base::ifelse(keep, cf_labels[idx], "CF-1")
  result
}

# scico requires n >= 2; this wrapper handles n == 1.
safe_scico <- function(n, palette = "batlow") {
  scico(max(n, 2L), palette = palette)[seq_len(n)]
}

prep_intensities <- function(data, intensity_cols) {
  long <- data.table::melt(
    data,
    id.vars = "Peptide",
    measure.vars = intensity_cols,
    variable.name = "Sample",
    value.name = "Intensity"
  )
  setDT(long)
  long[, Condition := sub("^B\\d+_", "", Sample)]
  long_avg <- long[,
    .(Intensity = mean(Intensity, na.rm = TRUE)),
    by = .(Peptide, Condition)
  ]
  wide <- dcast(long_avg, Peptide ~ Condition, value.var = "Intensity")
  mat <- as.data.frame(wide[, -"Peptide"])
  rownames(mat) <- wide[["Peptide"]]
  mat
}

# ── Build all plots ───────────────────────────────────────────────────────────

origin_plots <- list() # protein-of-origin  (supplementary)
dcf_plots <- list() # DCF assignment      (main row 1)
pf_plots <- list() # ProteoForge         (main row 2)

for (cname in names(complexes)) {
  paths <- complexes[[cname]]

  # ── Load raw data ──
  orig <- fread(paths$data)
  intensity_cols <- grep("^B\\d+_D\\d+$", names(orig), value = TRUE)

  # ── Protein-of-origin PCA ─────────────────────────────────────────────────
  pep_data <- orig[, c("Peptide", intensity_cols), with = FALSE]
  pep_data <- pep_data[rowSums(!is.na(pep_data[, ..intensity_cols])) > 0]

  mat_orig <- prep_intensities(pep_data, intensity_cols)
  res_orig <- run_pca(mat_orig)

  sc_orig <- res_orig$scores
  sc_orig$Acc <- rownames(sc_orig)
  sc_orig$Acc <- sub("^[^_]+_([^_]+)_.*$", "\\1", sc_orig$Acc)

  acc_levels <- sort(unique(sc_orig$Acc))
  pal_acc <- scico(length(acc_levels), palette = "batlow")
  acc_colors <- setNames(pal_acc, acc_levels)

  origin_plots[[cname]] <- ggplot(sc_orig, aes(PC1, PC2, color = Acc)) +
    geom_point(size = 2, alpha = 0.9) +
    scale_color_manual(values = acc_colors) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * res_orig$var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * res_orig$var_expl[2]),
      color = "Protein",
      title = paths$title
    ) +
    coord_fixed() +
    pca_theme

  # ── DCF PCA ───────────────────────────────────────────────────────────────
  if (file.exists(paths$dcf)) {
    dcf_res <- fread(paths$dcf)

    merged_dcf <- merge(
      dcf_res[, .(Peptide, dCF)],
      orig[, c("Peptide", intensity_cols), with = FALSE],
      by = "Peptide"
    )

    long_dcf <- data.table::melt(
      merged_dcf,
      id.vars = c("Peptide", "dCF"),
      measure.vars = intensity_cols,
      variable.name = "Sample",
      value.name = "Intensity"
    )
    setDT(long_dcf)
    long_dcf[, Condition := sub("^B\\d+_", "", Sample)]
    long_avg_dcf <- long_dcf[,
      .(Intensity = mean(Intensity, na.rm = TRUE)),
      by = .(Peptide, dCF, Condition)
    ]
    wide_dcf <- dcast(
      long_avg_dcf,
      Peptide + dCF ~ Condition,
      value.var = "Intensity"
    )

    cond_cols <- setdiff(names(wide_dcf), c("Peptide", "dCF"))
    mat_dcf <- as.data.frame(wide_dcf[, cond_cols, with = FALSE])
    rownames(mat_dcf) <- paste0(wide_dcf$Peptide, "_", wide_dcf$dCF)

    res_dcf <- run_pca(mat_dcf)
    sc_dcf <- res_dcf$scores
    sc_dcf$ID <- relabel_to_cf(sub(".*_", "", rownames(sc_dcf)))
    # DCF should not produce NA IDs; guard anyway
    sc_dcf <- sc_dcf[!is.na(sc_dcf$ID), ]

    id_levels_dcf <- sort(unique(sc_dcf$ID))
    pal_dcf <- safe_scico(length(id_levels_dcf))
    dcf_colors <- setNames(pal_dcf, id_levels_dcf)

    dcf_plots[[cname]] <- ggplot(sc_dcf, aes(PC1, PC2, color = ID)) +
      geom_point(size = 2, alpha = 0.9) +
      scale_color_manual(values = dcf_colors) +
      labs(
        x = sprintf("PC1 (%.1f%%)", 100 * res_dcf$var_expl[1]),
        y = sprintf("PC2 (%.1f%%)", 100 * res_dcf$var_expl[2]),
        color = "Assignment",
        title = paths$title
      ) +
      coord_fixed() +
      pca_theme
  }

  # ── ProteoForge PCA ───────────────────────────────────────────────────────
  if (file.exists(paths$pf)) {
    pf_res <- as.data.table(read_feather(paths$pf))

    pf_wide <- dcast(
      pf_res,
      peptide_id + ClusterID ~ Sample,
      value.var = "ms1adj",
      fun.aggregate = mean
    )
    pf_wide[, peptide_id := paste0(peptide_id, "_", ClusterID)]

    int_cols_pf <- grep("^B\\d+_D\\d+$", names(pf_wide), value = TRUE)
    long_pf <- data.table::melt(
      pf_wide,
      id.vars = "peptide_id",
      measure.vars = int_cols_pf,
      variable.name = "Sample",
      value.name = "Intensity"
    )
    setDT(long_pf)
    long_pf[, Condition := sub("^B\\d+_", "", Sample)]
    long_avg_pf <- long_pf[,
      .(Intensity = mean(Intensity, na.rm = TRUE)),
      by = .(peptide_id, Condition)
    ]
    data_pf <- dcast(
      long_avg_pf,
      peptide_id ~ Condition,
      value.var = "Intensity"
    )
    data_pf <- data_pf[!duplicated(peptide_id)]

    mat_pf <- as.data.frame(data_pf[, -"peptide_id"])
    rownames(mat_pf) <- data_pf[["peptide_id"]]

    res_pf <- run_pca(mat_pf)
    sc_pf <- res_pf$scores
    raw_ids <- sub(".*_", "", base::rownames(sc_pf))

    # NA ClusterID = unassigned peptide → CF-1
    raw_ids[is.na(raw_ids) | raw_ids == "NA"] <- "CF-1"

    cluster_counts <- table(raw_ids[raw_ids != "CF-1"])
    singletons <- names(cluster_counts[cluster_counts <= 1])
    raw_ids[raw_ids %in% singletons] <- "CF-1"

    sc_pf$ID <- relabel_to_cf(raw_ids)
    # CF-1 (singletons/unassigned) gets grey; CF0+ get the batlow palette
    cf_ids <- base::sort(base::unique(sc_pf$ID[sc_pf$ID != "CF-1"]))
    pal_pf <- safe_scico(length(cf_ids))
    pf_colors <- c(setNames(pal_pf, cf_ids), "CF-1" = "#AAAAAA")
    id_levels_pf <- c(cf_ids, "CF-1")

    sc_pf$ID <- factor(sc_pf$ID, levels = id_levels_pf)
    pf_plots[[cname]] <- ggplot(sc_pf, aes(PC1, PC2, color = ID)) +
      geom_point(size = 2, alpha = 0.9) +
      scale_color_manual(values = pf_colors) +
      labs(
        x = sprintf("PC1 (%.1f%%)", 100 * res_pf$var_expl[1]),
        y = sprintf("PC2 (%.1f%%)", 100 * res_pf$var_expl[2]),
        color = "Assignment",
        title = paths$title
      ) +
      coord_fixed() +
      pca_theme
  }
}

# ── Main figure: 2×3 grid (row = method, col = complex) ──────────────────────
#
# Row 1: ComplexoFinder (DCF)
# Row 2: ProteoForge
# Each cell has its own per-plot legend (legends vary by complex).

strip_legend <- function(p) p + theme(legend.position = "none")

row_label <- function(label) {
  ggdraw() +
    draw_label(label, fontface = "bold", size = 10, angle = 90, vjust = 0.5)
}

dcf_row <- plot_grid(
  plotlist = lapply(dcf_plots, function(p) {
    p + theme(legend.position = "right")
  }),
  nrow = 1,
  align = "hv",
  axis = "tblr"
)

pf_row <- plot_grid(
  plotlist = lapply(pf_plots, function(p) p + theme(legend.position = "right")),
  nrow = 1,
  align = "hv",
  axis = "tblr"
)

row_labels <- plot_grid(
  row_label("ComplexoFinder"),
  row_label("ProteoForge"),
  ncol = 1,
  rel_heights = c(1, 1)
)

grid_body <- plot_grid(
  dcf_row,
  pf_row,
  ncol = 1,
  rel_heights = c(1, 1),
  labels = c("A", "B"),
  label_size = 14,
  label_fontface = "bold"
)

main_figure <- plot_grid(
  row_labels,
  grid_body,
  nrow = 1,
  rel_widths = c(0.04, 1)
)

ggsave(
  filename = file.path(output_dir, "CBO_complexes_PCA_main.pdf"),
  plot = main_figure,
  width = 12,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_dir, "CBO_complexes_PCA_main.png"),
  plot = main_figure,
  width = 12,
  height = 6,
  dpi = 300
)

# ── Supplementary figure: protein-of-origin coloring ─────────────────────────

supp_row <- plot_grid(
  plotlist = lapply(origin_plots, function(p) {
    p + theme(legend.position = "right")
  }),
  nrow = 1,
  align = "hv",
  axis = "tblr"
)

ggsave(
  filename = file.path(output_dir, "CBO_complexes_PCA_supp.pdf"),
  plot = supp_row,
  width = 12,
  height = 3,
  dpi = 300
)
ggsave(
  filename = file.path(output_dir, "CBO_complexes_PCA_supp.png"),
  plot = supp_row,
  width = 12,
  height = 3,
  dpi = 300
)
