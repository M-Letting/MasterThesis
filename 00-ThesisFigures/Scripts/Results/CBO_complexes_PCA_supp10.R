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

# ── Complex definitions (remaining 10 not shown in main figure) ───────────────

complexes <- list(
  "Multiaminoacyl" = list(
    title = "Multiaminoacyl-tRNA Synthetase Complex",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Multiaminoacyl-tRNA_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_Multiaminoacyl-tRNA_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Multiaminoacyl-tRNA_result.feather"
  ),
  "26S-proteasome" = list(
    title = "26S Proteasome",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/26S-proteasome_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_26S-proteasome_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_26S-proteasome_result.feather"
  ),
  "60S-ribosome" = list(
    title = "60S Cytosolic Large Ribosomal Subunit",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/60S-cytosolic-large-ribosomal-subunit_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_60S-cytosolic-large-ribosomal-subunit_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_60S-cytosolic-large-ribosomal-subunit_result.feather"
  ),
  "BSSADCRC" = list(
    title = "BSSADCRC Complex",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/BSSADCRC_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_BSSADCRC_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_BSSADCRC_result.feather"
  ),
  "Dynactin" = list(
    title = "Dynactin Complex",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Dynactin-complex_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_Dynactin-complex_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Dynactin-complex_result.feather"
  ),
  "Dynein-1" = list(
    title = "Dynein-1 Complex (variant 1)",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Dynein-1-complex-variant-1_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_Dynein-1-complex-variant-1_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Dynein-1-complex-variant-1_result.feather"
  ),
  "IFT-B" = list(
    title = "Intraflagellar Transport Complex B",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Intraflagellar-transport-complex-B_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_Intraflagellar-transport-complex-B_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Intraflagellar-transport-complex-B_result.feather"
  ),
  "Laminin-213" = list(
    title = "Laminin-213",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Laminin-213_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_Laminin-213_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Laminin-213_result.feather"
  ),
  "Spliceosomal-B" = list(
    title = "Major Spliceosomal B Complex",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Major-Spliceosomal-B_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_Major-Spliceosomal-B_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Major-Spliceosomal-B_result.feather"
  ),
  "Nuclear-pore" = list(
    title = "Nuclear Pore Complex",
    data = "Q-Benchmark/02-DCF/00-Data/02-Prepared/Nuclear-pore-complex_WideLogNorm.csv",
    dcf = "Q-Benchmark/02-DCF/00-Data/03-Results/DCF/DCF_Complex_Nuclear-pore-complex_results.csv",
    pf = "Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge/ProteoForge_Nuclear-pore-complex_result.feather"
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

# ── Helpers ───────────────────────────────────────────────────────────────────

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
  base::ifelse(keep, cf_labels[idx], "CF-1")
}

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

dcf_plots <- list()
pf_plots <- list()

for (cname in names(complexes)) {
  paths <- complexes[[cname]]

  orig <- fread(paths$data)
  intensity_cols <- grep("^B\\d+_D\\d+$", names(orig), value = TRUE)

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
    sc_dcf <- sc_dcf[!is.na(sc_dcf$ID), ]

    id_levels_dcf <- sort(unique(sc_dcf$ID))
    dcf_colors <- setNames(safe_scico(length(id_levels_dcf)), id_levels_dcf)

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

    raw_ids[is.na(raw_ids) | raw_ids == "NA"] <- "CF-1"

    cluster_counts <- table(raw_ids[raw_ids != "CF-1"])
    singletons <- names(cluster_counts[cluster_counts <= 1])
    raw_ids[raw_ids %in% singletons] <- "CF-1"

    sc_pf$ID <- relabel_to_cf(raw_ids)
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

# ── Assemble supplementary figure: 2 panes (A above B), each 2×5 grid ────────
#
# Each pane: row 1 = ComplexoFinder (DCF), row 2 = ProteoForge
# Pane A: first 5 complexes | Pane B: last 5 complexes

cnames <- names(complexes)
group_a <- cnames[1:5]
group_b <- cnames[6:10]

row_label <- function(label) {
  ggdraw() +
    draw_label(label, fontface = "bold", size = 10, angle = 90, vjust = 0.5)
}

make_pane <- function(group) {
  dcf_row <- plot_grid(
    plotlist = lapply(dcf_plots[group], function(p) {
      p + theme(legend.position = "right")
    }),
    nrow = 1,
    align = "hv",
    axis = "tblr"
  )
  pf_row <- plot_grid(
    plotlist = lapply(pf_plots[group], function(p) {
      p + theme(legend.position = "right")
    }),
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
    rel_heights = c(1, 1)
  )
  plot_grid(
    row_labels,
    grid_body,
    nrow = 1,
    rel_widths = c(0.04, 1)
  )
}

pane_a <- make_pane(group_a)
pane_b <- make_pane(group_b)

supp_figure <- plot_grid(
  pane_a,
  pane_b,
  ncol = 1,
  labels = c("", ""),
  label_size = 14,
  label_fontface = "bold"
)

ggsave(
  filename = file.path(output_dir, "CBO_complexes_PCA_supp10.pdf"),
  plot = supp_figure,
  width = 20,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output_dir, "CBO_complexes_PCA_supp10.png"),
  plot = supp_figure,
  width = 20,
  height = 12,
  dpi = 300
)
