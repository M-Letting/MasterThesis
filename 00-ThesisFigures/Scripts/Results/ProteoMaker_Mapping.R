library(data.table)
library(ggplot2)
library(cowplot)
library(scico)
library(arrow)

setwd(here::here())

# ── Paths ─────────────────────────────────────────────────────────────────────

base_dir <- "Q-Benchmark/02-DCF"
output_dir <- "00-ThesisFigures/Figures/PM_mapping"

ProteoForge_NA_results <- list(
  LowNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNA_result.feather"
  ),
  MedNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNA_result.feather"
  ),
  HighNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNA_result.feather"
  )
)
DCF_NA_results <- list(
  LowNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerLowNA_results.csv"
  ),
  MedNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerMedNA_results.csv"
  ),
  HighNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerHighNA_results.csv"
  )
)
DCFBaseline_NA_results <- list(
  LowNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNA_results.csv"
  ),
  MedNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNA_results.csv"
  ),
  HighNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNA_results.csv"
  )
)
ProteoForge_NoNA_results <- list(
  LowNoNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNoNA_result.feather"
  ),
  MedNoNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNoNA_result.feather"
  ),
  HighNoNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNoNA_result.feather"
  )
)
DCF_NoNA_results <- list(
  LowNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerLowNoNA_results.csv"
  ),
  MedNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerMedNoNA_results.csv"
  ),
  HighNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerHighNoNA_results.csv"
  )
)
DCFBaseline_NoNA_results <- list(
  LowNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNoNA_results.csv"
  ),
  MedNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNoNA_results.csv"
  ),
  HighNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNoNA_results.csv"
  )
)
ProteoMaker_ref <- list(
  LowNoNA = file.path(
    base_dir,
    "00-Data/02-Prepared/ProteoMakerLowNoNA_wide_processed.csv"
  ),
  MedNoNA = file.path(
    base_dir,
    "00-Data/02-Prepared/ProteoMakerMedNoNA_wide_processed.csv"
  ),
  HighNoNA = file.path(
    base_dir,
    "00-Data/02-Prepared/ProteoMakerHighNoNA_wide_processed.csv"
  )
)

# ── Colours (matching FindDCF_Performance.R) ──────────────────────────────────

all_method_levels <- c(
  "ComplexoFinder",
  "ComplexoFinderBaseline",
  "ProteoForge",
  "PeCorA",
  "COPF"
)
pal <- scico(5, palette = "batlow", direction = 1)
all_method_colors <- setNames(pal, all_method_levels)

plot_methods <- c(
  "ComplexoFinder",
  "ComplexoFinderBaseline",
  "ProteoForge"
)
method_colors <- all_method_colors[plot_methods]

# ── Helper ────────────────────────────────────────────────────────────────────

count_non_singleton_proteoforms <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & ids != ""]
  if (!length(ids)) {
    return(0L)
  }
  split_ids <- unlist(strsplit(ids, "\\|"), use.names = FALSE)
  split_ids <- split_ids[split_ids != ""]
  if (!length(split_ids)) {
    return(0L)
  }
  as.integer(length(unique(split_ids)))
}

# ── Data loading ──────────────────────────────────────────────────────────────
# ref_name_fn: function(dataset_name) -> key into ProteoMaker_ref
# For NA datasets: ref is the matching NoNA file (sub "NA" -> "NoNA")
# For NoNA datasets: ref is the file itself (identity)

load_counts <- function(
  pf_results,
  dcf_results,
  dcfbl_results,
  ref_name_fn
) {
  dt <- data.table(
    Dataset = character(),
    Method = character(),
    Protein = character(),
    FoundProteoformCount = integer(),
    ProteoMakerProteoformCount = integer()
  )

  # ProteoForge
  for (name in names(pf_results)) {
    path <- pf_results[[name]]
    if (!file.exists(path)) {
      warning("File not found: ", path)
      next
    }

    pf_res <- as.data.table(read_feather(path))[Sample == "C1_R1"]
    clust_sizes <- pf_res[, .(n = .N), by = .(protein_id, ClusterID)]
    pf_counts <- clust_sizes[
      n > 1,
      .(FoundProteoformCount = .N),
      by = protein_id
    ]
    pf_counts <- pf_counts[!grepl("\\|", protein_id)]

    ref_data <- fread(ProteoMaker_ref[[ref_name_fn(name)]])
    for (prot in pf_counts$protein_id) {
      if (!(prot %in% ref_data$Accession)) {
        next
      }
      pf_counts[
        protein_id == prot,
        ProteoMakerProteoformCount := count_non_singleton_proteoforms(
          ref_data[Accession == prot, Proteoform_ID]
        )
      ]
    }
    pf_counts[, `:=`(Dataset = name, Method = "ProteoForge")]
    setnames(pf_counts, "protein_id", "Protein")
    dt <- rbind(
      dt,
      pf_counts[, .(
        Dataset,
        Method,
        Protein,
        FoundProteoformCount,
        ProteoMakerProteoformCount
      )],
      fill = TRUE
    )
  }

  # ComplexoFinder (DCF)
  for (name in names(dcf_results)) {
    path <- dcf_results[[name]]
    if (!file.exists(path)) {
      warning("File not found: ", path)
      next
    }

    dcf_res <- fread(path)
    dcf_sizes <- dcf_res[
      !is.na(dCF) & dCF != "dCF-1",
      .N,
      by = .(Accession, dCF)
    ]
    dcf_counts <- dcf_sizes[N > 1, .(FoundProteoformCount = .N), by = Accession]
    dcf_counts <- dcf_counts[!grepl("\\|", Accession)]

    ref_data <- fread(ProteoMaker_ref[[ref_name_fn(name)]])
    for (prot in dcf_counts$Accession) {
      if (!(prot %in% ref_data$Accession)) {
        next
      }
      dcf_counts[
        Accession == prot,
        ProteoMakerProteoformCount := count_non_singleton_proteoforms(
          ref_data[Accession == prot, Proteoform_ID]
        )
      ]
    }
    dcf_counts[, `:=`(Dataset = name, Method = "ComplexoFinder")]
    setnames(dcf_counts, "Accession", "Protein")
    dt <- rbind(
      dt,
      dcf_counts[, .(
        Dataset,
        Method,
        Protein,
        FoundProteoformCount,
        ProteoMakerProteoformCount
      )],
      fill = TRUE
    )
  }

  # ComplexoFinderBaseline (DCFBaseline)
  for (name in names(dcfbl_results)) {
    path <- dcfbl_results[[name]]
    if (!file.exists(path)) {
      warning("File not found: ", path)
      next
    }

    dcf_res <- fread(path)
    dcf_sizes <- dcf_res[
      !is.na(dCF) & dCF != "dCF-1",
      .N,
      by = .(Accession, dCF)
    ]
    dcf_counts <- dcf_sizes[N > 1, .(FoundProteoformCount = .N), by = Accession]
    dcf_counts <- dcf_counts[!grepl("\\|", Accession)]

    ref_data <- fread(ProteoMaker_ref[[ref_name_fn(name)]])
    for (prot in dcf_counts$Accession) {
      if (!(prot %in% ref_data$Accession)) {
        next
      }
      dcf_counts[
        Accession == prot,
        ProteoMakerProteoformCount := count_non_singleton_proteoforms(
          ref_data[Accession == prot, Proteoform_ID]
        )
      ]
    }
    dcf_counts[, `:=`(Dataset = name, Method = "ComplexoFinderBaseline")]
    setnames(dcf_counts, "Accession", "Protein")
    dt <- rbind(
      dt,
      dcf_counts[, .(
        Dataset,
        Method,
        Protein,
        FoundProteoformCount,
        ProteoMakerProteoformCount
      )],
      fill = TRUE
    )
  }

  dt
}

# ── Shared theme ──────────────────────────────────────────────────────────────

base_theme <- theme_classic(base_size = 10) +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
  )

prop_colors <- c(Under = "#d73027", Exact = "#4dac26", Over = "#4575b4")

x_label_short <- c(
  ComplexoFinder = "ComplexoFinder",
  ComplexoFinderBaseline = "Baseline",
  ProteoForge = "ProteoForge"
)

# ── Plotting functions ────────────────────────────────────────────────────────

make_shared_method_legend <- function(counts_dt) {
  get_legend(
    ggplot(counts_dt, aes(x = diff, fill = Method)) +
      geom_bar(alpha = 0.9, color = "black", linewidth = 0.25) +
      scale_fill_manual(values = method_colors, name = "Method") +
      theme(
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm")
      ) +
      guides(fill = guide_legend(nrow = 1))
  )
}

make_prop_legend <- function(prop_dt) {
  get_legend(
    ggplot(prop_dt, aes(x = Dataset, fill = Category)) +
      geom_bar(position = "fill") +
      scale_fill_manual(values = prop_colors, name = "Recovery") +
      theme(
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm")
      ) +
      guides(fill = guide_legend(nrow = 1))
  )
}

make_bar_panel <- function(counts_dt, method, dataset, bar_y_max) {
  sub <- counts_dt[Method == method & Dataset == dataset]
  col <- method_colors[[as.character(method)]]
  mad_diff <- mean(abs(sub$diff - mean(sub$diff, na.rm = TRUE)), na.rm = TRUE)
  mean_label <- sprintf("MAD = %s", formatC(round(mad_diff, 2), format = "f", digits = 2))

  ggplot(sub, aes(x = diff)) +
    geom_bar(
      fill = col,
      alpha = 0.9,
      color = "black",
      linewidth = 0.2,
      width = 0.8
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = mean_label,
      hjust = -0.05, vjust = 1.4,
      size = 3,
      color = "grey20"
    ) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
    scale_y_continuous(
      limits = c(0, bar_y_max * 1.08),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(y = "Number of proteins") +
    base_theme +
    theme(axis.title.x = element_blank())
}

make_violin_panel <- function(
  counts_dt,
  dataset,
  violin_y_range,
  show_y_title
) {
  sub <- counts_dt[Dataset == dataset]

  ggplot(sub, aes(x = Method, y = diff, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_violin(alpha = 0.75, color = NA, trim = TRUE) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = "white",
      color = "grey30",
      linewidth = 0.5
    ) +
    scale_fill_manual(values = method_colors) +
    scale_x_discrete(labels = x_label_short) +
    scale_y_continuous(limits = violin_y_range) +
    labs(
      title = dataset,
      x = NULL,
      y = if (show_y_title) "Found − ProteoMaker count" else NULL
    ) +
    base_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
      axis.title.y = if (show_y_title) {
        element_text(size = 9)
      } else {
        element_blank()
      }
    )
}

make_prop_panel <- function(prop_dt, dataset, show_y_title) {
  sub <- prop_dt[Dataset == dataset]

  ggplot(sub, aes(x = Method, fill = Category)) +
    geom_bar(position = "fill", alpha = 0.9, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = prop_colors) +
    scale_x_discrete(labels = x_label_short) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = dataset,
      x = NULL,
      y = if (show_y_title) "Proportion of proteins" else NULL
    ) +
    base_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
      axis.title.y = if (show_y_title) {
        element_text(size = 9)
      } else {
        element_blank()
      }
    )
}

# ── Main pipeline ─────────────────────────────────────────────────────────────

run_pipeline <- function(
  counts_raw,
  dataset_levels,
  dataset_labels,
  file_suffix
) {
  counts_dt <- copy(counts_raw)
  counts_dt[, Method := factor(Method, levels = plot_methods)]
  counts_dt[,
    Dataset := factor(Dataset, levels = dataset_levels, labels = dataset_labels)
  ]
  counts_dt[, diff := FoundProteoformCount - ProteoMakerProteoformCount]
  counts_dt <- counts_dt[is.finite(diff)]

  shared_legend <- make_shared_method_legend(counts_dt)

  # ── Plot A: 3×3 bar grid (rows = methods, cols = datasets) ──────────────────

  bar_y_max <- counts_dt[, max(table(diff)), by = .(Method, Dataset)][, max(V1)]
  row_lbl_w <- 0.08

  bar_panels <- list()
  for (i in seq_along(plot_methods)) {
    for (j in seq_along(dataset_labels)) {
      bar_panels[[paste(i, j)]] <- make_bar_panel(
        counts_dt,
        plot_methods[[i]],
        dataset_labels[[j]],
        bar_y_max
      )
    }
  }

  make_row_label <- function(method) {
    ggdraw() +
      draw_label(method, fontface = "bold", size = 9, angle = 90, hjust = 0.5)
  }

  make_bar_row <- function(idx) {
    row_plot <- plot_grid(
      plotlist = bar_panels[grep(paste0("^", idx, " "), names(bar_panels))],
      nrow = 1
    )
    plot_grid(
      make_row_label(plot_methods[[idx]]),
      row_plot,
      nrow = 1,
      rel_widths = c(row_lbl_w, 1)
    )
  }

  col_header_spacer <- plot_grid(
    ggdraw(),
    plot_grid(
      plotlist = lapply(dataset_labels, function(d) {
        ggdraw() + draw_label(d, fontface = "bold", size = 11, hjust = 0.5)
      }),
      nrow = 1
    ),
    nrow = 1,
    rel_widths = c(row_lbl_w, 1)
  )
  x_title_row <- plot_grid(
    ggdraw(),
    ggdraw() +
      draw_label(
        "Difference between number of identified proteoforms and ProteoMaker proteoforms",
        size = 9,
        hjust = 0.5
      ),
    nrow = 1,
    rel_widths = c(row_lbl_w, 1)
  )

  p_bars_final <- plot_grid(
    plot_grid(
      plotlist = c(
        list(col_header_spacer),
        lapply(seq_along(plot_methods), make_bar_row),
        list(x_title_row)
      ),
      ncol = 1,
      rel_heights = c(0.15, rep(1, length(plot_methods)), 0.10)
    ),
    shared_legend,
    ncol = 1,
    rel_heights = c(1, 0.07)
  )

  # ── Plot B: Violin + boxplot ─────────────────────────────────────────────────

  violin_y_range <- counts_dt[, range(diff, na.rm = TRUE)]

  violin_panels <- mapply(
    make_violin_panel,
    dataset = dataset_labels,
    show_y_title = c(TRUE, FALSE, FALSE),
    MoreArgs = list(counts_dt = counts_dt, violin_y_range = violin_y_range),
    SIMPLIFY = FALSE
  )

  p_violin_final <- plot_grid(
    plot_grid(plotlist = violin_panels, nrow = 1),
    shared_legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )

  # ── Plot C: Stacked proportional bar ─────────────────────────────────────────

  prop_dt <- counts_dt[,
    .(
      Category = fcase(
        diff < 0 , "Under" , diff == 0 , "Exact" , diff > 0 , "Over"
      )
    ),
    by = .(Dataset, Method, Protein)
  ]
  prop_dt[, Category := factor(Category, levels = c("Over", "Exact", "Under"))]

  prop_panels <- mapply(
    make_prop_panel,
    dataset = dataset_labels,
    show_y_title = c(TRUE, FALSE, FALSE),
    MoreArgs = list(prop_dt = prop_dt),
    SIMPLIFY = FALSE
  )

  p_prop_final <- plot_grid(
    plot_grid(plotlist = prop_panels, nrow = 1),
    make_prop_legend(prop_dt),
    ncol = 1,
    rel_heights = c(1, 0.12)
  )

  # ── Save ──────────────────────────────────────────────────────────────────────

  ggsave(
    file.path(
      output_dir,
      paste0("ProteoMaker_Mapping_Bars_", file_suffix, ".pdf")
    ),
    p_bars_final,
    width = 10,
    height = 5
  )
  ggsave(
    file.path(
      output_dir,
      paste0("ProteoMaker_Mapping_Bars_", file_suffix, ".png")
    ),
    p_bars_final,
    width = 10,
    height = 5,
    dpi = 300
  )
  ggsave(
    file.path(
      output_dir,
      paste0("ProteoMaker_Mapping_Violin_", file_suffix, ".pdf")
    ),
    p_violin_final,
    width = 10,
    height = 5
  )
  ggsave(
    file.path(
      output_dir,
      paste0("ProteoMaker_Mapping_Violin_", file_suffix, ".png")
    ),
    p_violin_final,
    width = 10,
    height = 5,
    dpi = 300
  )
  ggsave(
    file.path(
      output_dir,
      paste0("ProteoMaker_Mapping_Proportional_", file_suffix, ".pdf")
    ),
    p_prop_final,
    width = 10,
    height = 4.5
  )
  ggsave(
    file.path(
      output_dir,
      paste0("ProteoMaker_Mapping_Proportional_", file_suffix, ".png")
    ),
    p_prop_final,
    width = 10,
    height = 4.5,
    dpi = 300
  )

  cat(
    "Saved ProteoMaker_Mapping_*_",
    file_suffix,
    ".pdf/.png to ",
    output_dir,
    "\n",
    sep = ""
  )
}

# ── Run for NA datasets ───────────────────────────────────────────────────────

counts_NA <- load_counts(
  pf_results = ProteoForge_NA_results,
  dcf_results = DCF_NA_results,
  dcfbl_results = DCFBaseline_NA_results,
  ref_name_fn = function(name) sub("NA$", "NoNA", name)
)

run_pipeline(
  counts_raw = counts_NA,
  dataset_levels = c("LowNA", "MedNA", "HighNA"),
  dataset_labels = c("Low NA", "Medium NA", "High NA"),
  file_suffix = "NA"
)

# ── Run for NoNA datasets ─────────────────────────────────────────────────────

counts_NoNA <- load_counts(
  pf_results = ProteoForge_NoNA_results,
  dcf_results = DCF_NoNA_results,
  dcfbl_results = DCFBaseline_NoNA_results,
  ref_name_fn = function(name) name
)

run_pipeline(
  counts_raw = counts_NoNA,
  dataset_levels = c("LowNoNA", "MedNoNA", "HighNoNA"),
  dataset_labels = c("Low NoNA", "Medium NoNA", "High NoNA"),
  file_suffix = "NoNA"
)
