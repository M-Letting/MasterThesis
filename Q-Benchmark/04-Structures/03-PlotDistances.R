# Libraries
library(data.table)
library(ggplot2)

# Set working directory to project root
setwd(here::here())

# Set output directory for figures
output_dir <- "Q-Benchmark/04-Structures/00-Plots"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

################################################################################
## Function to plot distance ###################################################
################################################################################

plot_distance <- function(file_path, title) {
  data <- fread(file_path)
  data <- data[!is.na(Distance_Angstrom)]
  if (!"SameCoordFrame" %in% names(data)) {
    data[, SameCoordFrame := TRUE]
  }
  data <- data[SameCoordFrame == TRUE]
  data[, SameCoordFrame := NULL]

  # Build unique peptide labels: Gene_Residue_Datatype
  make_label <- function(gene, residue, datatype) {
    paste(gene, residue, datatype, sep = "_")
  }

  data[, label1 := make_label(Gene1, Residue1, Datatype1)]
  data[, label2 := make_label(Gene2, Residue2, Datatype2)]

  # Collect all unique peptides with their gene membership
  peps1 <- unique(data[, .(label = label1, gene = Gene1, residue = Residue1)])
  peps2 <- unique(data[, .(label = label2, gene = Gene2, residue = Residue2)])
  all_peps <- unique(rbind(peps1, peps2))

  # Order: group by gene (gene order determined by first appearance in data),
  # then by residue number within each gene
  gene_order <- unique(c(data$Gene1, data$Gene2))
  all_peps[, gene := factor(gene, levels = gene_order)]
  setorder(all_peps, gene, residue)
  pep_levels <- all_peps$label

  # Symmetrize: mirror every row so distance(A,B) == distance(B,A)
  mirrored <- data[, .(label1 = label2, label2 = label1, Distance_Angstrom)]
  # Diagonal: every peptide has distance 0 to itself
  self_pairs <- data.table(
    label1 = pep_levels,
    label2 = pep_levels,
    Distance_Angstrom = 0
  )
  data <- rbind(
    data[, .(label1, label2, Distance_Angstrom)],
    mirrored,
    self_pairs
  )

  data[, label1 := factor(label1, levels = pep_levels)]
  data[, label2 := factor(label2, levels = pep_levels)]

  # Collapse duplicates by taking minimum distance (closest observed structure)
  data <- data[,
    .(Distance_Angstrom = min(Distance_Angstrom)),
    by = .(label1, label2)
  ]

  # Complete the full grid so unmeasured pairs render as NA (not invisible)
  full_grid <- CJ(label1 = pep_levels, label2 = pep_levels)
  full_grid[, label1 := factor(label1, levels = pep_levels)]
  full_grid[, label2 := factor(label2, levels = pep_levels)]
  data <- merge(full_grid, data, by = c("label1", "label2"), all.x = TRUE)

  # Gene boundary positions (for separator lines)
  n_peps <- length(pep_levels)
  gene_sizes <- all_peps[, .N, by = gene][order(match(gene, gene_order))]
  x_breaks <- head(cumsum(gene_sizes$N), -1) + 0.5
  # y-axis is reversed: coordinate of pep_levels[i] is n_peps - i + 1
  y_breaks <- n_peps - head(cumsum(gene_sizes$N), -1) + 0.5

  # Gene label: one tick per gene, placed at the midpoint level name
  mid_idx <- round(
    cumsum(c(0, head(gene_sizes$N, -1))) + gene_sizes$N / 2 + 0.5
  )
  mid_idx <- pmin(pmax(mid_idx, 1L), n_peps)
  gene_labels <- all_peps$label[mid_idx]

  p <- ggplot(data, aes(x = label1, y = label2, fill = Distance_Angstrom)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("#590007", "#932E05", "#C27142", "#DBAB8E", "#EBE5E0"),
      name = "Distance (Å)",
      na.value = "white"
    ) +
    # Gene boundary lines
    geom_vline(xintercept = x_breaks, colour = "black", linewidth = 0.6) +
    geom_hline(yintercept = y_breaks, colour = "black", linewidth = 0.6) +
    # Gene labels on top (x) and left (y)
    scale_x_discrete(
      breaks = gene_labels,
      labels = gene_order,
      position = "top"
    ) +
    scale_y_discrete(
      breaks = gene_labels,
      labels = gene_order,
      limits = rev(pep_levels)
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, size = 7, vjust = 0),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 15),
      legend.position = "bottom"
    )

  return(p)
}

################################################################################
## Define input ################################################################
################################################################################

CF_prefix <- "Q-Benchmark/04-Structures/00-Data/Results/ComplexoFinder"
PF_prefix <- "Q-Benchmark/04-Structures/00-Data/Results/ProteoForge"

CF_files <- list(
  "26S-proteasome" = list(
    C1 = file.path(
      CF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_cluster_1_ptm_distances.csv"
    ),
    C2 = file.path(
      CF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_cluster_2_ptm_distances.csv"
    ),
    C3 = file.path(
      CF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_cluster_3_ptm_distances.csv"
    ),
    C4 = file.path(
      CF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_cluster_4_ptm_distances.csv"
    ),
    C5 = file.path(
      CF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_cluster_5_ptm_distances.csv"
    )
  ),
  "Brain-SWI-SNF" = list(
    C1 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_1_ptm_distances.csv"
    ),
    C2 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_2_ptm_distances.csv"
    ),
    C3 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_3_ptm_distances.csv"
    ),
    C4 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_4_ptm_distances.csv"
    ),
    C5 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_5_ptm_distances.csv"
    ),
    C6 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_6_ptm_distances.csv"
    ),
    C7 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_7_ptm_distances.csv"
    )
  ),
  "Dynein-1-complex-variant-1" = list(
    C1 = file.path(
      CF_prefix,
      "Dynein_1_complex_variant_4/distanceFiles/Dynein_1_complex_variant_4_cluster_1_ptm_distances.csv"
    ),
    C2 = file.path(
      CF_prefix,
      "Dynein_1_complex_variant_4/distanceFiles/Dynein_1_complex_variant_4_cluster_2_ptm_distances.csv"
    ),
    C3 = file.path(
      CF_prefix,
      "Dynein_1_complex_variant_4/distanceFiles/Dynein_1_complex_variant_4_cluster_3_ptm_distances.csv"
    ),
    C4 = file.path(
      CF_prefix,
      "Dynein_1_complex_variant_4/distanceFiles/Dynein_1_complex_variant_4_cluster_4_ptm_distances.csv"
    )
  ),
  "eITF3" = list(
    C1 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_1_ptm_distances.csv"
    ),
    C2 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_2_ptm_distances.csv"
    ),
    C3 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_3_ptm_distances.csv"
    ),
    C4 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_4_ptm_distances.csv"
    ),
    C5 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_5_ptm_distances.csv"
    ),
    C6 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_6_ptm_distances.csv"
    ),
    C7 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_7_ptm_distances.csv"
    ),
    C8 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_8_ptm_distances.csv"
    ),
    C9 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_9_ptm_distances.csv"
    ),
    C10 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_10_ptm_distances.csv"
    )
  )
)

PF_files <- list(
  "26S-proteasome" = list(
    C1 = file.path(
      PF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_cluster_1_ptm_distances.csv"
    ),
    C2 = file.path(
      PF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_cluster_2_ptm_distances.csv"
    )
  ),
  "Brain-SWI-SNF" = list(
    C1 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_1_ptm_distances.csv"
    ),
    C2 = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_cluster_2_ptm_distances.csv"
    )
  ),
  "Dynein-1-complex-variant-1" = list(
    C2 = file.path(
      CF_prefix,
      "Dynein_1_complex_variant_4/distanceFiles/Dynein_1_complex_variant_4_cluster_2_ptm_distances.csv"
    )
  ),
  "eITF3" = list(
    C1 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_1_ptm_distances.csv"
    ),
    C2 = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_cluster_2_ptm_distances.csv"
    )
  )
)

# Check that all files exist
# for (method in list(CF_files, PF_files)) {
for (method in list(CF_files)) {
  for (complex in names(method)) {
    for (cluster in names(method[[complex]])) {
      if (cluster == "prefix") {
        next
      }
      file_path <- method[[complex]][[cluster]]
      if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
      }
    }
  }
}

################################################################################
## Create combined figures for each complex ####################################
################################################################################

for (complex in names(CF_files)) {
  cf_plots <- list()
  pf_plots <- list()

  for (i in names(CF_files[[complex]])) {
    file_path <- CF_files[[complex]][[i]]
    title <- paste0("Cluster ", i)
    cf_plots[[i]] <- plot_distance(file_path, title)
  }

  for (i in names(PF_files[[complex]])) {
    if (i == "prefix") {
      next
    }
    file_path <- PF_files[[complex]][[i]]
    title <- paste0("PF - Cluster ", i)
    pf_plots[[i]] <- plot_distance(file_path, title)
  }

  # Combine CF and save plots
  combined_cf <- cowplot::plot_grid(
    plotlist = cf_plots,
    ncol = length(cf_plots),
    align = "hv"
  )

  ggsave(
    filename = file.path(output_dir, paste0(complex, "_CF_Distances.pdf")),
    plot = combined_cf,
    width = 4 * length(cf_plots),
    height = 4.37,
    units = "in"
  )

  # # Combine PF and save plots
  # combined_pf <- cowplot::plot_grid(
  #   plotlist = pf_plots,
  #   ncol = length(pf_plots),
  #   align = "hv"
  # )

  # ggsave(
  #   filename = file.path(output_dir, paste0(complex, "_PF_Distances.pdf")),
  #   plot = combined_pf,
  #   width = 4 * length(pf_plots),
  #   height = 4.37,
  #   units = "in"
  # )
}
