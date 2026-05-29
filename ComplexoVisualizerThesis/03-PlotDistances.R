# Libraries
library(data.table)
library(ggplot2)

# Set working directory to project root
setwd(here::here())

output_dir <- "ComplexoVisualiserThesis/00-Plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

results_dir <- "ComplexoVisualiserThesis/00-Data/Results"

################################################################################
## Plotting function ###########################################################
################################################################################

plot_distance <- function(file_path, title) {
  data <- fread(file_path)
  data <- data[!is.na(Distance_Angstrom)]
  if (!"SameCoordFrame" %in% names(data)) {
    data[, SameCoordFrame := TRUE]
  }
  data <- data[SameCoordFrame == TRUE]
  data[, SameCoordFrame := NULL]

  make_label <- function(gene, residue, datatype) {
    paste(gene, residue, datatype, sep = "_")
  }

  data[, label1 := make_label(Gene1, Residue1, Datatype1)]
  data[, label2 := make_label(Gene2, Residue2, Datatype2)]

  peps1 <- unique(data[, .(label = label1, gene = Gene1, residue = Residue1)])
  peps2 <- unique(data[, .(label = label2, gene = Gene2, residue = Residue2)])
  all_peps <- unique(rbind(peps1, peps2))

  gene_order <- unique(c(data$Gene1, data$Gene2))
  all_peps[, gene := factor(gene, levels = gene_order)]
  setorder(all_peps, gene, residue)
  pep_levels <- all_peps$label

  mirrored <- data[, .(label1 = label2, label2 = label1, Distance_Angstrom)]
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

  data <- data[,
    .(Distance_Angstrom = min(Distance_Angstrom)),
    by = .(label1, label2)
  ]

  full_grid <- CJ(label1 = pep_levels, label2 = pep_levels)
  full_grid[, label1 := factor(label1, levels = pep_levels)]
  full_grid[, label2 := factor(label2, levels = pep_levels)]
  data <- merge(full_grid, data, by = c("label1", "label2"), all.x = TRUE)

  n_peps <- length(pep_levels)
  gene_sizes <- all_peps[, .N, by = gene][order(match(gene, gene_order))]
  x_breaks <- head(cumsum(gene_sizes$N), -1) + 0.5
  y_breaks <- n_peps - head(cumsum(gene_sizes$N), -1) + 0.5

  mid_idx <- round(
    cumsum(c(0, head(gene_sizes$N, -1))) + gene_sizes$N / 2 + 0.5
  )
  mid_idx <- pmin(pmax(mid_idx, 1L), n_peps)
  gene_labels <- all_peps$label[mid_idx]

  ggplot(data, aes(x = label1, y = label2, fill = Distance_Angstrom)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("#590007", "#932E05", "#C27142", "#DBAB8E", "#EBE5E0"),
      name    = "Distance (Å)",
      na.value = "white"
    ) +
    geom_vline(xintercept = x_breaks, colour = "black", linewidth = 0.6) +
    geom_hline(yintercept = y_breaks, colour = "black", linewidth = 0.6) +
    scale_x_discrete(
      breaks   = gene_labels,
      labels   = gene_order,
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
      axis.text.x  = element_text(angle = 45, hjust = 0, size = 7, vjust = 0),
      axis.text.y  = element_text(size = 7),
      panel.grid   = element_blank(),
      plot.title   = element_text(hjust = 0.5, size = 15),
      legend.position = "bottom"
    )
}

################################################################################
## Define input files ##########################################################
################################################################################

# Each entry: list of per-cluster distance CSVs.
# Cluster files are discovered automatically from the results directory so that
# adding new complexes or changing cluster counts requires no edits here.

discover_cluster_files <- function(complex_folder) {
  dist_dir <- file.path(results_dir, complex_folder, "distanceFiles")
  files <- list.files(
    dist_dir,
    pattern = "_cluster_[0-9]+_ptm_distances\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0L) return(NULL)
  # Sort by cluster number
  nums <- as.integer(gsub(".*_cluster_([0-9]+)_ptm_distances\\.csv$", "\\1", files))
  files <- files[order(nums)]
  setNames(as.list(files), paste0("C", sort(nums)))
}

# Map short label -> full complex name (must match cx$name in 02-VisualiseStructures.R)
# Folder names are derived with the same gsub used there, so they always stay in sync.
complex_names <- c(
  "26S-proteasome"             = "26S proteasome complex",
  "Brain-SWI-SNF"              = "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  "Dynein-1-complex-variant-1" = "Dynein-1 complex, variant 1",
  "eIF3"                       = "Eukaryotic translation initiation factor 3 complex"
)

complexes <- setNames(
  as.list(gsub("[^A-Za-z0-9]+", "_", complex_names)),
  names(complex_names)
)

################################################################################
## Plot ########################################################################
################################################################################

for (label in names(complexes)) {
  folder  <- complexes[[label]]
  cl_files <- discover_cluster_files(folder)

  if (is.null(cl_files)) {
    message("No cluster distance files found for ", label, " — skipping.")
    next
  }

  plots <- lapply(names(cl_files), function(cl) {
    plot_distance(cl_files[[cl]], paste0("Cluster ", sub("^C", "", cl)))
  })

  combined <- cowplot::plot_grid(
    plotlist = plots,
    ncol     = length(plots),
    align    = "hv"
  )

  out_path <- file.path(output_dir, paste0(label, "_Distances.pdf"))
  ggsave(
    filename = out_path,
    plot     = combined,
    width    = 4 * length(plots),
    height   = 4.37,
    units    = "in"
  )
  message("Written: ", out_path)
}
