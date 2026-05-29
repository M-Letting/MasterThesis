# Libraries
suppressMessages(library(arrow, warn.conflicts = FALSE))
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(ggplot2, warn.conflicts = FALSE))
suppressMessages(library(cowplot, warn.conflicts = FALSE))
suppressMessages(library(scico, warn.conflicts = FALSE))
suppressMessages(library(pcaMethods, warn.conflicts = FALSE))
suppressMessages(library(viridis, warn.conflicts = FALSE))

# Change working directory
setwd(here::here("Q-Benchmark/02-DCF"))

# Output directory
output_dir <- "./00-Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

##################
## PCA function ##
##################

plot_peptide_pca <- function(
  protein_data,
  intensity_cols = NULL,
  peptide_col = NULL,
  proteoform_col = "Proteoform_ID",
  scale_features = TRUE,
  dataset_type = c("proteomaker", "complex")
) {
  dataset_type <- match.arg(dataset_type)
  dt <- as.data.table(protein_data)

  if (is.null(intensity_cols)) {
    intensity_cols <- grep(
      "^C\\d+_R\\d+$|^B\\d+_D\\d+$",
      names(dt),
      value = TRUE
    )
  }
  if (!length(intensity_cols)) {
    stop("No intensity columns found.")
  }

  if (is.null(peptide_col)) {
    candidates <- c(
      "Peptidoform",
      "peptide_id",
      "Peptide",
      "Peptide.Modified.Sequence"
    )
    peptide_col <- candidates[candidates %in% names(dt)][1]
  }
  if (is.na(peptide_col)) {
    stop("No peptide id column found.")
  }

  ## Classification
  if (proteoform_col %in% names(dt)) {
    classify_pf <- function(x) {
      if (is.na(x) || x == "") {
        return("unknown")
      }
      vals <- unlist(strsplit(as.character(x), "\\|"))

      if (length(vals) == 1 && vals == "1") {
        "only_1"
      } else if ("1" %in% vals) {
        "contains_1"
      } else {
        "no_1"
      }
    }

    dt[, class := vapply(get(proteoform_col), classify_pf, character(1))]
  } else {
    dt[, class := "unknown"]
  }

  peptide_class <- unique(dt[, .(peptide = get(peptide_col), class)])

  ## PCA prep
  long <- melt(
    dt,
    id.vars = peptide_col,
    measure.vars = intensity_cols,
    variable.name = "Sample",
    value.name = "Intensity"
  )
  # Ensure correct naming of columns
  setnames(
    long,
    old = names(long),
    new = c("Peptidoform", "Sample", "Intensity")
  )
  setDT(long)

  if (dataset_type == "proteomaker") {
    long[, Condition := sub("_R\\d+$", "", Sample)]
  } else {
    long[, Condition := sub("^B\\d+_", "", Sample)]
  }

  long_avg <- long[,
    .(Intensity = mean(Intensity, na.rm = TRUE)),
    by = .(peptide = get(peptide_col), Condition)
  ]

  wide_avg <- dcast(long_avg, peptide ~ Condition, value.var = "Intensity")

  X <- as.matrix(wide_avg[, -"peptide"])
  rownames(X) <- wide_avg[["peptide"]]

  keep <- stats::complete.cases(X)
  X <- X[keep, , drop = FALSE]

  pca <- stats::prcomp(X, center = TRUE, scale. = scale_features)

  scores <- as.data.table(pca$x[, 1:2, drop = FALSE], keep.rownames = "peptide")
  scores <- merge(scores, peptide_class, by = "peptide", all.x = TRUE)

  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

  ## Colors
  pal <- scico::scico(3, palette = "batlow")
  class_levels <- c("only_1", "contains_1", "no_1", "unknown")

  class_colors <- c(
    setNames(pal, c("only_1", "contains_1", "no_1")),
    unknown = "#7f7f7f"
  )

  ## Dummy data to force full legend
  legend_dummy <- data.table(
    PC1 = NA,
    PC2 = NA,
    class = factor(class_levels, levels = class_levels)
  )

  ## Plot
  p <- ggplot(
    scores,
    aes(x = PC1, y = PC2, color = factor(class, levels = class_levels))
  ) +
    geom_point(size = 2.5, alpha = 0.95) +
    geom_point(
      data = legend_dummy,
      aes(x = PC1, y = PC2, color = class),
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = class_colors, drop = FALSE) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
      color = "Proteoform class"
    ) +
    coord_fixed() +
    theme_bw()

  return(list(plot = p))
}

#############################################
## Run PCA on top proteins across datasets ##
#############################################

ProteoMaker_datasets <- list(
  "LowNoNA" = "./00-Data/02-Prepared/ProteoMakerLowNoNA_wide_processed.csv",
  "MedNoNA" = "./00-Data/02-Prepared/ProteoMakerMedNoNA_wide_processed.csv",
  "HighNoNA" = "./00-Data/02-Prepared/ProteoMakerHighNoNA_wide_processed.csv"
)

pca_plots <- list()

# Top proteins
high_no_na_data <- fread(ProteoMaker_datasets[["HighNoNA"]])
top_proteins_high_no_na <- high_no_na_data[,
  .N,
  by = Accession
][order(-N)][1:5, Accession]

# Loop
for (name in names(ProteoMaker_datasets)) {
  data <- fread(ProteoMaker_datasets[[name]])
  intensity_cols <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)

  dataset_plot <- list()

  for (protein in top_proteins_high_no_na) {
    protein_data <- data[
      Accession == protein,
      c("Peptidoform", "Proteoform_ID", intensity_cols),
      with = FALSE
    ]

    res <- plot_peptide_pca(
      protein_data = protein_data,
      intensity_cols = intensity_cols,
      peptide_col = "Peptidoform",
      scale_features = FALSE
    )

    p <- res$plot +
      ggtitle(paste(name, "-", protein)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        aspect.ratio = 1
      )

    dataset_plot[[protein]] <- p
  }

  pca_plots[[name]] <- dataset_plot
}

#########################################
## Shared legend handling of PCA plots ##
#########################################

plot_list <- unlist(pca_plots, recursive = FALSE)

# Extract legend from first plot
legend <- get_legend(
  plot_list[[1]] + theme(legend.position = "bottom")
)

# Remove legends from all plots
plot_list <- lapply(plot_list, function(p) {
  p + theme(legend.position = "none")
})

# Combine plots
plots_grid <- plot_grid(
  plotlist = plot_list,
  ncol = 5
)

# Final combined plot with shared legend
combined_plot <- plot_grid(
  plots_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Add title to the combined plot
title <- ggdraw() +
  draw_label(
    "PCA of top proteins across ProteoMaker datasets",
    fontface = "bold",
    size = 20
  )

combined_plot <- plot_grid(
  title,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

# Save the combined plot
ggsave(
  filename = file.path(output_dir, "07-ProteoMakerTopProteinsGrouping.pdf"),
  plot = combined_plot,
  width = 12,
  height = 9,
  dpi = 300
)

################################################################################
## PCA plots with coloring based on grouping by methods ########################
################################################################################

plot_peptide_pca_method <- function(
  protein_data,
  intensity_cols = NULL,
  peptide_col = NULL,
  proteoform_col = "Proteoform_ID",
  scale_features = TRUE,
  dataset_type = c("proteomaker", "complex")
) {
  dataset_type <- match.arg(dataset_type)
  dt <- as.data.table(protein_data)

  # Intensity columns
  if (is.null(intensity_cols)) {
    intensity_cols <- grep(
      "^C\\d+_R\\d+$|^B\\d+_D\\d+$",
      names(dt),
      value = TRUE
    )
  }
  if (!length(intensity_cols)) {
    stop("No intensity columns found.")
  }

  # Peptide column
  if (is.null(peptide_col)) {
    candidates <- c(
      "Peptidoform",
      "peptide_id",
      "Peptide",
      "Peptide.Modified.Sequence"
    )
    peptide_col <- candidates[candidates %in% names(dt)][1]
  }
  if (is.na(peptide_col)) {
    stop("No peptide id column found.")
  }

  # Determined class based on proteoform values
  if (!proteoform_col %in% names(dt)) {
    dt[, class := "unknown"]
  } else {
    dt[, class := as.character(get(proteoform_col))]
    dt[is.na(class) | class == "", class := "unknown"]
  }

  peptide_class <- unique(dt[, .(peptide = get(peptide_col), class)])

  # PCA prep
  long <- melt(
    dt,
    id.vars = peptide_col,
    measure.vars = intensity_cols,
    variable.name = "Sample",
    value.name = "Intensity"
  )
  # Ensure correct naming of columns
  setnames(long, old = names(long), new = c(peptide_col, "Sample", "Intensity"))
  setDT(long)

  if (dataset_type == "proteomaker") {
    long[, Condition := sub("_R\\d+$", "", Sample)]
  } else {
    long[, Condition := sub("^B\\d+_", "", Sample)]
  }

  long_avg <- long[,
    .(Intensity = mean(Intensity, na.rm = TRUE)),
    by = .(peptide = get(peptide_col), Condition)
  ]

  wide_avg <- dcast(long_avg, peptide ~ Condition, value.var = "Intensity")

  X <- as.matrix(wide_avg[, -"peptide"])
  rownames(X) <- wide_avg[["peptide"]]

  keep <- stats::complete.cases(X)
  X <- X[keep, , drop = FALSE]

  pca <- stats::prcomp(X, center = TRUE, scale. = scale_features)

  scores <- as.data.table(pca$x[, 1:2, drop = FALSE], keep.rownames = "peptide")
  scores <- merge(scores, peptide_class, by = "peptide", all.x = TRUE)

  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

  # Determine levels and grouping (stable ordering)
  class_levels <- sort(unique(scores$class))

  n_colors <- length(class_levels)

  # Generate palette dynamically
  pal <- scico::scico(n_colors, palette = "batlow")

  class_colors <- setNames(pal, class_levels)

  # Dummy data to force full legend
  legend_dummy <- data.table(
    PC1 = NA,
    PC2 = NA,
    class = factor(class_levels, levels = class_levels)
  )

  # Plot
  p <- ggplot(
    scores,
    aes(x = PC1, y = PC2, color = factor(class, levels = class_levels))
  ) +
    geom_point(size = 2.5, alpha = 0.95) +
    geom_point(
      data = legend_dummy,
      aes(x = PC1, y = PC2, color = class),
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = class_colors, drop = FALSE) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
      color = proteoform_col
    ) +
    coord_fixed() +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.title = element_blank()
    )

  return(list(plot = p))
}

##########################
## ProteoForge datasets ##
##########################

# ProteoForge_results <- list(
#   "LowNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNoNA_result.feather",
#   "MedNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNoNA_result.feather",
#   "HighNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNoNA_result.feather"
# )

# ProteoForge_NoNA_PCA <- list()

# for (name in names(ProteoForge_results)) {
#   if (!file.exists(ProteoForge_results[[name]])) {
#     warning(paste("File not found:", ProteoForge_results[[name]]))
#     next
#   }

#   # Get ProteoForge grouping results
#   PF_res <- as.data.table(read_feather(ProteoForge_results[[name]]))

#   # Subset ProteoForge results to only Sample = "C1_R1"
#   PF_res <- PF_res[Sample == "C1_R1"]

#   # Subset to protein_id, peptide_id, and ClusterID
#   PF_res <- PF_res[, .(peptide_id, protein_id, ClusterID)]

#   # Get original ProteoMaker data for this dataset to merge
#   protein_data <- fread(ProteoMaker_datasets[[name]])
#   intensity_cols <- grep("^C\\d+_R\\d+$", names(protein_data), value = TRUE)
#   protein_data <- protein_data[,
#     c("Accession", "Peptidoform", intensity_cols),
#     with = FALSE
#   ]

#   # Subset to top_proteins_high_no_na
#   protein_data <- protein_data[Accession %in% top_proteins_high_no_na]

#   # Merge ProteoForge results with original data based on peptide_id and protein_id
#   merged_data <- merge(
#     PF_res,
#     protein_data,
#     by.x = c("peptide_id", "protein_id"),
#     by.y = c("Peptidoform", "Accession")
#   )

#   # Rename ClusterID to Proteoform_ID for PCA function compatibility
#   setnames(merged_data, "ClusterID", "Proteoform_ID")

#   dataset_plots <- list()

#   for (protein in top_proteins_high_no_na) {
#     protein_data <- merged_data[
#       protein_id == protein,
#       c("peptide_id", "Proteoform_ID", intensity_cols),
#       with = FALSE
#     ]

#     res <- plot_peptide_pca_method(
#       protein_data = protein_data,
#       intensity_cols = intensity_cols,
#       peptide_col = "peptide_id",
#       proteoform_col = "Proteoform_ID",
#       scale_features = FALSE
#     )

#     p <- res$plot +
#       ggtitle(paste(name, "-", protein)) +
#       theme(
#         plot.title = element_text(hjust = 0.5, face = "bold"),
#         aspect.ratio = 1
#       )

#     dataset_plots[[protein]] <- p
#   }

#   ProteoForge_NoNA_PCA[[name]] <- dataset_plots
# }

# # Plot in grid
# plot_list <- unlist(ProteoForge_NoNA_PCA, recursive = FALSE)

# combined_plot1 <- plot_grid(
#   plotlist = plot_list,
#   ncol = 5
# )

# combined_plot1 <- plot_grid(
#   ggdraw() +
#     draw_label(
#       "PCA of top proteins across ProteoMaker datasets colored by ProteoForge NoNA grouping",
#       fontface = "bold",
#       size = 20
#     ),
#   combined_plot1,
#   ncol = 1,
#   rel_heights = c(0.08, 1)
# )

# combined_plot1

# Repeat with NA grouping from ProteoForge results
ProteoForge_results <- list(
  "LowNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNA_result.feather",
  "MedNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNA_result.feather",
  "HighNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNA_result.feather"
)

ProteoForge_NA_PCA <- list()

for (name in names(ProteoForge_results)) {
  if (!file.exists(ProteoForge_results[[name]])) {
    warning(paste("File not found:", ProteoForge_results[[name]]))
    next
  }

  # Get ProteoForge grouping results
  PF_res <- as.data.table(read_feather(ProteoForge_results[[name]]))

  # Subset ProteoForge results to only Sample = "C1_R1"
  PF_res <- PF_res[Sample == "C1_R1"]

  # Subset to protein_id, peptide_id, and ClusterID
  PF_res <- PF_res[, .(peptide_id, protein_id, ClusterID)]

  # Get original ProteoMaker data for this dataset to merge
  alt_name <- sub("NA$", "NoNA", name)
  protein_data <- fread(ProteoMaker_datasets[[alt_name]])
  intensity_cols <- grep("^C\\d+_R\\d+$", names(protein_data), value = TRUE)
  protein_data <- protein_data[,
    c("Accession", "Peptidoform", intensity_cols),
    with = FALSE
  ]

  # Subset to top_proteins_high_no_na
  protein_data <- protein_data[Accession %in% top_proteins_high_no_na]

  # Merge ProteoForge results with original data based on peptide_id and protein_id
  merged_data <- merge(
    PF_res,
    protein_data,
    by.x = c("peptide_id", "protein_id"),
    by.y = c("Peptidoform", "Accession")
  )

  # Rename ClusterID to Proteoform_ID for PCA function compatibility
  setnames(merged_data, "ClusterID", "Proteoform_ID")

  dataset_plots <- list()

  for (protein in top_proteins_high_no_na) {
    protein_data <- merged_data[
      protein_id == protein,
      c("peptide_id", "Proteoform_ID", intensity_cols),
      with = FALSE
    ]

    res <- plot_peptide_pca_method(
      protein_data = protein_data,
      intensity_cols = intensity_cols,
      peptide_col = "peptide_id",
      proteoform_col = "Proteoform_ID",
      scale_features = FALSE
    )

    p <- res$plot +
      ggtitle(paste(name, "-", protein)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        aspect.ratio = 1
      )
    dataset_plots[[protein]] <- p
  }

  ProteoForge_NA_PCA[[name]] <- dataset_plots
}

plot_list <- unlist(ProteoForge_NA_PCA, recursive = FALSE)

combined_plot2 <- plot_grid(
  plotlist = plot_list,
  ncol = 5
)

combined_plot2 <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of top proteins across ProteoMaker datasets colored by ProteoForge NA grouping",
      fontface = "bold",
      size = 20
    ),
  combined_plot2,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

# Save the combined plot
ggsave(
  filename = file.path(
    output_dir,
    "07-ProteoMakerTopProteinsProteoForgeGrouping-NA.pdf"
  ),
  plot = combined_plot2,
  width = 12,
  height = 9,
  dpi = 300
)

#########
## DCF ##
#########

DCF_NA_results <- list(
  "LowNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerLowNA_results.csv",
  "MedNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerMedNA_results.csv",
  "HighNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerHighNA_results.csv"
)

DCF_NA_PCA <- list()

for (name in names(DCF_NA_results)) {
  if (!file.exists(DCF_NA_results[[name]])) {
    warning(paste("File not found:", DCF_NA_results[[name]]))
    next
  }

  dcf_res <- fread(DCF_NA_results[[name]])
  dcf_res[, Peptidoform := sub(".*_(?=[^_]+$)", "", Peptide_ID, perl = TRUE)]
  dcf_res <- dcf_res[, .(Accession, Peptidoform, dCF)]

  alt_name <- sub("NA$", "NoNA", name)
  protein_data <- fread(ProteoMaker_datasets[[alt_name]])
  intensity_cols <- grep("^C\\d+_R\\d+$", names(protein_data), value = TRUE)
  protein_data <- protein_data[,
    c("Accession", "Peptidoform", intensity_cols),
    with = FALSE
  ]
  protein_data <- protein_data[Accession %in% top_proteins_high_no_na]

  merged_data <- merge(
    dcf_res,
    protein_data,
    by = c("Peptidoform", "Accession")
  )

  dataset_plots <- list()

  for (protein in top_proteins_high_no_na) {
    pdata <- merged_data[
      Accession == protein,
      c("Peptidoform", "dCF", intensity_cols),
      with = FALSE
    ]

    res <- plot_peptide_pca_method(
      protein_data = pdata,
      intensity_cols = intensity_cols,
      peptide_col = "Peptidoform",
      proteoform_col = "dCF",
      scale_features = FALSE
    )

    p <- res$plot +
      ggtitle(paste(name, "-", protein)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        aspect.ratio = 1
      )
    dataset_plots[[protein]] <- p
  }

  DCF_NA_PCA[[name]] <- dataset_plots
}

plot_list <- unlist(DCF_NA_PCA, recursive = FALSE)

combined_DCF_plot <- plot_grid(plotlist = plot_list, ncol = 5)

combined_DCF_plot <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of top proteins across ProteoMaker datasets colored by DCF NA grouping",
      fontface = "bold",
      size = 20
    ),
  combined_DCF_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

ggsave(
  filename = file.path(
    output_dir,
    "07-ProteoMakerTopProteinsDCFGrouping-NA.pdf"
  ),
  plot = combined_DCF_plot,
  width = 12,
  height = 9,
  dpi = 300
)

##################
## DCFBaseline  ##
##################

DCFBaseline_NA_results <- list(
  "LowNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNA_results.csv",
  "MedNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNA_results.csv",
  "HighNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNA_results.csv"
)

DCFBaseline_NA_PCA <- list()

for (name in names(DCFBaseline_NA_results)) {
  if (!file.exists(DCFBaseline_NA_results[[name]])) {
    warning(paste("File not found:", DCFBaseline_NA_results[[name]]))
    next
  }

  dcf_res <- fread(DCFBaseline_NA_results[[name]])
  dcf_res[, Peptidoform := sub(".*_(?=[^_]+$)", "", Peptide_ID, perl = TRUE)]
  dcf_res <- dcf_res[, .(Accession, Peptidoform, dCF)]

  alt_name <- sub("NA$", "NoNA", name)
  protein_data <- fread(ProteoMaker_datasets[[alt_name]])
  intensity_cols <- grep("^C\\d+_R\\d+$", names(protein_data), value = TRUE)
  protein_data <- protein_data[,
    c("Accession", "Peptidoform", intensity_cols),
    with = FALSE
  ]
  protein_data <- protein_data[Accession %in% top_proteins_high_no_na]

  merged_data <- merge(
    dcf_res,
    protein_data,
    by = c("Peptidoform", "Accession")
  )

  dataset_plots <- list()

  for (protein in top_proteins_high_no_na) {
    pdata <- merged_data[
      Accession == protein,
      c("Peptidoform", "dCF", intensity_cols),
      with = FALSE
    ]

    res <- plot_peptide_pca_method(
      protein_data = pdata,
      intensity_cols = intensity_cols,
      peptide_col = "Peptidoform",
      proteoform_col = "dCF",
      scale_features = FALSE
    )

    p <- res$plot +
      ggtitle(paste(name, "-", protein)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        aspect.ratio = 1
      )
    dataset_plots[[protein]] <- p
  }

  DCFBaseline_NA_PCA[[name]] <- dataset_plots
}

plot_list <- unlist(DCFBaseline_NA_PCA, recursive = FALSE)

combined_DCFBaseline_plot <- plot_grid(plotlist = plot_list, ncol = 5)

combined_DCFBaseline_plot <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of top proteins across ProteoMaker datasets colored by DCFBaseline NA grouping",
      fontface = "bold",
      size = 20
    ),
  combined_DCFBaseline_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

ggsave(
  filename = file.path(
    output_dir,
    "07-ProteoMakerTopProteinsDCFBaselineGrouping-NA.pdf"
  ),
  plot = combined_DCFBaseline_plot,
  width = 12,
  height = 9,
  dpi = 300
)

#########
## RPC ##
#########

RPC_results <- list(
  "LowNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerLowNA_results.csv",
  "MedNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerMedNA_results.csv",
  "HighNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerHighNA_results.csv"
)

RPC_NA_PCA <- list()

for (name in names(RPC_results)) {
  if (!file.exists(RPC_results[[name]])) {
    warning(paste("File not found:", RPC_results[[name]]))
    next
  }

  RPC_res <- fread(RPC_results[[name]])
  RPC_res <- RPC_res[, .(Accession, Peptide_ID, dCF)]
  RPC_res[, Peptide_ID := sub(".*_", "", Peptide_ID)]

  alt_name <- sub("NA$", "NoNA", name)
  protein_data <- fread(ProteoMaker_datasets[[alt_name]])
  intensity_cols <- grep("^C\\d+_R\\d+$", names(protein_data), value = TRUE)
  protein_data <- protein_data[,
    c("Accession", "Peptidoform", intensity_cols),
    with = FALSE
  ]
  protein_data <- protein_data[Accession %in% top_proteins_high_no_na]
  merged_data <- merge(
    RPC_res,
    protein_data,
    by.x = c("Peptide_ID", "Accession"),
    by.y = c("Peptidoform", "Accession")
  )
  dataset_plots <- list()
  for (protein in top_proteins_high_no_na) {
    protein_data <- merged_data[
      Accession == protein,
      c("Peptide_ID", "dCF", intensity_cols),
      with = FALSE
    ]
    res <- plot_peptide_pca_method(
      protein_data = protein_data,
      intensity_cols = intensity_cols,
      peptide_col = "Peptide_ID",
      proteoform_col = "dCF",
      scale_features = FALSE
    )
    p <- res$plot +
      ggtitle(paste(name, "-", protein)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        aspect.ratio = 1
      )
    dataset_plots[[protein]] <- p
  }
  RPC_NA_PCA[[name]] <- dataset_plots
}

plot_list <- unlist(RPC_NA_PCA, recursive = FALSE)
combined_plot2_1 <- plot_grid(plotlist = plot_list, ncol = 5)
combined_plot2_1 <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of top proteins across ProteoMaker datasets colored by RPC NA grouping",
      fontface = "bold",
      size = 20
    ),
  combined_plot2_1,
  ncol = 1,
  rel_heights = c(0.08, 1)
)
ggsave(
  filename = file.path(
    output_dir,
    "07-ProteoMakerTopProteinsRPCGrouping-NA.pdf"
  ),
  plot = combined_plot2_1,
  width = 12,
  height = 9,
  dpi = 300
)

################################################################################
## Compare mapping of protein complex data #####################################
################################################################################

# Datasets for complex data
complex_datasets <- list(
  "B-WICH" = "./00-Data/02-Prepared/B-Wich_WideLogNorm.csv",
  "CRD" = "./00-Data/02-Prepared/CRD-mediated_WideLogNorm.csv",
  "Multiaminoacyl" = "./00-Data/02-Prepared/Multiaminoacyl-tRNA_WideLogNorm.csv",
  "26S-proteasome" = "./00-Data/02-Prepared/26S-proteasome_WideLogNorm.csv",
  "60S-ribosome" = "./00-Data/02-Prepared/60S-cytosolic-large-ribosomal-subunit_WideLogNorm.csv",
  "BSSADCRC" = "./00-Data/02-Prepared/BSSADCRC_WideLogNorm.csv",
  "Dynactin" = "./00-Data/02-Prepared/Dynactin-complex_WideLogNorm.csv",
  "Dynein-1" = "./00-Data/02-Prepared/Dynein-1-complex-variant-1_WideLogNorm.csv",
  "eIF3" = "./00-Data/02-Prepared/Eukaryotic-translation-initiation-factor-3-complex_WideLogNorm.csv",
  "IFT-B" = "./00-Data/02-Prepared/Intraflagellar-transport-complex-B_WideLogNorm.csv",
  "Laminin-213" = "./00-Data/02-Prepared/Laminin-213_WideLogNorm.csv",
  "Spliceosomal-B" = "./00-Data/02-Prepared/Major-Spliceosomal-B_WideLogNorm.csv",
  "Nuclear-pore" = "./00-Data/02-Prepared/Nuclear-pore-complex_WideLogNorm.csv"
)

# Complex PCA plots
complex_PCA_plots <- list()

for (name in names(complex_datasets)) {
  if (!file.exists(complex_datasets[[name]])) {
    warning(paste("File not found:", complex_datasets[[name]]))
    next
  }

  # Get RPC grouping results for complex data
  data <- fread(complex_datasets[[name]])

  intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)

  # Subset to "Peptide" and intensity columns
  data <- data[, c("Peptide", intensity_cols), with = FALSE]

  # Remove rows with all NA values in intensity columns
  data <- data[rowSums(!is.na(data[, ..intensity_cols])) > 0]

  # Calculate average for each condition for each peptide
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
  wide_avg <- dcast(long_avg, Peptide ~ Condition, value.var = "Intensity")

  # Remove Peptide column
  X <- as.data.frame((wide_avg))
  rownames(X) <- X$Peptide
  X$Peptide <- NULL

  # Set NaN values to NA
  for (col in names(X)) {
    X[[col]][is.nan(X[[col]])] <- NA
  }

  # Remove rows with > 50% NA values
  na_threshold <- ncol(X) * 1.0
  X <- as.data.frame(X[rowSums(is.na(X)) <= na_threshold, ])

  pca <- pcaMethods::pca(
    X,
    method = "ppca",
    nPcs = 2,
    scale = "none",
    center = TRUE
  )
  scores <- as.data.frame(pca@scores)
  scores$Acc <- rownames(scores)
  scores$Acc <- sub("^[^_]+_([^_]+)_.*$", "\\1", scores$Acc)
  var_expl <- pca@R2

  # Colors
  pal <- scico::scico(length(unique(scores$Acc)), palette = "batlow")
  acc_levels <- sort(unique(scores$Acc))
  acc_colors <- setNames(pal, acc_levels)

  # Plot
  p <- ggplot(scores, aes(x = PC1, y = PC2, color = Acc)) +
    geom_point(size = 2.5, alpha = 0.95) +
    scale_color_manual(values = acc_colors) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2])
    ) +
    coord_fixed() +
    theme_bw() +
    ggtitle(paste(name, "- Complex Data")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.title = element_blank()
    )
  complex_PCA_plots[[name]] <- p
}

# Combine complex PCA plots
combined_complex_plot <- plot_grid(
  plotlist = complex_PCA_plots,
  ncol = 5
)

combined_complex_plot <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of complex datasets",
      fontface = "bold",
      size = 20
    ),
  combined_complex_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

# Save the combined plot
ggsave(
  filename = file.path(output_dir, "07-ComplexDataPCA.pdf"),
  plot = combined_complex_plot,
  width = 20,
  height = 12,
  dpi = 300
)

########################################
# ProteoForge mapping for complex data #
########################################

# ProteoForge results for complex data
ProteoForge_complex_results <- list(
  "B-WICH" = "./00-Data/03-Results/ProteoForge/ProteoForge_B-Wich_result.feather",
  "CRD" = "./00-Data/03-Results/ProteoForge/ProteoForge_CRD-mediated_result.feather",
  "Multiaminoacyl" = "./00-Data/03-Results/ProteoForge/ProteoForge_Multiaminoacyl-tRNA_result.feather",
  "26S-proteasome" = "./00-Data/03-Results/ProteoForge/ProteoForge_26S-proteasome_result.feather",
  "60S-ribosome" = "./00-Data/03-Results/ProteoForge/ProteoForge_60S-cytosolic-large-ribosomal-subunit_result.feather",
  "BSSADCRC" = "./00-Data/03-Results/ProteoForge/ProteoForge_BSSADCRC_result.feather",
  "Dynactin" = "./00-Data/03-Results/ProteoForge/ProteoForge_Dynactin-complex_result.feather",
  "Dynein-1" = "./00-Data/03-Results/ProteoForge/ProteoForge_Dynein-1-complex-variant-1_result.feather",
  "eIF3" = "./00-Data/03-Results/ProteoForge/ProteoForge_Eukaryotic-translation-initiation-factor-3-complex_result.feather",
  "IFT-B" = "./00-Data/03-Results/ProteoForge/ProteoForge_Intraflagellar-transport-complex-B_result.feather",
  "Laminin-213" = "./00-Data/03-Results/ProteoForge/ProteoForge_Laminin-213_result.feather",
  "Spliceosomal-B" = "./00-Data/03-Results/ProteoForge/ProteoForge_Major-Spliceosomal-B_result.feather",
  "Nuclear-pore" = "./00-Data/03-Results/ProteoForge/ProteoForge_Nuclear-pore-complex_result.feather"
)

ProteoForge_complex_PCA <- list()

for (name in names(ProteoForge_complex_results)) {
  if (!file.exists(ProteoForge_complex_results[[name]])) {
    warning(paste("File not found:", ProteoForge_complex_results[[name]]))
    next
  }

  # Get ProteoForge grouping results
  PF_res <- as.data.table(read_feather(ProteoForge_complex_results[[name]]))

  # Convert to wide format
  PF_res <- dcast(
    PF_res,
    peptide_id + ClusterID ~ Sample,
    value.var = "ms1adj",
    fun.aggregate = mean
  )

  # Combine peptide_id and ClusterID for merging
  PF_res[, peptide_id := paste0(peptide_id, "_", ClusterID)]

  data <- PF_res

  # Get intensity columns
  intensity_cols <- grep("^B\\d+_D\\d+$", names(data), value = TRUE)

  # Calculate average for each condition for each peptide
  long <- data.table::melt(
    data,
    id.vars = "peptide_id",
    measure.vars = intensity_cols,
    variable.name = "Sample",
    value.name = "Intensity"
  )

  long[, Condition := sub("^B\\d+_", "", Sample)]
  long_avg <- long[,
    .(Intensity = mean(Intensity, na.rm = TRUE)),
    by = .(peptide_id, Condition)
  ]
  data <- dcast(long_avg, peptide_id ~ Condition, value.var = "Intensity")

  # Remove all rows in PF_res with duplicate peptide_id
  data <- data[!duplicated(peptide_id)]

  data <- as.data.frame(data)

  rownames(data) <- data$peptide_id
  data$peptide_id <- NULL

  # Set NaN values to NA
  for (col in names(data)) {
    data[[col]][is.nan(data[[col]])] <- NA
  }

  # Run PCA
  pca <- pcaMethods::pca(
    data,
    method = "ppca",
    nPcs = 2,
    scale = "none",
    center = TRUE
  )
  scores <- as.data.frame(pca@scores)
  scores$ID <- rownames(scores)
  # Set ID to just ClusterID (value after last "_")
  scores$ID <- sub(".*_", "", scores$ID)

  # Set cluster ID of all clusters with 1 or fewer peptides to "singleton"
  cluster_counts <- table(scores$ID)
  singleton_clusters <- names(cluster_counts[cluster_counts <= 1])
  scores$ID[scores$ID %in% singleton_clusters] <- "singleton"

  var_expl <- pca@R2

  # Colors
  pal <- scico::scico(length(unique(scores$ID)), palette = "batlow")
  id_levels <- sort(unique(scores$ID))
  id_colors <- setNames(pal, id_levels)

  # Plot
  p <- ggplot(scores, aes(x = PC1, y = PC2, color = ID)) +
    geom_point(size = 2.5, alpha = 0.95) +
    scale_color_manual(values = id_colors) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
      color = "ProteoForge ClusterID"
    ) +
    coord_fixed() +
    theme_bw() +
    ggtitle(paste(name, "- ProteoForge Clusters")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.title = element_blank()
    )
  ProteoForge_complex_PCA[[name]] <- p
}

# Combine ProteoForge complex PCA plots
combined_PF_complex_plot <- plot_grid(
  plotlist = ProteoForge_complex_PCA,
  ncol = 5
)
combined_PF_complex_plot <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of complex datasets colored by ProteoForge clusters",
      fontface = "bold",
      size = 20
    ),
  combined_PF_complex_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

# Save the combined plot
ggsave(
  filename = file.path(output_dir, "07-ComplexDataProteoForgePCA.pdf"),
  plot = combined_PF_complex_plot,
  width = 20,
  height = 12,
  dpi = 300
)

##############################
# DCF mapping for complex data
##############################

DCF_complex_results <- list(
  "B-WICH" = "./00-Data/03-Results/DCF/DCF_Complex_B-Wich_results.csv",
  "CRD" = "./00-Data/03-Results/DCF/DCF_Complex_CRD-mediated_results.csv",
  "Multiaminoacyl" = "./00-Data/03-Results/DCF/DCF_Complex_Multiaminoacyl-tRNA_results.csv",
  "26S-proteasome" = "./00-Data/03-Results/DCF/DCF_Complex_26S-proteasome_results.csv",
  "60S-ribosome" = "./00-Data/03-Results/DCF/DCF_Complex_60S-cytosolic-large-ribosomal-subunit_results.csv",
  "BSSADCRC" = "./00-Data/03-Results/DCF/DCF_Complex_BSSADCRC_results.csv",
  "Dynactin" = "./00-Data/03-Results/DCF/DCF_Complex_Dynactin-complex_results.csv",
  "Dynein-1" = "./00-Data/03-Results/DCF/DCF_Complex_Dynein-1-complex-variant-1_results.csv",
  "eIF3" = "./00-Data/03-Results/DCF/DCF_Complex_Eukaryotic-translation-initiation-factor-3-complex_results.csv",
  "IFT-B" = "./00-Data/03-Results/DCF/DCF_Complex_Intraflagellar-transport-complex-B_results.csv",
  "Laminin-213" = "./00-Data/03-Results/DCF/DCF_Complex_Laminin-213_results.csv",
  "Spliceosomal-B" = "./00-Data/03-Results/DCF/DCF_Complex_Major-Spliceosomal-B_results.csv",
  "Nuclear-pore" = "./00-Data/03-Results/DCF/DCF_Complex_Nuclear-pore-complex_results.csv"
)

DCF_complex_PCA <- list()

for (name in names(DCF_complex_results)) {
  if (!file.exists(DCF_complex_results[[name]])) {
    warning(paste("File not found:", DCF_complex_results[[name]]))
    next
  }

  dcf_res <- fread(DCF_complex_results[[name]])

  # Join with original complex data to get intensities
  orig <- fread(complex_datasets[[name]])
  intensity_cols <- grep("^B\\d+_D\\d+$", names(orig), value = TRUE)

  merged <- merge(
    dcf_res[, .(Peptide, dCF)],
    orig[, c("Peptide", intensity_cols), with = FALSE],
    by = "Peptide"
  )

  long <- data.table::melt(
    merged,
    id.vars = c("Peptide", "dCF"),
    measure.vars = intensity_cols,
    variable.name = "Sample",
    value.name = "Intensity"
  )
  long[, Condition := sub("^B\\d+_", "", Sample)]
  long_avg <- long[,
    .(Intensity = mean(Intensity, na.rm = TRUE)),
    by = .(Peptide, dCF, Condition)
  ]
  wide_avg <- dcast(
    long_avg,
    Peptide + dCF ~ Condition,
    value.var = "Intensity"
  )

  cond_cols <- setdiff(names(wide_avg), c("Peptide", "dCF"))
  mat <- as.data.frame(wide_avg[, cond_cols, with = FALSE])
  rownames(mat) <- paste0(wide_avg$Peptide, "_", wide_avg$dCF)
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
  scores <- as.data.frame(pca@scores)
  scores$ID <- sub(".*_", "", rownames(scores))
  var_expl <- pca@R2

  pal <- scico::scico(length(unique(scores$ID)), palette = "batlow")
  id_levels <- sort(unique(scores$ID))
  id_colors <- setNames(pal, id_levels)

  p <- ggplot(scores, aes(x = PC1, y = PC2, color = ID)) +
    geom_point(size = 2.5, alpha = 0.95) +
    scale_color_manual(values = id_colors) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
      color = "DCF dCF"
    ) +
    coord_fixed() +
    theme_bw() +
    ggtitle(paste(name, "- DCF Clusters")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.title = element_blank()
    )
  DCF_complex_PCA[[name]] <- p
}

combined_DCF_complex_plot <- plot_grid(
  plotlist = Filter(Negate(is.null), DCF_complex_PCA),
  ncol = 5
)
combined_DCF_complex_plot <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of complex datasets colored by DCF clusters",
      fontface = "bold",
      size = 20
    ),
  combined_DCF_complex_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

ggsave(
  filename = file.path(output_dir, "07-ComplexDataDCFPCA.pdf"),
  plot = combined_DCF_complex_plot,
  width = 20,
  height = 12,
  dpi = 300
)

########################################
# DCFBaseline mapping for complex data #
########################################

DCFBaseline_complex_results <- list(
  "B-WICH" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_B-Wich_results.csv",
  "CRD" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_CRD-mediated_results.csv",
  "Multiaminoacyl" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_Multiaminoacyl-tRNA_results.csv",
  "26S-proteasome" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_26S-proteasome_results.csv",
  "60S-ribosome" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_60S-cytosolic-large-ribosomal-subunit_results.csv",
  "BSSADCRC" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_BSSADCRC_results.csv",
  "Dynactin" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_Dynactin-complex_results.csv",
  "Dynein-1" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_Dynein-1-complex-variant-1_results.csv",
  "eIF3" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_Eukaryotic-translation-initiation-factor-3-complex_results.csv",
  "IFT-B" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_Intraflagellar-transport-complex-B_results.csv",
  "Laminin-213" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_Laminin-213_results.csv",
  "Spliceosomal-B" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_Major-Spliceosomal-B_results.csv",
  "Nuclear-pore" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_Complex_Nuclear-pore-complex_results.csv"
)

DCFBaseline_complex_PCA <- list()

for (name in names(DCFBaseline_complex_results)) {
  if (!file.exists(DCFBaseline_complex_results[[name]])) {
    warning(paste("File not found:", DCFBaseline_complex_results[[name]]))
    next
  }

  dcf_res <- fread(DCFBaseline_complex_results[[name]])

  # Join with original complex data to get intensities
  orig <- fread(complex_datasets[[name]])
  intensity_cols <- grep("^B\\d+_D\\d+$", names(orig), value = TRUE)

  merged <- merge(
    dcf_res[, .(Peptide, dCF)],
    orig[, c("Peptide", intensity_cols), with = FALSE],
    by = "Peptide"
  )

  long <- data.table::melt(
    merged,
    id.vars = c("Peptide", "dCF"),
    measure.vars = intensity_cols,
    variable.name = "Sample",
    value.name = "Intensity"
  )
  long[, Condition := sub("^B\\d+_", "", Sample)]
  long_avg <- long[,
    .(Intensity = mean(Intensity, na.rm = TRUE)),
    by = .(Peptide, dCF, Condition)
  ]
  wide_avg <- dcast(
    long_avg,
    Peptide + dCF ~ Condition,
    value.var = "Intensity"
  )

  cond_cols <- setdiff(names(wide_avg), c("Peptide", "dCF"))
  mat <- as.data.frame(wide_avg[, cond_cols, with = FALSE])
  rownames(mat) <- paste0(wide_avg$Peptide, "_", wide_avg$dCF)
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
  scores <- as.data.frame(pca@scores)
  scores$ID <- sub(".*_", "", rownames(scores))
  var_expl <- pca@R2

  pal <- scico::scico(length(unique(scores$ID)), palette = "batlow")
  id_levels <- sort(unique(scores$ID))
  id_colors <- setNames(pal, id_levels)

  p <- ggplot(scores, aes(x = PC1, y = PC2, color = ID)) +
    geom_point(size = 2.5, alpha = 0.95) +
    scale_color_manual(values = id_colors) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
      color = "DCFBaseline dCF"
    ) +
    coord_fixed() +
    theme_bw() +
    ggtitle(paste(name, "- DCFBaseline Clusters")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.title = element_blank()
    )
  DCFBaseline_complex_PCA[[name]] <- p
}

combined_DCFBaseline_complex_plot <- plot_grid(
  plotlist = DCFBaseline_complex_PCA,
  ncol = 5
)
combined_DCFBaseline_complex_plot <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of complex datasets colored by DCFBaseline clusters",
      fontface = "bold",
      size = 20
    ),
  combined_DCFBaseline_complex_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

ggsave(
  filename = file.path(output_dir, "07-ComplexDataDCFBaselinePCA.pdf"),
  plot = combined_DCFBaseline_complex_plot,
  width = 20,
  height = 12,
  dpi = 300
)

# ###############################
# # COPF mapping for complex data
# ###############################

# COPF_complex_results <- list(
#   "B-WICH" = "./00-Data/03-Results/COPF/COPF_Complex_B-Wich_results.csv",
#   "CRD" = "./00-Data/03-Results/COPF/COPF_Complex_CRD-mediated_results.csv",
#   "Multiaminoacyl" = "./00-Data/03-Results/COPF/COPF_Complex_Multiaminoacyl-tRNA_results.csv",
#   "26S-proteasome" = "./00-Data/03-Results/COPF/COPF_Complex_26S-proteasome_results.csv",
#   "60S-ribosome" = "./00-Data/03-Results/COPF/COPF_Complex_60S-cytosolic-large-ribosomal-subunit_results.csv",
#   "BSSADCRC" = "./00-Data/03-Results/COPF/COPF_Complex_BSSADCRC_results.csv",
#   "Dynactin" = "./00-Data/03-Results/COPF/COPF_Complex_Dynactin-complex_results.csv",
#   "Dynein-1" = "./00-Data/03-Results/COPF/COPF_Complex_Dynein-1-complex-variant-1_results.csv",
#   "eIF3" = "./00-Data/03-Results/COPF/COPF_Complex_Eukaryotic-translation-initiation-factor-3-complex_results.csv",
#   "IFT-B" = "./00-Data/03-Results/COPF/COPF_Complex_Intraflagellar-transport-complex-B_results.csv",
#   "Laminin-213" = "./00-Data/03-Results/COPF/COPF_Complex_Laminin-213_results.csv",
#   "Spliceosomal-B" = "./00-Data/03-Results/COPF/COPF_Complex_Major-Spliceosomal-B_results.csv",
#   "Nuclear-pore" = "./00-Data/03-Results/COPF/COPF_Complex_Nuclear-pore-complex_results.csv"
# )

# COPF_complex_PCA <- list()

# for (name in names(COPF_complex_results)) {
#   if (!file.exists(COPF_complex_results[[name]])) {
#     warning(paste("File not found:", COPF_complex_results[[name]]))
#     next
#   }

#   copf_res <- fread(COPF_complex_results[[name]])

#   orig <- fread(complex_datasets[[name]])
#   intensity_cols <- grep("^B\\d+_D\\d+$", names(orig), value = TRUE)

#   merged <- merge(
#     copf_res[, .(peptide_id, proteoform)],
#     orig[, c("Peptide", intensity_cols), with = FALSE],
#     by.x = "peptide_id",
#     by.y = "Peptide"
#   )

#   long <- data.table::melt(
#     merged,
#     id.vars = c("peptide_id", "proteoform"),
#     measure.vars = intensity_cols,
#     variable.name = "Sample",
#     value.name = "Intensity"
#   )
#   long[, Condition := sub("^B\\d+_", "", Sample)]
#   long_avg <- long[,
#     .(Intensity = mean(Intensity, na.rm = TRUE)),
#     by = .(peptide_id, proteoform, Condition)
#   ]
#   wide_avg <- dcast(
#     long_avg,
#     peptide_id + proteoform ~ Condition,
#     value.var = "Intensity"
#   )

#   cond_cols <- setdiff(names(wide_avg), c("peptide_id", "proteoform"))
#   mat <- as.data.frame(wide_avg[, cond_cols, with = FALSE])
#   rownames(mat) <- paste0(wide_avg$peptide_id, "_", wide_avg$proteoform)
#   for (col in names(mat)) {
#     mat[[col]][is.nan(mat[[col]])] <- NA
#   }

#   pca <- pcaMethods::pca(
#     mat,
#     method = "ppca",
#     nPcs = 2,
#     scale = "none",
#     center = TRUE
#   )
#   scores <- as.data.frame(pca@scores)
#   scores$ID <- sub(".*_", "", rownames(scores))
#   var_expl <- pca@R2

#   pal <- scico::scico(length(unique(scores$ID)), palette = "batlow")
#   id_levels <- sort(unique(scores$ID))
#   id_colors <- setNames(pal, id_levels)

#   p <- ggplot(scores, aes(x = PC1, y = PC2, color = ID)) +
#     geom_point(size = 2.5, alpha = 0.95) +
#     scale_color_manual(values = id_colors) +
#     labs(
#       x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
#       y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
#       color = "COPF proteoform"
#     ) +
#     coord_fixed() +
#     theme_bw() +
#     ggtitle(paste(name, "- COPF Proteoforms")) +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       aspect.ratio = 1,
#       legend.position = "bottom",
#       legend.title = element_blank()
#     )
#   COPF_complex_PCA[[name]] <- p
# }

# combined_COPF_complex_plot <- plot_grid(plotlist = COPF_complex_PCA, ncol = 5)
# combined_COPF_complex_plot <- plot_grid(
#   ggdraw() +
#     draw_label(
#       "PCA of complex datasets colored by COPF proteoforms",
#       fontface = "bold",
#       size = 20
#     ),
#   combined_COPF_complex_plot,
#   ncol = 1,
#   rel_heights = c(0.05, 1)
# )

# ggsave(
#   filename = file.path(output_dir, "07-ComplexDataCOPFPCA.pdf"),
#   plot = combined_COPF_complex_plot,
#   width = 20,
#   height = 12,
#   dpi = 300
# )

################################
# RPC mapping for complex data #
################################

RCP_datasets <- list(
  "B-WICH_2-2" = "./00-Data/03-Results/RPC/RPC2-2_Complex_B-Wich_results.csv",
  "CRD_2-2" = "./00-Data/03-Results/RPC/RPC2-2_Complex_CRD-mediated_results.csv",
  "Multiaminoacyl_2-2" = "./00-Data/03-Results/RPC/RPC2-2_Complex_Multiaminoacyl-tRNA_results.csv"
)

RPC_complex_PCA <- list()

for (name in names(RCP_datasets)) {
  if (!file.exists(RCP_datasets[[name]])) {
    warning(paste("File not found:", RCP_datasets[[name]]))
    next
  }
  RPC_res <- fread(RCP_datasets[[name]])
  intensity_cols <- grep("^B\\d+_D\\d+$", names(RPC_res), value = TRUE)
  RPC_res <- RPC_res[, c("Peptide", "dCF", intensity_cols), with = FALSE]
  long <- data.table::melt(
    RPC_res,
    id.vars = c("Peptide", "dCF"),
    measure.vars = intensity_cols,
    variable.name = "Sample",
    value.name = "Intensity"
  )
  long[, Condition := sub("^B\\d+_", "", Sample)]
  long_avg <- long[,
    .(Intensity = mean(Intensity, na.rm = TRUE)),
    by = .(Peptide, dCF, Condition)
  ]
  wide_avg <- dcast(
    long_avg,
    Peptide + dCF ~ Condition,
    value.var = "Intensity"
  )
  for (col in names(wide_avg)) {
    wide_avg[[col]][is.nan(wide_avg[[col]])] <- NA
  }
  wide_avg <- as.data.frame(wide_avg)
  wide_avg$ID <- paste0(wide_avg$Peptide, "_", wide_avg$dCF)
  rownames(wide_avg) <- wide_avg$ID
  wide_avg$ID <- NULL
  wide_avg$Peptide <- NULL
  wide_avg$dCF <- NULL
  pca <- pcaMethods::pca(
    wide_avg,
    method = "ppca",
    nPcs = 2,
    scale = "none",
    center = TRUE
  )
  scores <- as.data.frame(pca@scores)
  scores$ID <- sub(".*_", "", rownames(scores))
  var_expl <- pca@R2
  pal <- scico::scico(length(unique(scores$ID)), palette = "batlow")
  id_colors <- setNames(pal, sort(unique(scores$ID)))
  p <- ggplot(scores, aes(x = PC1, y = PC2, color = ID)) +
    geom_point(size = 2.5, alpha = 0.95) +
    scale_color_manual(values = id_colors) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
      color = "RPC dCF"
    ) +
    coord_fixed() +
    theme_bw() +
    ggtitle(paste(name, "- RPC Clusters")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.title = element_blank()
    )
  RPC_complex_PCA[[name]] <- p
}

combined_RPC_complex_plot <- plot_grid(plotlist = RPC_complex_PCA, ncol = 3)
combined_RPC_complex_plot <- plot_grid(
  ggdraw() +
    draw_label(
      "PCA of complex datasets colored by RPC2-2 dCF",
      fontface = "bold",
      size = 20
    ),
  combined_RPC_complex_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)
ggsave(
  filename = file.path(output_dir, "07-ComplexDataRPC2-2PCA.pdf"),
  plot = combined_RPC_complex_plot,
  width = 12,
  height = 6,
  dpi = 300
)

################################################################################
## Tile plots of number of proteoforms vs ProteoMaker grouping #################
################################################################################

# Data paths
ProteoForge_results <- list(
  "LowNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNA_result.feather",
  "MedNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNA_result.feather",
  "HighNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNA_result.feather"
)
DCF_tile_results <- list(
  "LowNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerLowNA_results.csv",
  "MedNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerMedNA_results.csv",
  "HighNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerHighNA_results.csv"
)
DCFBaseline_tile_results <- list(
  "LowNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNA_results.csv",
  "MedNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNA_results.csv",
  "HighNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNA_results.csv"
)
RPC_tile_results <- list(
  "LowNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerLowNA_results.csv",
  "MedNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerMedNA_results.csv",
  "HighNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerHighNA_results.csv"
)
# COPF_tile_results <- list(
#   "LowNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerLowNA_results.csv",
#   "MedNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerMedNA_results.csv",
#   "HighNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerHighNA_results.csv"
# )

# Initialize data table to hold counts
counts_dt <- data.table(
  Dataset = character(),
  Method = character(),
  ProteoMakerProteoformCount = integer(),
  FoundProteoformCount = integer()
)

count_non_singleton_proteoforms <- function(proteoform_id_values) {
  id_values <- as.character(proteoform_id_values)
  id_values <- id_values[!is.na(id_values) & id_values != ""]

  if (!length(id_values)) {
    return(0L)
  }

  split_values <- unlist(strsplit(id_values, split = "\\|"), use.names = FALSE)
  split_values <- split_values[split_values != ""]

  if (!length(split_values)) {
    return(0L)
  }

  as.integer(length(unique(split_values)))
}

# Get counts of identified proteoforms per protein in ProteoForge results
for (name in names(ProteoForge_results)) {
  if (!file.exists(ProteoForge_results[[name]])) {
    warning(paste("File not found:", ProteoForge_results[[name]]))
    next
  }

  # Get ProteoForge grouping results
  PF_res <- as.data.table(read_feather(ProteoForge_results[[name]]))

  # Subset to Sample = "C1_R1"
  PF_res <- PF_res[Sample == "C1_R1"]

  # Count non-singleton clusters (>1 peptide) per protein, matching DCF's dCF-1 exclusion
  cluster_sizes <- PF_res[, .(n_peptides = .N), by = .(protein_id, ClusterID)]
  PF_counts <- cluster_sizes[
    n_peptides > 1,
    .(FoundProteoformCount = .N),
    by = protein_id
  ]

  # Remove proteins with ambigiguous mapping (Protein contains |)
  PF_counts <- PF_counts[!grepl("\\|", protein_id)]

  # Get original ProteoMaker data
  alt_name <- sub("NA$", "NoNA", name)
  protein_data <- fread(ProteoMaker_datasets[[alt_name]])

  # Count number of unique Proteoform_IDs per protein in original data
  for (protein in PF_counts$protein_id) {
    if (!(protein %in% protein_data$Accession)) {
      warning(paste(
        "Protein",
        protein,
        "not found in ProteoMaker data for dataset",
        alt_name
      ))
      next
    }

    # Subset to this protein
    prot_data <- protein_data[Accession == protein]

    # Count Proteoform_IDs that occur more than once
    proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)

    # Store count in PF_counts
    PF_counts[
      protein_id == protein,
      ProteoMakerProteoformCount := proteoform_count
    ]
  }

  # Store counts in data table
  PF_counts[, Dataset := name]
  PF_counts[, Method := "ProteoForge"]
  setnames(PF_counts, "protein_id", "Protein")
  PF_counts <- PF_counts[, .(
    Dataset,
    Method,
    Protein,
    FoundProteoformCount,
    ProteoMakerProteoformCount
  )]
  counts_dt <- rbind(counts_dt, PF_counts, fill = TRUE)
}

# Get counts of identified proteoforms per protein in DCF results
for (name in names(DCF_tile_results)) {
  if (!file.exists(DCF_tile_results[[name]])) {
    warning(paste("File not found:", DCF_tile_results[[name]]))
    next
  }
  dcf_res <- fread(DCF_tile_results[[name]])

  # Count dCF groups with >1 peptide per protein
  dcf_sizes <- dcf_res[!is.na(dCF) & dCF != "dCF-1", .N, by = .(Accession, dCF)]
  dcf_counts <- dcf_sizes[N > 1, .(FoundProteoformCount = .N), by = Accession]
  dcf_counts <- dcf_counts[!grepl("\\|", Accession)]

  alt_name <- sub("NA$", "NoNA", name)
  protein_data <- fread(ProteoMaker_datasets[[alt_name]])

  for (protein in dcf_counts$Accession) {
    if (!(protein %in% protein_data$Accession)) {
      warning(paste(
        "Protein",
        protein,
        "not found in ProteoMaker data for dataset",
        alt_name
      ))
      next
    }
    prot_data <- protein_data[Accession == protein]
    proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)
    dcf_counts[
      Accession == protein,
      ProteoMakerProteoformCount := proteoform_count
    ]
  }

  dcf_counts[, Dataset := name]
  dcf_counts[, Method := "DCF"]
  setnames(dcf_counts, "Accession", "Protein")
  dcf_counts <- dcf_counts[, .(
    Dataset,
    Method,
    Protein,
    FoundProteoformCount,
    ProteoMakerProteoformCount
  )]
  counts_dt <- rbind(counts_dt, dcf_counts, fill = TRUE)
}

# Get counts of identified proteoforms per protein in DCFBaseline results
for (name in names(DCFBaseline_tile_results)) {
  if (!file.exists(DCFBaseline_tile_results[[name]])) {
    warning(paste("File not found:", DCFBaseline_tile_results[[name]]))
    next
  }
  dcf_res <- fread(DCFBaseline_tile_results[[name]])

  dcf_sizes <- dcf_res[!is.na(dCF) & dCF != "dCF-1", .N, by = .(Accession, dCF)]
  dcf_counts <- dcf_sizes[N > 1, .(FoundProteoformCount = .N), by = Accession]
  dcf_counts <- dcf_counts[!grepl("\\|", Accession)]

  alt_name <- sub("NA$", "NoNA", name)
  protein_data <- fread(ProteoMaker_datasets[[alt_name]])

  for (protein in dcf_counts$Accession) {
    if (!(protein %in% protein_data$Accession)) {
      warning(paste(
        "Protein",
        protein,
        "not found in ProteoMaker data for dataset",
        alt_name
      ))
      next
    }
    prot_data <- protein_data[Accession == protein]
    proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)
    dcf_counts[
      Accession == protein,
      ProteoMakerProteoformCount := proteoform_count
    ]
  }

  dcf_counts[, Dataset := name]
  dcf_counts[, Method := "DCFBaseline"]
  setnames(dcf_counts, "Accession", "Protein")
  dcf_counts <- dcf_counts[, .(
    Dataset,
    Method,
    Protein,
    FoundProteoformCount,
    ProteoMakerProteoformCount
  )]
  counts_dt <- rbind(counts_dt, dcf_counts, fill = TRUE)
}

# Get counts of identififed proteoforms per protein in RPC results
for (name in names(RPC_tile_results)) {
  if (!file.exists(RPC_tile_results[[name]])) {
    warning(paste("File not found:", RPC_tile_results[[name]]))
    next
  }
  rpc_res <- fread(RPC_tile_results[[name]])

  rpc_sizes <- rpc_res[!is.na(dCF) & dCF != "dCF-1", .N, by = .(Accession, dCF)]
  rpc_counts <- rpc_sizes[N > 1, .(FoundProteoformCount = .N), by = Accession]
  rpc_counts <- rpc_counts[!grepl("\\|", Accession)]

  alt_name <- sub("NA$", "NoNA", name)
  protein_data <- fread(ProteoMaker_datasets[[alt_name]])

  for (protein in rpc_counts$Accession) {
    if (!(protein %in% protein_data$Accession)) {
      warning(paste(
        "Protein",
        protein,
        "not found in ProteoMaker data for dataset",
        alt_name
      ))
      next
    }
    prot_data <- protein_data[Accession == protein]
    proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)
    rpc_counts[
      Accession == protein,
      ProteoMakerProteoformCount := proteoform_count
    ]
  }
  rpc_counts[, Dataset := name]
  rpc_counts[, Method := "RPC"]
  setnames(rpc_counts, "Accession", "Protein")
  rpc_counts <- rpc_counts[, .(
    Dataset,
    Method,
    Protein,
    FoundProteoformCount,
    ProteoMakerProteoformCount
  )]
  counts_dt <- rbind(counts_dt, rpc_counts, fill = TRUE)
}

# Get counts of identified proteoforms per protein in COPF results
# for (name in names(COPF_tile_results)) {
#   if (!file.exists(COPF_tile_results[[name]])) {
#     warning(paste("File not found:", COPF_tile_results[[name]]))
#     next
#   }
#   copf_res <- fread(COPF_tile_results[[name]])

#   # uniqueN(proteoform) per protein: 1 = canonical only, 2 = has discordant proteoform
#   copf_counts <- copf_res[
#     !is.na(proteoform),
#     .(FoundProteoformCount = uniqueN(proteoform)),
#     by = protein_id
#   ]
#   copf_counts <- copf_counts[!grepl("\\|", protein_id)]

#   alt_name <- sub("NA$", "NoNA", name)
#   protein_data <- fread(ProteoMaker_datasets[[alt_name]])

#   for (protein in copf_counts$protein_id) {
#     if (!(protein %in% protein_data$Accession)) {
#       warning(paste(
#         "Protein",
#         protein,
#         "not found in ProteoMaker data for dataset",
#         alt_name
#       ))
#       next
#     }
#     prot_data <- protein_data[Accession == protein]
#     proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)
#     copf_counts[
#       protein_id == protein,
#       ProteoMakerProteoformCount := proteoform_count
#     ]
#   }

#   copf_counts[, Dataset := name]
#   copf_counts[, Method := "COPF"]
#   setnames(copf_counts, "protein_id", "Protein")
#   copf_counts <- copf_counts[, .(
#     Dataset,
#     Method,
#     Protein,
#     FoundProteoformCount,
#     ProteoMakerProteoformCount
#   )]
#   counts_dt <- rbind(counts_dt, copf_counts, fill = TRUE)
# }

# Plot scatterplot of FoundProteoformCount vs ProteoMakerProteoformCount
plot_scatter_counts <- function(data, dataset_name, method_name, color) {
  # Subset counts_dt to this dataset and method
  subset_data <- data[
    Dataset == dataset_name &
      Method == method_name &
      is.finite(ProteoMakerProteoformCount) &
      is.finite(FoundProteoformCount)
  ]

  # Use an individual square axis per subplot
  max_val <- max(
    c(subset_data$ProteoMakerProteoformCount, subset_data$FoundProteoformCount),
    na.rm = TRUE
  )

  # Count observations in each integer square
  tile_counts <- subset_data[,
    .(Count = .N),
    by = .(
      ProteoMakerProteoformCount = as.integer(ProteoMakerProteoformCount),
      FoundProteoformCount = as.integer(FoundProteoformCount)
    )
  ]

  # Complete grid so empty squares are shown as zero
  tile_grid <- CJ(
    ProteoMakerProteoformCount = 0:max_val,
    FoundProteoformCount = 0:max_val
  )

  tile_data <- merge(
    tile_grid,
    tile_counts,
    by = c("ProteoMakerProteoformCount", "FoundProteoformCount"),
    all.x = TRUE
  )
  tile_data[is.na(Count), Count := 0L]

  p <- ggplot(
    tile_data,
    aes(x = ProteoMakerProteoformCount, y = FoundProteoformCount, fill = Count)
  ) +
    geom_tile(color = "gray95", linewidth = 0.2) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray",
      linewidth = 1
    ) +
    scale_fill_viridis_c(option = "C", trans = "sqrt") +
    scale_x_continuous(
      limits = c(-0.5, max_val + 0.5),
      breaks = seq(0, max_val, by = 5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(-0.5, max_val + 0.5),
      breaks = seq(0, max_val, by = 5),
      expand = c(0, 0)
    ) +
    labs(
      title = paste(dataset_name),
      x = "ProteoMaker Proteoform Count",
      y = "Found Proteoform Count",
      fill = "Protein count"
    ) +
    coord_equal() +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1
    )
  return(p)
}

# Plot histogram of ProteoMakerProteoformCount - FoundProteoformCount
plot_hist_diff <- function(data, method_name, color, y_max = NULL) {
  subset_data <- data[Method == method_name]
  subset_data[, diff := FoundProteoformCount - ProteoMakerProteoformCount]

  p <- ggplot(
    subset_data,
    aes(x = diff)
  ) +
    geom_bar(
      fill = color,
      color = "white",
      alpha = 0.9
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      title = paste(method_name, "diff"),
      x = "Found Count - ProteoMaker Count",
      y = "Count"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1
    )

  if (!is.null(y_max)) {
    p <- p + coord_cartesian(ylim = c(0, y_max))
  }
  return(p)
}

# Create tile plots for each dataset and method
scatter_plots <- list()
hist_plots <- list()
method_levels <- c("ProteoForge", "DCF", "DCFBaseline", "RPC") # COPF
colors <- scico::scico(4, palette = "batlow", end = 0.7)
for (method in method_levels) {
  for (dataset in c("LowNA", "MedNA", "HighNA")) {
    plot_name <- paste(dataset, method, sep = "_")
    scatter_plots[[plot_name]] <- plot_scatter_counts(
      data = counts_dt,
      dataset_name = dataset,
      method_name = method,
      color = colors[which(method_levels == method)]
    )
  }
  hist_plots[[method]] <- plot_hist_diff(
    data = counts_dt,
    method_name = method,
    color = colors[which(method_levels == method)]
  )
}

make_method_row <- function(method, scatter_plots, hist_plots) {
  row <- plot_grid(
    plotlist = c(
      scatter_plots[grepl(paste0("_", method, "$"), names(scatter_plots))],
      list(hist_plots[[method]])
    ),
    ncol = 4
  )
  plot_grid(
    ggdraw() + draw_label(method, fontface = "bold", size = 20),
    row,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
}

# Combine all rows into one plot
final_combined_plot <- plot_grid(
  make_method_row("ProteoForge", scatter_plots, hist_plots),
  make_method_row("DCF", scatter_plots, hist_plots),
  make_method_row("DCFBaseline", scatter_plots, hist_plots),
  make_method_row("RPC", scatter_plots, hist_plots),
  ncol = 1
)

ggsave(
  filename = file.path(
    output_dir,
    "07-ProteoMakerMappingCompare.pdf"
  ),
  plot = final_combined_plot,
  width = 16,
  height = 12,
  dpi = 300
)

########################################
# Diff bar plots: Found - ProteoMaker  #
########################################

diff_methods <- c("ProteoForge", "DCF", "DCFBaseline")
diff_colors <- scico::scico(3, palette = "batlow", end = 0.7)

diff_y_max <- counts_dt[
  Method %in%
    diff_methods &
    is.finite(FoundProteoformCount) &
    is.finite(ProteoMakerProteoformCount),
  .(diff = FoundProteoformCount - ProteoMakerProteoformCount),
  by = Method
][, max(table(diff)), by = Method][, max(V1)]

diff_bar_plots <- mapply(
  FUN = plot_hist_diff,
  method_name = diff_methods,
  color = diff_colors,
  MoreArgs = list(data = counts_dt, y_max = diff_y_max),
  SIMPLIFY = FALSE
)

diff_combined_plot <- plot_grid(
  plot_grid(plotlist = diff_bar_plots, ncol = 3),
  ncol = 1
)

ggsave(
  filename = file.path(output_dir, "08-ProteoMakerDiffBarPlots.pdf"),
  plot = diff_combined_plot,
  width = 12,
  height = 4,
  dpi = 300
)

########################################################
# Diff bar plots: MedNA + HighNA, ProteoForge/DCF/DCFBaseline
########################################################

counts_dt_MedHighNA <- counts_dt[
  Dataset %in% c("MedNA", "HighNA") & Method %in% diff_methods
]

diff_y_max_MedHighNA <- counts_dt_MedHighNA[
  is.finite(FoundProteoformCount) & is.finite(ProteoMakerProteoformCount),
  .(diff = FoundProteoformCount - ProteoMakerProteoformCount),
  by = Method
][, max(table(diff)), by = Method][, max(V1)]

diff_bar_plots_MedHighNA <- mapply(
  FUN = plot_hist_diff,
  method_name = diff_methods,
  color = diff_colors,
  MoreArgs = list(data = counts_dt_MedHighNA, y_max = diff_y_max_MedHighNA),
  SIMPLIFY = FALSE
)

ggsave(
  filename = file.path(output_dir, "08b-ProteoMakerDiffBarPlots-MedHighNA.pdf"),
  plot = plot_grid(plotlist = diff_bar_plots_MedHighNA, ncol = 3),
  width = 12,
  height = 4,
  dpi = 300
)

final_combined_plot

################################################################################
## Diff bar plots for NoNA datasets ############################################
################################################################################

ProteoForge_NoNA_tile_results <- list(
  "LowNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNoNA_result.feather",
  "MedNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNoNA_result.feather",
  "HighNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNoNA_result.feather"
)
DCF_NoNA_tile_results <- list(
  "LowNoNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerLowNoNA_results.csv",
  "MedNoNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerMedNoNA_results.csv",
  "HighNoNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerHighNoNA_results.csv"
)
DCFBaseline_NoNA_tile_results <- list(
  "LowNoNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNoNA_results.csv",
  "MedNoNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNoNA_results.csv",
  "HighNoNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNoNA_results.csv"
)
RPC_NoNA_tile_results <- list(
  "LowNoNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerLowNoNA_results.csv",
  "MedNoNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerMedNoNA_results.csv",
  "HighNoNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerHighNoNA_results.csv"
)

counts_dt_NoNA <- data.table(
  Dataset = character(),
  Method = character(),
  ProteoMakerProteoformCount = integer(),
  FoundProteoformCount = integer()
)

for (name in names(ProteoForge_NoNA_tile_results)) {
  if (!file.exists(ProteoForge_NoNA_tile_results[[name]])) {
    warning(paste("File not found:", ProteoForge_NoNA_tile_results[[name]]))
    next
  }
  PF_res <- as.data.table(read_feather(ProteoForge_NoNA_tile_results[[name]]))
  PF_res <- PF_res[Sample == "C1_R1"]
  cluster_sizes <- PF_res[, .(n_peptides = .N), by = .(protein_id, ClusterID)]
  PF_counts <- cluster_sizes[
    n_peptides > 1,
    .(FoundProteoformCount = .N),
    by = protein_id
  ]
  PF_counts <- PF_counts[!grepl("\\|", protein_id)]

  protein_data <- fread(ProteoMaker_datasets[[name]])

  for (protein in PF_counts$protein_id) {
    if (!(protein %in% protein_data$Accession)) {
      warning(paste(
        "Protein",
        protein,
        "not found in ProteoMaker data for dataset",
        name
      ))
      next
    }
    prot_data <- protein_data[Accession == protein]
    proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)
    PF_counts[
      protein_id == protein,
      ProteoMakerProteoformCount := proteoform_count
    ]
  }

  PF_counts[, Dataset := name]
  PF_counts[, Method := "ProteoForge"]
  setnames(PF_counts, "protein_id", "Protein")
  PF_counts <- PF_counts[, .(
    Dataset,
    Method,
    Protein,
    FoundProteoformCount,
    ProteoMakerProteoformCount
  )]
  counts_dt_NoNA <- rbind(counts_dt_NoNA, PF_counts, fill = TRUE)
}

for (name in names(DCF_NoNA_tile_results)) {
  if (!file.exists(DCF_NoNA_tile_results[[name]])) {
    warning(paste("File not found:", DCF_NoNA_tile_results[[name]]))
    next
  }
  dcf_res <- fread(DCF_NoNA_tile_results[[name]])
  dcf_sizes <- dcf_res[!is.na(dCF) & dCF != "dCF-1", .N, by = .(Accession, dCF)]
  dcf_counts <- dcf_sizes[N > 1, .(FoundProteoformCount = .N), by = Accession]
  dcf_counts <- dcf_counts[!grepl("\\|", Accession)]

  protein_data <- fread(ProteoMaker_datasets[[name]])

  for (protein in dcf_counts$Accession) {
    if (!(protein %in% protein_data$Accession)) {
      warning(paste(
        "Protein",
        protein,
        "not found in ProteoMaker data for dataset",
        name
      ))
      next
    }
    prot_data <- protein_data[Accession == protein]
    proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)
    dcf_counts[
      Accession == protein,
      ProteoMakerProteoformCount := proteoform_count
    ]
  }

  dcf_counts[, Dataset := name]
  dcf_counts[, Method := "DCF"]
  setnames(dcf_counts, "Accession", "Protein")
  dcf_counts <- dcf_counts[, .(
    Dataset,
    Method,
    Protein,
    FoundProteoformCount,
    ProteoMakerProteoformCount
  )]
  counts_dt_NoNA <- rbind(counts_dt_NoNA, dcf_counts, fill = TRUE)
}

for (name in names(DCFBaseline_NoNA_tile_results)) {
  if (!file.exists(DCFBaseline_NoNA_tile_results[[name]])) {
    warning(paste("File not found:", DCFBaseline_NoNA_tile_results[[name]]))
    next
  }
  dcf_res <- fread(DCFBaseline_NoNA_tile_results[[name]])
  dcf_sizes <- dcf_res[!is.na(dCF) & dCF != "dCF-1", .N, by = .(Accession, dCF)]
  dcf_counts <- dcf_sizes[N > 1, .(FoundProteoformCount = .N), by = Accession]
  dcf_counts <- dcf_counts[!grepl("\\|", Accession)]

  protein_data <- fread(ProteoMaker_datasets[[name]])

  for (protein in dcf_counts$Accession) {
    if (!(protein %in% protein_data$Accession)) {
      warning(paste(
        "Protein",
        protein,
        "not found in ProteoMaker data for dataset",
        name
      ))
      next
    }
    prot_data <- protein_data[Accession == protein]
    proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)
    dcf_counts[
      Accession == protein,
      ProteoMakerProteoformCount := proteoform_count
    ]
  }

  dcf_counts[, Dataset := name]
  dcf_counts[, Method := "DCFBaseline"]
  setnames(dcf_counts, "Accession", "Protein")
  dcf_counts <- dcf_counts[, .(
    Dataset,
    Method,
    Protein,
    FoundProteoformCount,
    ProteoMakerProteoformCount
  )]
  counts_dt_NoNA <- rbind(counts_dt_NoNA, dcf_counts, fill = TRUE)
}

for (name in names(RPC_NoNA_tile_results)) {
  if (!file.exists(RPC_NoNA_tile_results[[name]])) {
    warning(paste("File not found:", RPC_NoNA_tile_results[[name]]))
    next
  }
  rpc_res <- fread(RPC_NoNA_tile_results[[name]])
  rpc_sizes <- rpc_res[!is.na(dCF) & dCF != "dCF-1", .N, by = .(Accession, dCF)]
  rpc_counts <- rpc_sizes[N > 1, .(FoundProteoformCount = .N), by = Accession]
  rpc_counts <- rpc_counts[!grepl("\\|", Accession)]

  protein_data <- fread(ProteoMaker_datasets[[name]])

  for (protein in rpc_counts$Accession) {
    if (!(protein %in% protein_data$Accession)) {
      warning(paste(
        "Protein",
        protein,
        "not found in ProteoMaker data for dataset",
        name
      ))
      next
    }
    prot_data <- protein_data[Accession == protein]
    proteoform_count <- count_non_singleton_proteoforms(prot_data$Proteoform_ID)
    rpc_counts[
      Accession == protein,
      ProteoMakerProteoformCount := proteoform_count
    ]
  }

  rpc_counts[, Dataset := name]
  rpc_counts[, Method := "RPC"]
  setnames(rpc_counts, "Accession", "Protein")
  rpc_counts <- rpc_counts[, .(
    Dataset,
    Method,
    Protein,
    FoundProteoformCount,
    ProteoMakerProteoformCount
  )]
  counts_dt_NoNA <- rbind(counts_dt_NoNA, rpc_counts, fill = TRUE)
}

diff_y_max_NoNA <- counts_dt_NoNA[
  Method %in%
    diff_methods &
    is.finite(FoundProteoformCount) &
    is.finite(ProteoMakerProteoformCount),
  .(diff = FoundProteoformCount - ProteoMakerProteoformCount),
  by = Method
][, max(table(diff)), by = Method][, max(V1)]

diff_bar_plots_NoNA <- mapply(
  FUN = plot_hist_diff,
  method_name = diff_methods,
  color = diff_colors,
  MoreArgs = list(data = counts_dt_NoNA, y_max = diff_y_max_NoNA),
  SIMPLIFY = FALSE
)

diff_combined_plot_NoNA <- plot_grid(
  plot_grid(plotlist = diff_bar_plots_NoNA, ncol = 3),
  ncol = 1
)

ggsave(
  filename = file.path(output_dir, "09-ProteoMakerDiffBarPlotsNoNA.pdf"),
  plot = diff_combined_plot_NoNA,
  width = 12,
  height = 4,
  dpi = 300
)
