# Libraries
library(data.table)
library(scico)
library(ggplot2)
library(cowplot)

# Set working directory
setwd(here::here())

# Set output directory for figures
output_dir <- "00-ThesisFigures/Figures"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

datatype_colors <- c(
  "NM" = "#00441b",
  "NMpeptide" = "#02d054ff",
  "RmCys" = "#006d2c",
  "Deglyco" = "#99d38d",
  "Phospho" = "#be4e16",
  "FreeCys" = "#ff8000",
  "LysAc" = "#fec389"
)

# Function to make plot
cluster_from_deviation_with_plot <- function(
  DT,
  cond_regex = ".*_(D[0-9]+)$",
  agg_fun = function(x) mean(x, na.rm = TRUE),
  deep_split = 2,
  minClusterSize = 2,
  singleton_label = "dCF-1",
  scale_rows = FALSE,
  plot = TRUE,
  plot_distance = FALSE,
  plot_title = "Clustered Deviation Heatmap"
) {
  # ---- dependencies ----
  if (!requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    stop("Please install dynamicTreeCut")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table")
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Please install ComplexHeatmap")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Please install circlize")
  }

  DT <- data.table::as.data.table(DT)

  # ------------------------------------------------------------
  # 1. Validate input
  # ------------------------------------------------------------
  if (!"Identifier" %in% names(DT)) {
    stop("Input must contain 'Identifier'")
  }

  intensity_cols <- setdiff(names(DT), "Identifier")

  # ------------------------------------------------------------
  # 2. Parse conditions
  # ------------------------------------------------------------
  cond_labels <- sub(cond_regex, "\\1", intensity_cols, perl = TRUE)

  if (any(cond_labels == intensity_cols)) {
    stop("Failed to parse conditions with cond_regex")
  }

  cond_levels <- unique(cond_labels)

  # ------------------------------------------------------------
  # 3. Aggregate replicates → Bd
  # ------------------------------------------------------------
  Bd <- sapply(cond_levels, function(cc) {
    cols <- intensity_cols[cond_labels == cc]

    apply(
      as.matrix(DT[, ..cols]),
      1,
      function(x) {
        vals <- x[is.finite(x)]
        if (length(vals) == 0) {
          return(NA_real_)
        }
        agg_fun(vals)
      }
    )
  })

  Bd <- as.matrix(Bd)
  rownames(Bd) <- DT$Identifier

  # Optional row scaling
  if (scale_rows) {
    Bd <- t(scale(t(Bd)))
  }

  # ------------------------------------------------------------
  # 4. DISTANCE + CLUSTERING (UNCHANGED)
  # ------------------------------------------------------------
  Kd <- ncol(Bd)
  n <- nrow(Bd)

  distM <- matrix(0, n, n, dimnames = list(rownames(Bd), rownames(Bd)))

  for (ii in seq_len(n - 1)) {
    xi <- Bd[ii, ]
    for (jj in (ii + 1):n) {
      xj <- Bd[jj, ]

      ok <- !is.na(xi) & !is.na(xj)
      ks <- sum(ok)

      if (ks >= 2) {
        dd <- sqrt(sum((xi[ok] - xj[ok])^2))
        dd <- dd * sqrt(Kd / ks)
      } else {
        dd <- Inf
      }

      distM[ii, jj] <- dd
      distM[jj, ii] <- dd
    }
  }

  if (any(!is.finite(distM))) {
    finite_vals <- distM[is.finite(distM)]
    big <- if (length(finite_vals) > 0) max(finite_vals) else 1
    distM[!is.finite(distM)] <- big * 10
  }

  hc <- stats::hclust(stats::as.dist(distM), method = "ward.D2")

  if (n >= 3) {
    cl <- dynamicTreeCut::cutreeDynamic(
      dendro = hc,
      distM = distM,
      method = "hybrid",
      deepSplit = deep_split,
      minClusterSize = minClusterSize,
      verbose = 0
    )
  } else {
    cl <- seq_len(n)
  }

  cl <- as.integer(cl)
  names(cl) <- rownames(Bd)

  # ------------------------------------------------------------
  # 5. LABEL ASSIGNMENT (UNCHANGED)
  # ------------------------------------------------------------
  cl_dt <- data.table::data.table(
    Identifier = names(cl),
    cluster_id = cl
  )

  sizes <- cl_dt[, .N, by = cluster_id]

  good_ids <- sort(unique(
    sizes[cluster_id > 0 & N >= 2, cluster_id]
  ))

  cl_dt[, label := singleton_label]

  if (length(good_ids) > 0) {
    map <- data.table::data.table(
      cluster_id = good_ids,
      label = paste0("dCF", seq_along(good_ids))
    )

    cl_dt <- merge(
      cl_dt,
      map,
      by = "cluster_id",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("", ".m")
    )

    cl_dt[!is.na(label.m), label := label.m]
    cl_dt[, label.m := NULL]
  }

  labels <- stats::setNames(cl_dt$label, cl_dt$Identifier)

  # ------------------------------------------------------------
  # 6. ORDER rows by clustering
  # ------------------------------------------------------------
  ord <- hc$order
  Bd_ord <- Bd[ord, , drop = FALSE]
  cl_ord <- cl[ord]
  lab_ord <- labels[rownames(Bd_ord)]

  # ------------------------------------------------------------
  # 7. ComplexHeatmap plot (RIGHT annotation + colored dendrogram)
  # ------------------------------------------------------------
  if (plot) {
    library(ComplexHeatmap)
    library(circlize)
    library(dendextend)

    # ------------------------------------------------------------
    # Cluster color mapping (shared by sidebar and dendrogram)
    # ------------------------------------------------------------
    cluster_levels <- unique(lab_ord)

    cluster_cols <- setNames(
      c(
        "#E41A1C",
        "#377EB8",
        "#4DAF4A",
        "#FF7F00",
        "#984EA3",
        "#A65628",
        "#F781BF",
        "#FFFF33"
      ),
      cluster_levels
    )

    # ------------------------------------------------------------
    # Build dendrogram from hc (leaves now labeled by identifier)
    # ------------------------------------------------------------
    dend <- as.dendrogram(hc)

    # ------------------------------------------------------------
    # Color dendrogram branches to match sidebar colors
    # ------------------------------------------------------------
    leaf_labels_ordered <- labels[labels(dend)]
    cl_for_dend <- match(leaf_labels_ordered, cluster_levels)

    dend <- dendextend::color_branches(
      dend,
      clusters = cl_for_dend,
      col = unname(cluster_cols)
    )

    # Annotation in Bd's original row order; ComplexHeatmap reorders with dend
    row_labels <- labels[rownames(Bd)]

    # ------------------------------------------------------------
    # RIGHT-SIDE cluster annotation (legend hidden)
    # ------------------------------------------------------------
    right_ha <- rowAnnotation(
      Cluster = row_labels,
      col = list(Cluster = cluster_cols),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )

    # ------------------------------------------------------------
    # Build heatmap
    # ------------------------------------------------------------
    if (plot_distance) {
      distM_ord <- distM[ord, ord]

      dist_vals <- distM[is.finite(distM)]
      dist_q_breaks <- quantile(dist_vals, probs = seq(0, 1, by = 0.25))
      col_fun <- circlize::colorRamp2(
        dist_q_breaks,
        rev(scico::scico(5, palette = "vik", begin = 0.5, end = 1.0))
      )

      legend_at <- seq(0, max(dist_vals), length.out = 5)

      ht <- Heatmap(
        distM_ord,
        name = "Distance",
        col = col_fun,

        cluster_rows = dend,
        cluster_columns = dend,

        show_row_names = FALSE,
        show_column_names = FALSE,

        right_annotation = right_ha,
        column_title = plot_title,

        heatmap_legend_param = list(
          title = "Distance",
          at = legend_at,
          labels = round(legend_at, 2),
          color_bar = "continuous"
        )
      )
    } else {
      col_fun <- circlize::colorRamp2(
        c(-2, -1, 0, 1, 2),
        scico::scico(5, palette = "vik")
      )

      ht <- Heatmap(
        Bd,
        name = "Deviation",
        col = col_fun,

        cluster_rows = dend,
        cluster_columns = FALSE,

        show_row_names = FALSE,
        show_column_names = TRUE,

        column_title = plot_title
      )
    }

    pdf(
      file.path(output_dir, paste0(gsub(" ", "_", plot_title), ".pdf")),
      width = 5,
      height = 6
    )
    draw(ht, heatmap_legend_side = "bottom")
    dev.off()
  }

  return(list(
    labels = labels,
    cluster_id = cl,
    ht = if (plot) ht else NULL,
    hc = hc,
    distM = distM,
    Bd = Bd
  ))
}

######################################
## Discover complexoforms base idea ##
######################################

data <- readRDS(
  "ComplexoFinder/10-Results/B-WICH_chromatin_remodelling_complex/complexoform_results.rds"
)

data <- data$dCF$peptide_results

data <- data[which(data$dCF != "dCF0")]

intensity_cols <- grep("H9", names(data), value = TRUE)
data <- data[, c("Identifier", intensity_cols), with = FALSE]

result <- cluster_from_deviation_with_plot(
  data,
  agg_fun = function(x) mean(x, na.rm = TRUE),
  deep_split = 2,
  minClusterSize = 2,
  singleton_label = "dCF-1",
  plot_distance = TRUE,
  plot_title = "Determine Number of Complexoforms"
)

# Draw heatmap
if (!is.null(result$ht)) {
  pdf(
    file.path(output_dir, "Complexoform_Clustering.pdf"),
    width = 5,
    height = 6
  )
  draw(result$ht, heatmap_legend_side = "bottom")
  dev.off()
}
