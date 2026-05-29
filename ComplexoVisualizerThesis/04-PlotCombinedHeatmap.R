# Libraries
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(scico)

setwd(here::here())

source("ComplexoVisualizer/MakeCXC.R")

output_dir <- "ComplexoVisualiserThesis/00-Plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

results_dir  <- "ComplexoVisualiserThesis/00-Data/Results"
summaries_dir <- "ComplexoVisualiserThesis/00-Data/Summaries"

################################################################################
## Plotting function ###########################################################
################################################################################

plot_combined_distance_heatmap <- function(
  distance_csv,
  cluster_summary,
  mode = c("threshold", "highest"),
  threshold = 0.5,
  col_order = c("clustered", "by_assignment", "both_by_assignment", "by_protein"),
  title = "",
  cluster_method = "complete"
) {
  mode      <- match.arg(mode)
  col_order <- match.arg(col_order)

  dist_dt <- data.table::fread(distance_csv)
  cs      <- data.table::fread(cluster_summary)
  nm_types <- c("NM", "NMpeptide")
  cs_ptm  <- cs[!Datatype %in% nm_types]

  mem_cols <- grep("^membership of cluster", names(cs_ptm), value = TRUE)
  k <- length(mem_cols)

  if (k > 0L) {
    mem_mat <- as.matrix(cs_ptm[, .SD, .SDcols = mem_cols])
    cs_ptm[, max_mem  := apply(mem_mat, 1L, max)]
    cs_ptm[, cl_assign := apply(mem_mat, 1L, which.max)]
    cs_ptm[max_mem == 0, cl_assign := NA_integer_]
  } else {
    cs_ptm[, max_mem   := 1.0]
    cs_ptm[, cl_assign := ClusterID]
    k <- length(unique(cs_ptm$ClusterID))
  }

  if (mode == "threshold") cs_ptm <- cs_ptm[max_mem >= threshold]

  if (nrow(cs_ptm) < 2L) {
    warning("Fewer than 2 features after filtering (mode='", mode, "') — skipping.")
    return(invisible(NULL))
  }

  cs_ptm[, feat_id := paste0(Position, "\n", Datatype)]

  # Multiple modification forms of the same peptide share the same feat_id and
  # the same structural position in the distance file.  Keep only the form with
  # the highest max_mem so each (Position, Datatype) appears exactly once.
  cs_ptm <- cs_ptm[cs_ptm[, .I[which.max(max_mem)], by = feat_id]$V1]

  feats <- cs_ptm$feat_id
  n     <- length(feats)

  dist_mat <- matrix(NA_real_, n, n, dimnames = list(feats, feats))
  diag(dist_mat) <- 0.0

  pairs <- data.table::copy(dist_dt)
  if (!"SameCoordFrame" %in% names(pairs)) pairs[, SameCoordFrame := TRUE]
  pairs <- pairs[SameCoordFrame == TRUE]
  pairs[, SameCoordFrame := NULL]
  pairs[, feat_id1 := paste0(Position1, "\n", Datatype1)]
  pairs[, feat_id2 := paste0(Position2, "\n", Datatype2)]
  pairs <- pairs[feat_id1 %in% feats & feat_id2 %in% feats]
  pairs <- pairs[,
    .(Distance_Angstrom = min(Distance_Angstrom)),
    by = .(feat_id1, feat_id2)
  ]

  for (i in seq_len(nrow(pairs))) {
    p1 <- pairs$feat_id1[i]
    p2 <- pairs$feat_id2[i]
    d  <- pairs$Distance_Angstrom[i]
    dist_mat[p1, p2] <- d
    dist_mat[p2, p1] <- d
  }

  mapped <- apply(dist_mat, 1L, function(row) any(!is.na(row) & row > 0))
  if (!all(mapped)) {
    n_drop <- sum(!mapped)
    message(
      "Dropping ", n_drop, " feature(s) with no structure coverage: ",
      paste(feats[!mapped], collapse = ", ")
    )
    keep_idx <- which(mapped)
    feats    <- feats[keep_idx]
    cs_ptm   <- cs_ptm[keep_idx]
    dist_mat <- dist_mat[feats, feats]
    n        <- length(feats)
    if (n < 2L) {
      warning("Fewer than 2 features with structure coverage — skipping.")
      return(invisible(NULL))
    }
  }

  max_d      <- max(dist_mat, na.rm = TRUE)
  dist_clean <- dist_mat
  dist_clean[is.na(dist_clean)] <- max_d * 1.5
  hc     <- hclust(as.dist(dist_clean), method = cluster_method)
  hc_ord <- hc$order

  assign_ord  <- order(cs_ptm$cl_assign, na.last = TRUE)
  protein_ord <- order(cs_ptm[["Gene name"]], na.last = TRUE)

  row_ord <- switch(col_order,
    clustered          = hc_ord,
    by_assignment      = hc_ord,
    both_by_assignment = assign_ord,
    by_protein         = protein_ord
  )
  col_ord_idx <- switch(col_order,
    clustered          = hc_ord,
    by_assignment      = assign_ord,
    both_by_assignment = assign_ord,
    by_protein         = protein_ord
  )

  row_feats   <- feats[row_ord]
  col_feats   <- feats[col_ord_idx]
  ordered_mat <- dist_mat[row_feats, col_feats]
  row_cl      <- as.character(cs_ptm$cl_assign[row_ord])
  col_cl      <- as.character(cs_ptm$cl_assign[col_ord_idx])
  row_genes   <- cs_ptm[["Gene name"]][row_ord]
  col_genes   <- cs_ptm[["Gene name"]][col_ord_idx]
  row_dtypes  <- cs_ptm$Datatype[row_ord]
  col_dtypes  <- cs_ptm$Datatype[col_ord_idx]

  cl_colors        <- cluster_palette(k)
  names(cl_colors) <- as.character(seq_len(k))
  if (any(is.na(c(row_cl, col_cl)))) cl_colors <- c(cl_colors, "NA" = "#CCCCCC")

  unique_genes <- unique(c(row_genes, col_genes))
  gene_colors  <- setNames(
    scico(length(unique_genes), palette = "batlow", end = 0.9),
    unique_genes
  )

  dt_pal      <- datatype_palette()
  present_dts <- unique(c(row_dtypes, col_dtypes))
  dt_colors   <- dt_pal[intersect(names(dt_pal), present_dts)]
  unknown_dts <- setdiff(present_dts, names(dt_colors))
  if (length(unknown_dts) > 0L) dt_colors[unknown_dts] <- "#CCCCCC"

  ann_col_list <- list(Cluster = cl_colors, Gene = gene_colors, Datatype = dt_colors)

  top_ann <- HeatmapAnnotation(
    Cluster  = col_cl,
    Gene     = col_genes,
    Datatype = col_dtypes,
    col      = ann_col_list,
    show_legend          = c(TRUE, TRUE, TRUE),
    annotation_name_side = "left",
    annotation_legend_param = list(
      Cluster  = list(title = "Cluster"),
      Gene     = list(title = "Gene"),
      Datatype = list(title = "Datatype")
    )
  )

  left_ann <- rowAnnotation(
    Cluster  = row_cl,
    Gene     = row_genes,
    Datatype = row_dtypes,
    col      = ann_col_list,
    show_legend          = FALSE,
    annotation_name_side = "bottom"
  )

  col_fun <- colorRamp2(c(0, max_d), c("#590007", "#EBE5E0"))

  col_order_label <- switch(col_order,
    clustered          = "both axes clustered",
    by_assignment      = "rows clustered  |  columns by assignment",
    both_by_assignment = "both axes by cluster assignment",
    by_protein         = "both axes by protein"
  )
  mode_label <- if (mode == "threshold") {
    paste0(
      "membership threshold >= ", threshold,
      "  |  ", col_order_label,
      "  (n=", n, " features)"
    )
  } else {
    paste0(
      "highest-membership assignment",
      "  |  ", col_order_label,
      "  (n=", n, " features)"
    )
  }

  Heatmap(
    ordered_mat,
    name              = "Distance (Å)",
    col               = col_fun,
    na_col            = "white",
    top_annotation    = top_ann,
    left_annotation   = left_ann,
    column_title      = paste0(title, "\n", mode_label),
    show_row_names    = FALSE,
    show_column_names = FALSE,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_dend_reorder  = FALSE,
    column_title_gp   = grid::gpar(fontsize = 11)
  )
}

################################################################################
## Input definitions ###########################################################
################################################################################

# Map short label -> (summary file, full complex name for folder derivation)
# Folder names are derived with the same gsub used in 02-VisualiseStructures.R.
make_complex_entry <- function(label, summary_file, full_name) {
  folder <- gsub("[^A-Za-z0-9]+", "_", full_name)
  list(
    label   = label,
    summary = file.path(summaries_dir, summary_file),
    dist    = file.path(results_dir, folder, "distanceFiles",
                        paste0(folder, "_all_ptm_distances.csv"))
  )
}

complexes <- list(
  make_complex_entry(
    "26S-proteasome",
    "26S-proteasome_ComplexoFinderSummary.csv",
    "26S proteasome complex"
  ),
  make_complex_entry(
    "Brain-SWI-SNF",
    "Brain-SWI-SNF_ComplexoFinderSummary.csv",
    "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant"
  ),
  make_complex_entry(
    "Dynein-1-complex-variant-1",
    "Dynein-1-complex-variant-1_ComplexoFinderSummary.csv",
    "Dynein-1 complex, variant 1"
  ),
  make_complex_entry(
    "eIF3",
    "Eukaryotic-translation-initiation-factor-3-complex_ComplexoFinderSummary.csv",
    "Eukaryotic translation initiation factor 3 complex"
  )
)

################################################################################
## Generate heatmaps ###########################################################
################################################################################

for (cx in complexes) {
  if (!file.exists(cx$dist) || !file.exists(cx$summary)) {
    message("Skipping ", cx$label, " (file not found)")
    next
  }

  for (mode in c("threshold", "highest")) {
    for (col_ord in c("clustered", "by_assignment", "both_by_assignment", "by_protein")) {

      col_suffix <- switch(col_ord,
        clustered          = "",
        by_assignment      = "_rowclust",
        both_by_assignment = "_assignclust",
        by_protein         = "_protein"
      )

      ht <- plot_combined_distance_heatmap(
        distance_csv    = cx$dist,
        cluster_summary = cx$summary,
        mode            = mode,
        col_order       = col_ord,
        threshold       = 0.5,
        title           = cx$label
      )

      if (!is.null(ht)) {
        pdf_path <- file.path(
          output_dir,
          paste0(cx$label, "_combined_", mode, col_suffix, ".pdf")
        )
        pdf(pdf_path, width = 8, height = 7)
        draw(ht)
        dev.off()
        message("Written: ", pdf_path)
      }
    }
  }
}
