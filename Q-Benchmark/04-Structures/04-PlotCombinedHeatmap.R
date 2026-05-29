# Libraries
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(scico)

setwd(here::here())

source("ComplexoVisualizer/MakeCXC.R")

output_dir <- "Q-Benchmark/04-Structures/00-Plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

################################################################################
## Core plotting function #######################################################
################################################################################

#' Build a combined distance heatmap for all PTM features across clusters
#'
#' @param distance_csv  Path to the `_all_ptm_distances.csv` produced by
#'   `compute_ptm_distances_combined()` or its AF-subunit equivalent.
#' @param cluster_summary  Path to the cluster summary CSV (fuzzy membership
#'   columns named "membership of cluster N").
#' @param mode  `"threshold"`: include only features where max membership >=
#'   `threshold`; `"highest"`: include all features, assign each to its
#'   highest-membership cluster regardless of value.
#' @param threshold  Numeric cutoff used when `mode = "threshold"` (default 0.5).
#' @param title  Character string used as the heatmap column title.
#' @param cluster_method  Linkage method passed to `hclust` (default "complete").
#'
#' @return A `Heatmap` object (invisibly NULL if fewer than 2 features remain).
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
    # Hard ClusterID assignment (non-fuzzy)
    cs_ptm[, max_mem   := 1.0]
    cs_ptm[, cl_assign := ClusterID]
    k <- length(unique(cs_ptm$ClusterID))
  }

  if (mode == "threshold") cs_ptm <- cs_ptm[max_mem >= threshold]

  if (nrow(cs_ptm) < 2L) {
    warning("Fewer than 2 features after filtering (mode='", mode, "') — skipping.")
    return(invisible(NULL))
  }

  # Composite key: Position + Datatype — unique even when two modification types
  # share the same peptide position (e.g. RmCys and FreeCys at [25-34]).
  cs_ptm[, feat_id := paste0(Position, "\n", Datatype)]
  feats <- cs_ptm$feat_id
  n     <- length(feats)

  # ── Build n×n distance matrix ─────────────────────────────────────────────
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

  # ── Drop features with no measured distances (subunit absent from structure) ─
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

  # ── Hierarchical clustering (used when either axis is "clustered") ────────
  max_d      <- max(dist_mat, na.rm = TRUE)
  dist_clean <- dist_mat
  dist_clean[is.na(dist_clean)] <- max_d * 1.5
  hc      <- hclust(as.dist(dist_clean), method = cluster_method)
  hc_ord  <- hc$order

  assign_ord  <- order(cs_ptm$cl_assign, na.last = TRUE)
  protein_ord <- order(cs_ptm[["Gene name"]], na.last = TRUE)

  row_ord <- switch(col_order,
    clustered          = hc_ord,
    by_assignment      = hc_ord,
    both_by_assignment = assign_ord,
    by_protein         = protein_ord
  )
  col_ord <- switch(col_order,
    clustered          = hc_ord,
    by_assignment      = assign_ord,
    both_by_assignment = assign_ord,
    by_protein         = protein_ord
  )

  row_feats <- feats[row_ord]
  col_feats <- feats[col_ord]
  ordered_mat   <- dist_mat[row_feats, col_feats]
  row_cl        <- as.character(cs_ptm$cl_assign[row_ord])
  col_cl        <- as.character(cs_ptm$cl_assign[col_ord])
  row_genes     <- cs_ptm[["Gene name"]][row_ord]
  col_genes     <- cs_ptm[["Gene name"]][col_ord]
  row_dtypes    <- cs_ptm$Datatype[row_ord]
  col_dtypes    <- cs_ptm$Datatype[col_ord]

  # Cluster color palette (same as CXC scripts)
  cl_colors        <- cluster_palette(k)
  names(cl_colors) <- as.character(seq_len(k))
  all_cl <- c(row_cl, col_cl)
  if (any(is.na(all_cl))) cl_colors <- c(cl_colors, "NA" = "#CCCCCC")

  # Gene color palette — build from all genes so colors are consistent
  unique_genes <- unique(c(row_genes, col_genes))
  gene_colors  <- setNames(
    scico(length(unique_genes), palette = "batlow", end = 0.9),
    unique_genes
  )

  # Datatype color palette — subset to only the types present in this heatmap
  dt_pal       <- datatype_palette()
  present_dts  <- unique(c(row_dtypes, col_dtypes))
  dt_colors    <- dt_pal[intersect(names(dt_pal), present_dts)]
  unknown_dts  <- setdiff(present_dts, names(dt_colors))
  if (length(unknown_dts) > 0L) dt_colors[unknown_dts] <- "#CCCCCC"

  ann_col_list <- list(Cluster = cl_colors, Gene = gene_colors, Datatype = dt_colors)

  # Top annotation (columns)
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

  # Left annotation (rows) — legends suppressed to avoid duplication
  left_ann <- rowAnnotation(
    Cluster  = row_cl,
    Gene     = row_genes,
    Datatype = row_dtypes,
    col      = ann_col_list,
    show_legend          = FALSE,
    annotation_name_side = "bottom"
  )

  # Distance color scale
  col_fun <- colorRamp2(
    c(0, max_d),
    c("#590007", "#EBE5E0")
  )

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
## Input file lists #############################################################
################################################################################

CF_prefix <- "Q-Benchmark/04-Structures/00-Data/Results/ComplexoFinder"
PF_prefix <- "Q-Benchmark/04-Structures/00-Data/Results/ProteoForge"

complexes <- list(
  list(
    label   = "26S-proteasome",
    summary = "Q-Benchmark/04-Structures/00-Data/Summaries/26S-proteasome_ComplexoFinderSummary.csv",
    cf_dist = file.path(
      CF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_all_ptm_distances.csv"
    ),
    pf_summary = "Q-Benchmark/04-Structures/00-Data/Summaries/26S-proteasome_ProteoForgeSummary.csv",
    pf_dist = file.path(
      PF_prefix,
      "26S_proteasome_complex/distanceFiles/26S_proteasome_complex_all_ptm_distances.csv"
    )
  ),
  list(
    label   = "Brain-SWI-SNF",
    summary = "Q-Benchmark/04-Structures/00-Data/Summaries/Brain-SWI-SNF_ComplexoFinderSummary.csv",
    cf_dist = file.path(
      CF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_all_ptm_distances.csv"
    ),
    pf_summary = "Q-Benchmark/04-Structures/00-Data/Summaries/Brain-SWI-SNF_ProteoForgeSummary.csv",
    pf_dist = file.path(
      PF_prefix,
      "Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant/distanceFiles/Brain_specific_SWI_SNF_ATP_dependent_chromatin_remodeling_complex_ARID1A_SMARCA2_variant_all_ptm_distances.csv"
    )
  ),
  list(
    label   = "Dynein-1-complex-variant-1",
    summary = "Q-Benchmark/04-Structures/00-Data/Summaries/Dynein-1-complex-variant-1_ComplexoFinderSummary.csv",
    cf_dist = file.path(
      CF_prefix,
      "Dynein_1_complex_variant_4/distanceFiles/Dynein_1_complex_variant_4_all_ptm_distances.csv"
    ),
    pf_summary = "Q-Benchmark/04-Structures/00-Data/Summaries/Dynein-1-complex-variant-1_ProteoForgeSummary.csv",
    pf_dist = file.path(
      PF_prefix,
      "Dynein_1_complex_variant_4/distanceFiles/Dynein_1_complex_variant_4_all_ptm_distances.csv"
    )
  ),
  list(
    label   = "eITF3",
    summary = "Q-Benchmark/04-Structures/00-Data/Summaries/Eukaryotic-translation-initiation-factor-3-complex_ComplexoFinderSummary.csv",
    cf_dist = file.path(
      CF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_all_ptm_distances.csv"
    ),
    pf_summary = "Q-Benchmark/04-Structures/00-Data/Summaries/Eukaryotic-translation-initiation-factor-3-complex_ProteoForgeSummary.csv",
    pf_dist = file.path(
      PF_prefix,
      "Eukaryotic_translation_initiation_factor_3_complex/distanceFiles/Eukaryotic_translation_initiation_factor_3_complex_all_ptm_distances.csv"
    )
  )
)

################################################################################
## Generate heatmaps ############################################################
################################################################################

for (cx in complexes) {
  for (mode in c("threshold", "highest")) {
    for (col_ord in c("clustered", "by_assignment", "both_by_assignment", "by_protein")) {

      col_suffix <- switch(col_ord,
        clustered          = "",
        by_assignment      = "_rowclust",
        both_by_assignment = "_assignclust",
        by_protein         = "_protein"
      )

      # ── ComplexoFinder ────────────────────────────────────────────────────
      if (file.exists(cx$cf_dist) && file.exists(cx$summary)) {
        ht <- plot_combined_distance_heatmap(
          distance_csv    = cx$cf_dist,
          cluster_summary = cx$summary,
          mode            = mode,
          col_order       = col_ord,
          threshold       = 0.5,
          title           = paste0(cx$label, " — ComplexoFinder")
        )
        if (!is.null(ht)) {
          pdf_path <- file.path(
            output_dir,
            paste0(cx$label, "_CF_combined_", mode, col_suffix, ".pdf")
          )
          pdf(pdf_path, width = 8, height = 7)
          draw(ht)
          dev.off()
          message("Written: ", pdf_path)
        }
      } else {
        message("Skipping CF heatmap for ", cx$label, " (file not found)")
      }

      # ── ProteoForge ───────────────────────────────────────────────────────
      if (file.exists(cx$pf_dist) && file.exists(cx$pf_summary)) {
        ht <- plot_combined_distance_heatmap(
          distance_csv    = cx$pf_dist,
          cluster_summary = cx$pf_summary,
          mode            = mode,
          col_order       = col_ord,
          threshold       = 0.5,
          title           = paste0(cx$label, " — ProteoForge")
        )
        if (!is.null(ht)) {
          pdf_path <- file.path(
            output_dir,
            paste0(cx$label, "_PF_combined_", mode, col_suffix, ".pdf")
          )
          pdf(pdf_path, width = 8, height = 7)
          draw(ht)
          dev.off()
          message("Written: ", pdf_path)
        }
      } else {
        message("Skipping PF heatmap for ", cx$label, " (file not found)")
      }
    }
  }
}
