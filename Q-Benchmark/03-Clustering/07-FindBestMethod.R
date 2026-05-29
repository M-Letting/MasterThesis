# Libraries
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(scico)
library(cowplot)

# Set working directory
setwd(paste0(here::here(), "/Q-Benchmark/03-Clustering"))

# Set output directory for plots
plot_output_dir <- "00-Plots"
file_output_dir <- "00-Data/07-FindBestMethod"

if (!dir.exists(file_output_dir)) {
  dir.create(file_output_dir, recursive = TRUE)
}
if (!dir.exists(plot_output_dir)) {
  dir.create(plot_output_dir, recursive = TRUE)
}

# Start overall timer
start_time <- Sys.time()

################################################################################
## Calculation of combined comparison metrics for each datatype and mode #######
################################################################################

transform_metric <- function(x, transform) {
  x <- as.numeric(x)

  if (transform == "identity") {
    return(x)
  }

  if (transform == "log10") {
    x[x <= 0] <- NA_real_
    return(log10(x))
  }

  stop("Unknown transform: ", transform)
}

# Function to calculate combined score for a datatype and mode
# (ex. "Complex", "fixed")
# Input: data.table with metric specifications:
# Metric, file_path, value_col, transform, decrease_is_improvement
# Output: List with two data.tables: combined_dt and summary_dt
calculate_combined_score <- function(metric_specs) {
  combined_dt <- NULL
  merge_cols <- NULL

  # Read in all comparison results and combine into one data.table
  for (i in seq_len(nrow(metric_specs))) {
    name <- metric_specs$Metric[i]
    file_path <- metric_specs$file_path[i]
    value_col <- metric_specs$value_col[i]
    transform <- metric_specs$transform[i]
    decrease_is_improvement <- metric_specs$decrease_is_improvement[i]

    dt <- fread(file_path)

    # Keep only relevant columns
    dt <- dt[, unique(c(names(dt)[1:4], value_col)), with = FALSE]

    # Find columns (Dataset, Complex/Accession, Method, ConstraintThreshold)
    merge_cols_current <- c(
      "Dataset",
      "Complex",
      "Method",
      "ConstraintThreshold"
    )

    # If merge columns are not present, try with Accession instead of Complex
    if (!all(merge_cols_current %in% names(dt))) {
      merge_cols_current <- c(
        "Dataset",
        "Accession",
        "Method",
        "ConstraintThreshold"
      )
    }

    if (is.null(merge_cols)) {
      merge_cols <- merge_cols_current
    }

    # If name contains SSE, aggregate per method at dataset/accession level
    if (grepl("SSE", name)) {
      if (!value_col %in% names(dt)) {
        stop("Column '", value_col, "' not found in ", file_path)
      }

      agg_fun <- mean

      dt <- dt[,
        stats::setNames(list(agg_fun(get(value_col), na.rm = TRUE)), value_col),
        by = merge_cols_current
      ]
    }

    # Harmonize dataset labels (e.g. 04_HighNA -> HighNA)
    dt[, Dataset := gsub(".*_", "", Dataset)]

    # Rename metric column to the name of the datatype
    setnames(dt, old = value_col, new = name)
    dt[, (name) := as.numeric(get(name))]

    transformed_name <- paste0("Transformed_", name)
    dt[, (transformed_name) := transform_metric(get(name), transform)]

    # For each metric, calculate change from baseline (BlindDCF)
    id_col <- if ("Complex" %in% names(dt)) "Complex" else "Accession"

    baseline_name <- paste0("BaselineTransformed_", name)
    baseline_dt <- dt[
      Method == "BlindDCF" & is.na(ConstraintThreshold),
      c("Dataset", id_col, transformed_name),
      with = FALSE
    ]

    setnames(baseline_dt, old = transformed_name, new = baseline_name)
    baseline_dt <- unique(baseline_dt, by = c("Dataset", id_col))

    dt <- merge(
      dt,
      baseline_dt,
      by = c("Dataset", id_col),
      all.x = TRUE,
      sort = FALSE
    )

    comparison_name <- paste0("Change_", name)

    # Compare in transformed metric space to align with 06-PlotComparison
    dt[, (comparison_name) := get(transformed_name) - get(baseline_name)]

    # Reverse sign of comparison if decrease is improvement for this datatype
    if (decrease_is_improvement) {
      dt[, (comparison_name) := -get(comparison_name)]
    }

    # Combine with previous data.tables
    if (is.null(combined_dt)) {
      combined_dt <- dt
    } else {
      combined_dt <- merge(
        combined_dt,
        dt,
        by = merge_cols
      )
    }
  }

  # Get change columns
  change_cols <- grep("Change_", names(combined_dt), value = TRUE)

  # Calculate combined score as mean of baseline-relative changes
  combined_dt[,
    CombinedScore := rowMeans(.SD, na.rm = TRUE),
    .SDcols = change_cols
  ]
  combined_dt[is.nan(CombinedScore), CombinedScore := NA_real_]

  # Create summary table with mean combined score for each method and threshold
  summary_dt <- combined_dt[,
    .(MeanCombinedScore = mean(CombinedScore, na.rm = TRUE)),
    by = c("Method", "ConstraintThreshold")
  ]

  for (metric_name in metric_specs$Metric) {
    change_col <- paste0("Change_", metric_name)
    transformed_col <- paste0("Transformed_", metric_name)
    baseline_col <- paste0("BaselineTransformed_", metric_name)

    metric_summary <- combined_dt[,
      .(
        ChangeFromBaseline = mean(get(change_col), na.rm = TRUE),
        MethodValue = mean(get(transformed_col), na.rm = TRUE),
        BaselineValue = mean(get(baseline_col), na.rm = TRUE)
      ),
      by = c("Method", "ConstraintThreshold")
    ]

    setnames(
      metric_summary,
      old = c("ChangeFromBaseline", "MethodValue", "BaselineValue"),
      new = c(
        paste0("ChangeFromBaseline_", metric_name),
        paste0("MethodValue_", metric_name),
        paste0("BaselineValue_", metric_name)
      )
    )

    summary_dt <- merge(
      summary_dt,
      metric_summary,
      by = c("Method", "ConstraintThreshold")
    )
  }

  # Sort summary table by mean combined score
  summary_dt <- summary_dt[order(-MeanCombinedScore)]

  return(list(combined_dt = combined_dt, summary_dt = summary_dt))
}

##########################################################################
# Run comparison of clustering results for complexes for fixed threshold #
##########################################################################

metric_specs <- data.table(
  Metric = c("Fraction", "WeightedSSE", "Assigned"),
  file_path = c(
    "00-Data/05-CompareClustering/ComplexFractionEvaluation_fixed.csv",
    "00-Data/05-CompareClustering/ComplexWeightedSSEEvaluation_fixed.csv",
    "00-Data/05-CompareClustering/ComplexAssignedEvaluation_fixed.csv"
  ),
  value_col = c("Top_Fraction", "Weighted_SSE", "Assigned_Features"),
  transform = c("identity", "identity", "identity"),
  decrease_is_improvement = c(FALSE, TRUE, FALSE)
)

complex_comparison <- calculate_combined_score(metric_specs)

fwrite(
  complex_comparison$summary_dt,
  file.path(file_output_dir, "ComplexFixedSummaryScores.csv")
)

########################################################################
# Run comparison of clustering results for complexes for 1/k threshold #
########################################################################

metric_specs <- data.table(
  Metric = c("Fraction", "WeightedSSE", "Assigned"),
  file_path = c(
    "00-Data/05-CompareClustering/ComplexFractionEvaluation_one_over_k.csv",
    "00-Data/05-CompareClustering/ComplexWeightedSSEEvaluation_one_over_k.csv",
    "00-Data/05-CompareClustering/ComplexAssignedEvaluation_one_over_k.csv"
  ),
  value_col = c("Top_Fraction", "Weighted_SSE", "Assigned_Features"),
  transform = c("identity", "identity", "identity"),
  decrease_is_improvement = c(FALSE, TRUE, FALSE)
)

complex_comparison <- calculate_combined_score(metric_specs)

fwrite(
  complex_comparison$summary_dt,
  file.path(file_output_dir, "ComplexOneOverKSummaryScores.csv")
)

##########################################################################
# Run comparison of clustering for ProteoMaker datasets, fixed threshold #
##########################################################################

metric_specs <- data.table(
  Metric = c(
    "Purity",
    "Fraction",
    "SumUniqueProteoforms",
    "TotalPotentialProteoforms",
    "WeightedSSE",
    "Assigned"
  ),
  file_path = c(
    "00-Data/05-CompareClustering/ProteoMakerClusterPurityEvaluation_fixed.csv",
    "00-Data/05-CompareClustering/ProteoMakerFractionEvaluation_fixed.csv",
    "00-Data/05-CompareClustering/ProteoMakerSumUniqueProteoforms_fixed.csv",
    "00-Data/05-CompareClustering/ProteoMakerTotalPotentialProteoforms_fixed.csv",
    "00-Data/05-CompareClustering/ProteoMakerWeightedSSEEvaluation_fixed.csv",
    "00-Data/05-CompareClustering/ProteoMakerAssignedEvaluation_fixed.csv"
  ),
  value_col = c(
    "Weighted_Mean_Cluster_Purity",
    "Top_Fraction",
    "Sum_Unique_Proteoforms",
    "Sum_Potential_Proteoforms",
    "Weighted_SSE",
    "Assigned_Features"
  ),
  transform = c(
    "identity",
    "identity",
    "identity",
    "log10",
    "identity",
    "identity"
  ),
  decrease_is_improvement = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE)
)

proteomaker_comparison <- calculate_combined_score(metric_specs)

fwrite(
  proteomaker_comparison$summary_dt,
  file.path(file_output_dir, "ProteoMakerFixedSummaryScores.csv")
)

########################################################################
# Run comparison of clustering for ProteoMaker datasets, 1/k threshold #
########################################################################

metric_specs <- data.table(
  Metric = c(
    "Purity",
    "Fraction",
    "SumUniqueProteoforms",
    "TotalPotentialProteoforms",
    "WeightedSSE",
    "Assigned"
  ),
  file_path = c(
    "00-Data/05-CompareClustering/ProteoMakerClusterPurityEvaluation_one_over_k.csv",
    "00-Data/05-CompareClustering/ProteoMakerFractionEvaluation_one_over_k.csv",
    "00-Data/05-CompareClustering/ProteoMakerSumUniqueProteoforms_one_over_k.csv",
    "00-Data/05-CompareClustering/ProteoMakerTotalPotentialProteoforms_one_over_k.csv",
    "00-Data/05-CompareClustering/ProteoMakerWeightedSSEEvaluation_one_over_k.csv",
    "00-Data/05-CompareClustering/ProteoMakerAssignedEvaluation_one_over_k.csv"
  ),
  value_col = c(
    "Weighted_Mean_Cluster_Purity",
    "Top_Fraction",
    "Sum_Unique_Proteoforms",
    "Sum_Potential_Proteoforms",
    "Weighted_SSE",
    "Assigned_Features"
  ),
  transform = c(
    "identity",
    "identity",
    "identity",
    "log10",
    "identity",
    "identity"
  ),
  decrease_is_improvement = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE)
)

proteomaker_comparison <- calculate_combined_score(metric_specs)

fwrite(
  proteomaker_comparison$summary_dt,
  file.path(file_output_dir, "ProteoMakerOneOverKSummaryScores.csv")
)

################################################################################
## Compare methods across datatypes and modes ##################################
################################################################################

# Get average rank of each method + threshold across datatypes and modes
summary_files <- list.files(
  file_output_dir,
  pattern = "SummaryScores.csv",
  full.names = TRUE
)

# Read in summary files and calculate average rank of each method + threshold
all_summaries <- rbindlist(
  lapply(summary_files, function(file) {
    dt <- fread(file)
    dt[, Dataset := gsub("SummaryScores.csv", "", basename(file))]
    return(dt)
  }),
  fill = TRUE
)

all_summaries[,
  Rank := rank(-MeanCombinedScore, ties.method = "average"),
  by = .(Dataset)
]

average_ranks <- all_summaries[,
  .(AverageRank = mean(Rank)),
  by = .(Method, ConstraintThreshold)
]

average_ranks <- average_ranks[order(AverageRank)]

# Get method with average rank better than baseline (BlindDCF with no constraints)
baseline_rank <- average_ranks[
  Method == "BlindDCF" & is.na(ConstraintThreshold),
  AverageRank
]

average_ranks <- average_ranks[average_ranks$AverageRank < baseline_rank]

# Add all change from baseline values for each metric and dataset for these methods
best_methods <- merge(
  average_ranks,
  all_summaries,
  by = c("Method", "ConstraintThreshold")
)

# Reorder columns
best_methods <- best_methods[,
  c(
    "Method",
    "ConstraintThreshold",
    "Dataset",
    "MeanCombinedScore",
    "AverageRank",
    grep("ChangeFromBaseline_", names(best_methods), value = TRUE),
    grep("MethodValue_", names(best_methods), value = TRUE),
    grep("BaselineValue_", names(best_methods), value = TRUE)
  ),
  with = FALSE
]

fwrite(
  best_methods,
  file.path(file_output_dir, "BestMethodsSummary.csv")
)

################################################################################
## Plot comparison of best methods across datatypes and modes ##################
################################################################################

dt <- copy(best_methods)

dt[,
  MethodThreshold := fifelse(
    is.na(ConstraintThreshold),
    Method,
    sprintf("%s | %.2f", Method, ConstraintThreshold)
  )
]

# Define metric direction
higher_better <- c("Purity", "Fraction", "Assigned")
lower_better <- c(
  "SumUniqueProteoforms",
  "TotalPotentialProteoforms",
  "WeightedSSE"
)

direction_dt <- rbindlist(list(
  data.table(Metric = higher_better, HigherBetter = TRUE),
  data.table(Metric = lower_better, HigherBetter = FALSE)
))

# Keep top method-threshold combos
top_n <- min(10L, uniqueN(dt$MethodThreshold))
top_levels <- dt[order(AverageRank), unique(MethodThreshold)][1:top_n]
dt <- dt[MethodThreshold %in% top_levels]

# Build metric-wise table from MethodValue_* and BaselineValue_*
method_cols <- grep("^MethodValue_", names(dt), value = TRUE)
metrics <- sub("^MethodValue_", "", method_cols)
baseline_cols <- paste0("BaselineValue_", metrics)

valid <- baseline_cols %in% names(dt)
metrics <- metrics[valid]
method_cols <- method_cols[valid]
baseline_cols <- baseline_cols[valid]

id_cols <- intersect(
  c(
    "Method",
    "ConstraintThreshold",
    "Dataset",
    "MethodThreshold",
    "MeanCombinedScore",
    "AverageRank"
  ),
  names(dt)
)

improvement_dt <- rbindlist(
  lapply(seq_along(metrics), function(i) {
    metric <- metrics[i]
    mcol <- method_cols[i]
    bcol <- baseline_cols[i]

    tmp <- dt[, c(id_cols, mcol, bcol), with = FALSE]
    setnames(tmp, c(mcol, bcol), c("MethodValue", "BaselineValue"))
    tmp[, Metric := metric]
    tmp
  }),
  use.names = TRUE,
  fill = TRUE
)

improvement_dt <- merge(
  improvement_dt,
  direction_dt,
  by = "Metric",
  all.x = TRUE
)

# Safety check: fail fast if a metric direction is missing
if (improvement_dt[is.na(HigherBetter), .N] > 0) {
  stop("Missing direction (higher/lower better) for one or more metrics.")
}

improvement_dt[,
  Improvement := fifelse(
    HigherBetter,
    MethodValue - BaselineValue,
    BaselineValue - MethodValue
  )
]

improvement_dt <- improvement_dt[is.finite(Improvement)]

if (nrow(improvement_dt) == 0L) {
  warning(
    "No valid results forund; skipping improvement plot. ",
    "Re-run 05-CompareClustering.R / 07-RunComparison.R to regenerate CSVs."
  )
} else {
  score_order <- dt[,
    .(Rank = mean(AverageRank, na.rm = TRUE)),
    by = MethodThreshold
  ][order(Rank), MethodThreshold]

  improvement_dt[,
    MethodThreshold := factor(MethodThreshold, levels = score_order)
  ]

  p_improvement <- ggplot(
    improvement_dt,
    aes(
      x = MethodThreshold,
      y = Improvement,
      color = Method,
      shape = Dataset
    )
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    geom_point(
      size = 2.2,
      alpha = 0.8,
      position = position_jitter(width = 0.12, height = 0)
    ) +
    stat_summary(
      aes(x = MethodThreshold, y = Improvement, group = MethodThreshold),
      fun = mean,
      geom = "point",
      inherit.aes = FALSE,
      shape = 18,
      size = 1.5,
      color = "black"
    ) +
    facet_wrap(~Metric, scales = "free_y") +
    scale_color_scico_d(palette = "batlow") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = "right"
    ) +
    labs(
      y = "Signed improvement (> 0 is better)",
      color = "Method",
      shape = "Dataset",
      title = "Improvement vs Baseline by Metric"
    )

  ggsave(
    p_improvement,
    filename = file.path(plot_output_dir, "07-BestMethodImprovements.pdf"),
    width = 12,
    height = 6.5
  )
}

################################################################################
## Number of restrictions applied by each constraint method and threshold ######
################################################################################

# Extract per-dataset restriction counts from clustering result RDS files.
# PC files: threshold -> dataset -> list(vsclust_result, constraint_matrix)
# PM files: threshold -> dataset -> protein -> list(vsclust_result, constraint_matrix)
extract_restrictions_from_rds <- function(rds_files, is_pm) {
  result_list <- list()

  for (f in rds_files) {
    clm_match <- regmatches(
      basename(f),
      regexpr("_(ccc|corr|rmsd)\\.rds$", basename(f))
    )
    if (length(clm_match) == 0) {
      next
    }
    cannot_link_method <- gsub("^_|\\.rds$", "", clm_match)

    results <- readRDS(f)

    for (th in names(results)) {
      th_num <- as.numeric(th) / 10
      th_results <- results[[th]]

      for (ds in names(th_results)) {
        if (is_pm) {
          for (prot in names(th_results[[ds]])) {
            entry <- th_results[[ds]][[prot]]
            if (!is.null(entry$constraint_matrix)) {
              result_list[[length(result_list) + 1L]] <- data.table(
                CannotLinkMethod = cannot_link_method,
                DataType = "ProteoMaker",
                Dataset = ds,
                Protein = prot,
                Threshold = th_num,
                nRestrictions = sum(entry$constraint_matrix, na.rm = TRUE)
              )
            }
          }
        } else {
          entry <- th_results[[ds]]
          if (!is.null(entry$constraint_matrix)) {
            result_list[[length(result_list) + 1L]] <- data.table(
              CannotLinkMethod = cannot_link_method,
              DataType = "Complex",
              Dataset = ds,
              Protein = "All",
              Threshold = th_num,
              nRestrictions = sum(entry$constraint_matrix, na.rm = TRUE)
            )
          }
        }
      }
    }
  }

  if (length(result_list) == 0L) {
    return(NULL)
  }
  rbindlist(result_list)
}

summarise_restrictions <- function(dt) {
  dt[,
    .(
      TotalRestrictions = sum(nRestrictions, na.rm = TRUE),
      MeanRestrictionsPerProtein = mean(nRestrictions, na.rm = TRUE),
      MedianRestrictionsPerProtein = median(nRestrictions, na.rm = TRUE),
      NProteins = .N
    ),
    by = .(DataType, CannotLinkMethod, Threshold)
  ][order(DataType, CannotLinkMethod, Threshold)]
}

# Count violated cannot-link constraints per entry.
# A violation occurs when a feature is assigned to a cluster that is marked
# TRUE in its constraint_matrix row (i.e. it was constrained away from that cluster).
count_violations <- function(entry) {
  cm <- entry$constraint_matrix
  cl_obj <- entry$vsclust_result$ClustOut$Bestcl
  if (is.null(cm) || is.null(cl_obj)) {
    return(NA_integer_)
  }
  cluster_assign <- cl_obj$cluster
  common <- intersect(rownames(cm), names(cluster_assign))
  if (length(common) == 0L) {
    return(0L)
  }
  cm_sub <- cm[common, , drop = FALSE]
  cl_sub <- cluster_assign[common]
  col_names <- colnames(cm_sub)
  n_violated <- sum(
    mapply(
      function(feat, cl) {
        col <- as.character(cl)
        if (col %in% col_names) isTRUE(cm_sub[feat, col]) else FALSE
      },
      common,
      cl_sub
    ),
    na.rm = TRUE
  )
  as.integer(n_violated)
}

extract_violations_from_rds <- function(rds_files, is_pm) {
  result_list <- list()

  for (f in rds_files) {
    clm_match <- regmatches(
      basename(f),
      regexpr("_(ccc|corr|rmsd)\\.rds$", basename(f))
    )
    if (length(clm_match) == 0) {
      next
    }
    cannot_link_method <- gsub("^_|\\.rds$", "", clm_match)

    results <- readRDS(f)

    for (th in names(results)) {
      th_num <- as.numeric(th) / 10
      th_results <- results[[th]]

      for (ds in names(th_results)) {
        if (is_pm) {
          for (prot in names(th_results[[ds]])) {
            entry <- th_results[[ds]][[prot]]
            nv <- count_violations(entry)
            if (!is.na(nv)) {
              result_list[[length(result_list) + 1L]] <- data.table(
                CannotLinkMethod = cannot_link_method,
                DataType = "ProteoMaker",
                Dataset = ds,
                Protein = prot,
                Threshold = th_num,
                nViolations = nv
              )
            }
          }
        } else {
          entry <- th_results[[ds]]
          nv <- count_violations(entry)
          if (!is.na(nv)) {
            result_list[[length(result_list) + 1L]] <- data.table(
              CannotLinkMethod = cannot_link_method,
              DataType = "Complex",
              Dataset = ds,
              Protein = "All",
              Threshold = th_num,
              nViolations = nv
            )
          }
        }
      }
    }
  }

  if (length(result_list) == 0L) {
    return(NULL)
  }
  rbindlist(result_list)
}

summarise_violations <- function(dt) {
  dt[,
    .(
      TotalViolations = sum(nViolations, na.rm = TRUE),
      MeanViolationsPerProtein = mean(nViolations, na.rm = TRUE),
      MedianViolationsPerProtein = median(nViolations, na.rm = TRUE),
      NProteins = .N
    ),
    by = .(DataType, CannotLinkMethod, Threshold)
  ][order(DataType, CannotLinkMethod, Threshold)]
}

# DCF-derived constraints (dCF_to_constraints applied to cannot-link matrices)
dcf_dt <- rbindlist(
  list(
    extract_restrictions_from_rds(
      list.files(
        "00-Data/03-DCFConstraints",
        pattern = "^01-PC_DCF_ClusteringResults_(ccc|corr|rmsd)\\.rds$",
        full.names = TRUE
      ),
      is_pm = FALSE
    ),
    extract_restrictions_from_rds(
      list.files(
        "00-Data/03-DCFConstraints",
        pattern = "^01-PM_DCF_ClusteringResults_(ccc|corr|rmsd)\\.rds$",
        full.names = TRUE
      ),
      is_pm = TRUE
    )
  ),
  fill = TRUE
)

dcf_summary <- summarise_restrictions(dcf_dt)

fwrite(
  dcf_summary,
  file.path(file_output_dir, "RestrictionCountSummary_DCF.csv")
)

# VSClust-derived constraints (vsclust_to_restrictions applied to cannot-link matrices)
vsclust_dt <- rbindlist(
  list(
    extract_restrictions_from_rds(
      list.files(
        "00-Data/04-VSClustConstraints",
        pattern = "^01-PC_DCF_VSClust_ClusteringResults_(ccc|corr|rmsd)\\.rds$",
        full.names = TRUE
      ),
      is_pm = FALSE
    ),
    extract_restrictions_from_rds(
      list.files(
        "00-Data/04-VSClustConstraints",
        pattern = "^01-PM_DCF_VSClust_ClusteringResults_(ccc|corr|rmsd)\\.rds$",
        full.names = TRUE
      ),
      is_pm = TRUE
    )
  ),
  fill = TRUE
)

vsclust_summary <- summarise_restrictions(vsclust_dt)

fwrite(
  vsclust_summary,
  file.path(file_output_dir, "RestrictionCountSummary_VSClust.csv")
)

cat("\nDCF-derived restriction counts per cannot-link method and threshold:\n")
print(dcf_summary)
cat(
  "\nVSClust-derived restriction counts per cannot-link method and threshold:\n"
)
print(vsclust_summary)

################################################################################
## Number of violated constraints per constraint method and threshold ##########
################################################################################

dcf_viol_dt <- rbindlist(
  list(
    extract_violations_from_rds(
      list.files(
        "00-Data/03-DCFConstraints",
        pattern = "^01-PC_DCF_ClusteringResults_(ccc|corr|rmsd)\\.rds$",
        full.names = TRUE
      ),
      is_pm = FALSE
    ),
    extract_violations_from_rds(
      list.files(
        "00-Data/03-DCFConstraints",
        pattern = "^01-PM_DCF_ClusteringResults_(ccc|corr|rmsd)\\.rds$",
        full.names = TRUE
      ),
      is_pm = TRUE
    )
  ),
  fill = TRUE
)

dcf_viol_summary <- summarise_violations(dcf_viol_dt)

fwrite(
  dcf_viol_summary,
  file.path(file_output_dir, "ViolationCountSummary_DCF.csv")
)

vsclust_viol_dt <- rbindlist(
  list(
    extract_violations_from_rds(
      list.files(
        "00-Data/04-VSClustConstraints",
        pattern = "^01-PC_DCF_VSClust_ClusteringResults_(ccc|corr|rmsd)\\.rds$",
        full.names = TRUE
      ),
      is_pm = FALSE
    ),
    extract_violations_from_rds(
      list.files(
        "00-Data/04-VSClustConstraints",
        pattern = "^01-PM_DCF_VSClust_ClusteringResults_(ccc|corr|rmsd)\\.rds$",
        full.names = TRUE
      ),
      is_pm = TRUE
    )
  ),
  fill = TRUE
)

vsclust_viol_summary <- summarise_violations(vsclust_viol_dt)

fwrite(
  vsclust_viol_summary,
  file.path(file_output_dir, "ViolationCountSummary_VSClust.csv")
)

cat("\nDCF-derived violation counts per cannot-link method and threshold:\n")
print(dcf_viol_summary)
cat(
  "\nVSClust-derived violation counts per cannot-link method and threshold:\n"
)
print(vsclust_viol_summary)

################################################################################
## Plot number of restrictions applied by each constraint method and threshold #
################################################################################

top_n_rest <- min(10L, nrow(average_ranks))
top_methods_rest <- average_ranks[1:top_n_rest]

top_methods_rest[,
  CannotLinkMethod := fcase(
    grepl("_CCC$", Method)  , "ccc"  ,
    grepl("_Corr$", Method) , "corr" ,
    grepl("_RMSD$", Method) , "rmsd" ,
    default = NA_character_
  )
]
top_methods_rest[,
  ConstraintSource := fcase(
    grepl("^DCFDCF", Method)  , "DCF"     ,
    grepl("^VSClust", Method) , "VSClust" ,
    default = NA_character_
  )
]
top_methods_rest[,
  MethodThreshold := fifelse(
    is.na(ConstraintThreshold),
    Method,
    sprintf("%s | %.2f", Method, ConstraintThreshold)
  )
]

rest_totals <- rbindlist(list(
  fread(file.path(file_output_dir, "RestrictionCountSummary_DCF.csv"))[,
    ConstraintSource := "DCF"
  ],
  fread(file.path(file_output_dir, "RestrictionCountSummary_VSClust.csv"))[,
    ConstraintSource := "VSClust"
  ]
))[,
  .(TotalRestrictions = sum(TotalRestrictions, na.rm = TRUE)),
  by = .(CannotLinkMethod, Threshold, ConstraintSource)
]

rest_dt <- merge(
  top_methods_rest,
  rest_totals,
  by.x = c("CannotLinkMethod", "ConstraintThreshold", "ConstraintSource"),
  by.y = c("CannotLinkMethod", "Threshold", "ConstraintSource"),
  all.x = TRUE
)
rest_dt[is.na(TotalRestrictions), TotalRestrictions := 0L]
rest_dt[,
  MethodThreshold := factor(MethodThreshold, levels = MethodThreshold)
]

p_restrictions <- ggplot(
  rest_dt,
  aes(x = MethodThreshold, y = TotalRestrictions, fill = Method)
) +
  geom_col() +
  scale_fill_scico_d(palette = "batlow") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    y = "Total restrictions",
    fill = "Method",
    title = "Number of Applied Restrictions by Method"
  )

ggsave(
  p_restrictions,
  filename = file.path(plot_output_dir, "07-BestMethodRestrictions.pdf"),
  width = 8,
  height = 4
)

################################################################################
## Plot number of violated constraints by each constraint method and threshold #
################################################################################

viol_totals <- rbindlist(list(
  fread(file.path(file_output_dir, "ViolationCountSummary_DCF.csv"))[,
    ConstraintSource := "DCF"
  ],
  fread(file.path(file_output_dir, "ViolationCountSummary_VSClust.csv"))[,
    ConstraintSource := "VSClust"
  ]
))[,
  .(TotalViolations = sum(TotalViolations, na.rm = TRUE)),
  by = .(CannotLinkMethod, Threshold, ConstraintSource)
]

viol_dt <- merge(
  top_methods_rest,
  viol_totals,
  by.x = c("CannotLinkMethod", "ConstraintThreshold", "ConstraintSource"),
  by.y = c("CannotLinkMethod", "Threshold", "ConstraintSource"),
  all.x = TRUE
)
viol_dt[is.na(TotalViolations), TotalViolations := 0L]
viol_dt[,
  MethodThreshold := factor(
    MethodThreshold,
    levels = levels(rest_dt$MethodThreshold)
  )
]

p_violations <- ggplot(
  viol_dt,
  aes(x = MethodThreshold, y = TotalViolations, fill = Method)
) +
  geom_col() +
  scale_fill_scico_d(palette = "batlow") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    y = "Total violated constraints",
    fill = "Method",
    title = "Number of Violated Constraints by Method"
  )

ggsave(
  p_violations,
  filename = file.path(plot_output_dir, "07-BestMethodViolations.pdf"),
  width = 8,
  height = 4
)

################################################################################
# End timer
end_time <- Sys.time()
total_time <- round(end_time - start_time, 2)
cat("Total time taken for clustering comparison:", total_time, "seconds \n")
cat(paste(rep("#", 80), collapse = ""), "\n")
