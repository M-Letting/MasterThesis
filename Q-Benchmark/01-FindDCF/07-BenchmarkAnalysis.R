# ============================================================================
# BENCHMARK ANALYSIS - ROC, PR, F1, and MCC Curves
# ============================================================================

cat("\n")
cat("=", rep("=", 78), "\n", sep = "")
cat("           BENCHMARK ANALYSIS - PERFORMANCE EVALUATION            \n")
cat("=", rep("=", 78), "\n", sep = "")
cat("\n")

# Load required libraries
library(arrow)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scico)

# Configuration
result_dir <- "01-FindDCF/data/results"
figure_dir <- "01-FindDCF/figures"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
pthr <- 1e-3 # P-value threshold to highlight
thresholds <- 10^seq(0, -15, length.out = 10000)

# Data scenarios
data_ids <- c("1pep", "2pep", "050pep", "random")
data_names <- c(
  "1pep" = "1 Peptide",
  "2pep" = "2 Peptides",
  "050pep" = "%50 Peptides",
  "random" = "Random (1 to %50) Peptides"
)

# Methods
methods <- c("COPF", "PeCorA", "ProteoForge", "DCF", "DCFBaseline")

# Color palette
color_palette <- scico(5, palette = "batlowK", direction = -1)

method_colors <- c(
  "COPF" = color_palette[1],
  "PeCorA" = color_palette[2],
  "ProteoForge" = color_palette[3],
  "DCF" = color_palette[4],
  "DCFBaseline" = color_palette[5]
)

# Helper function to calculate performance metrics
calculate_metrics <- function(data, thresholds, label_col, pvalue_col) {
  results <- lapply(thresholds, function(thr) {
    # Classify as positive if p-value <= threshold
    predicted_positive <- data[[pvalue_col]] <= thr
    true_positive <- data[[label_col]] == 1 | data[[label_col]] == TRUE

    # Calculate confusion matrix
    TP <- sum(predicted_positive & true_positive, na.rm = TRUE)
    FP <- sum(predicted_positive & !true_positive, na.rm = TRUE)
    TN <- sum(!predicted_positive & !true_positive, na.rm = TRUE)
    FN <- sum(!predicted_positive & true_positive, na.rm = TRUE)

    # Calculate metrics
    TPR <- ifelse(TP + FN > 0, TP / (TP + FN), 0) # Sensitivity/Recall
    FPR <- ifelse(FP + TN > 0, FP / (FP + TN), 0)
    Precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
    Recall <- TPR

    # F1 Score
    F1 <- ifelse(
      Precision + Recall > 0,
      2 * (Precision * Recall) / (Precision + Recall),
      0
    )

    # Matthews Correlation Coefficient
    numerator <- (TP * TN) - (FP * FN)
    denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    MCC <- ifelse(denominator > 0, numerator / denominator, 0)

    data.frame(
      threshold = thr,
      TP = TP,
      FP = FP,
      TN = TN,
      FN = FN,
      TPR = TPR,
      FPR = FPR,
      Precision = Precision,
      Recall = Recall,
      F1 = F1,
      MCC = MCC
    )
  })

  do.call(rbind, results)
}

# Helper function to add complete curve endpoints
complete_curve_data <- function(df, curve_type, x_col, y_col) {
  if (curve_type == "ROC") {
    # Add (0,0) and (1,1) if not present
    if (nrow(df) > 0) {
      if (!any(df[[x_col]] == 0 & df[[y_col]] == 0)) {
        new_row <- df[1, ]
        new_row$FPR <- 0
        new_row$TPR <- 0
        new_row$threshold <- max(df$threshold, na.rm = TRUE) * 10
        df <- rbind(new_row, df)
      }
      if (!any(df[[x_col]] == 1 & df[[y_col]] == 1)) {
        new_row <- df[nrow(df), ]
        new_row$FPR <- 1
        new_row$TPR <- 1
        new_row$threshold <- min(df$threshold, na.rm = TRUE) / 10
        df <- rbind(df, new_row)
      }
    }
  } else if (curve_type == "PR") {
    # For PR curves, add endpoints if needed
    if (nrow(df) > 0) {
      if (!any(df[[x_col]] == 0)) {
        new_row <- df[which.min(df$Recall), ]
        new_row$Recall <- 0
        new_row$threshold <- max(df$threshold, na.rm = TRUE) * 10
        df <- rbind(new_row, df)
      }
      if (!any(df[[x_col]] == 1)) {
        new_row <- df[which.max(df$Recall), ]
        new_row$Recall <- 1
        new_row$threshold <- min(df$threshold, na.rm = TRUE) / 10
        df <- rbind(df, new_row)
      }
    }
  }
  df
}

cat("Step 1: Loading and processing data...\n")

# Load and process all results
results_data <- data.frame()

for (data_id in data_ids) {
  for (method in methods) {
    cat(sprintf("  Processing %s - %s...\n", method, data_id))

    # Determine p-value column and file name
    if (method == "COPF") {
      pval_col <- "proteoform_score_pval"
      file_name <- paste0(method, "_", data_id, "_result.feather")
    } else if (method %in% c("DCF", "DCFBaseline")) {
      pval_col <- "p_adj"
      file_name <- paste0(method, "_", data_id, "_result.feather")
    } else {
      pval_col <- "adj_pval"
      file_name <- paste0(method, "_", data_id, "_result.feather")
    }

    # Read data
    file_path <- file.path(result_dir, file_name)
    data <- read_feather(file_path)

    # For ProteoForge, remove duplicates
    if (method == "ProteoForge") {
      data <- data %>%
        select(protein_id, peptide_id, perturbed_peptide, all_of(pval_col)) %>%
        distinct()
    }

    # Calculate metrics
    metric_data <- calculate_metrics(
      data,
      thresholds,
      label_col = "perturbed_peptide",
      pvalue_col = pval_col
    )

    metric_data$method <- method
    metric_data$perturbation <- data_id

    results_data <- rbind(results_data, metric_data)
  }
}

# Replace perturbation IDs with display names
results_data$perturbation <- data_names[results_data$perturbation]
results_data$perturbation <- factor(
  results_data$perturbation,
  levels = unname(data_names)
)

cat(sprintf("\nData processed: %d rows\n", nrow(results_data)))

# Save processed data
output_file <- file.path(
  result_dir,
  "peptide_identification_performance_data.feather"
)
write_feather(results_data, output_file)
write.csv(results_data, sub(".feather", ".csv", output_file), row.names = FALSE)

cat("\n")
cat("Step 2: Creating ROC curves...\n")

# ============================================================================
# ROC Curves
# ============================================================================

plot_list <- list()

for (i in seq_along(data_names)) {
  pert <- unname(data_names)[i]

  cur_data <- results_data %>%
    filter(perturbation == pert) %>%
    group_by(method) %>%
    group_modify(~ complete_curve_data(.x, "ROC", "FPR", "TPR")) %>%
    ungroup() %>%
    arrange(method, FPR)

  # Calculate AUC for each method
  auc_data <- cur_data %>%
    group_by(method) %>%
    summarise(
      AUC = pracma::trapz(sort(FPR), TPR[order(FPR)]),
      .groups = "drop"
    )

  # Get point at threshold
  pthr_data <- cur_data %>%
    filter(abs(threshold - pthr) == min(abs(threshold - pthr)))

  pthr_005_data <- cur_data %>%
    group_by(method) %>%
    filter(abs(threshold - 0.05) == min(abs(threshold - 0.05))) %>%
    ungroup()

  # Create plot
  p <- ggplot(cur_data, aes(x = FPR, y = TPR, color = method)) +
    geom_line(linewidth = 1) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "black"
    ) +
    geom_point(
      data = pthr_data,
      aes(x = FPR, y = TPR),
      size = 5,
      shape = 21,
      fill = NA,
      stroke = 1.5,
      color = "black"
    ) +
    geom_point(
      data = pthr_005_data,
      aes(x = FPR, y = TPR),
      size = 3,
      shape = 21,
      fill = NA,
      stroke = 1.5,
      color = "red"
    ) +
    scale_color_manual(values = method_colors) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05)) +
    labs(
      title = pert,
      x = "FPR",
      y = "TPR"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.position = "none",
      panel.grid.major = element_line(
        color = "lightgrey",
        linetype = "dashed",
        linewidth = 0.3
      ),
      panel.grid.minor = element_line(
        color = "lightgrey",
        linetype = "dashed",
        linewidth = 0.2
      )
    )

  # Add AUC annotations
  for (j in seq_len(nrow(auc_data))) {
    p <- p +
      annotate(
        "text",
        x = 0.95,
        y = 0.05 + (j - 1) * 0.075,
        label = sprintf("%s: %.2f", auc_data$method[j], auc_data$AUC[j]),
        color = method_colors[auc_data$method[j]],
        hjust = 1,
        vjust = 0,
        fontface = "bold",
        size = 4
      )
  }

  plot_list[[i]] <- p
}

# Combine plots
roc_combined <- wrap_plots(plot_list, nrow = 1) +
  plot_annotation(
    caption = sprintf(
      "Black circles: p-value threshold %.0e | Red circles: p-value threshold 0.05",
      pthr
    ),
    theme = theme(
      plot.caption = element_text(hjust = 0.9, face = "italic", color = "gray")
    )
  )

print(roc_combined)

# Save ROC plots
ggsave(
  file.path(figure_dir, "ROC_curves_perturbations_methods.png"),
  roc_combined,
  width = 16,
  height = 4,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "ROC_curves_perturbations_methods.pdf"),
  roc_combined,
  width = 16,
  height = 4
)

cat("  ROC curves saved to figures/\n")

cat("\n")
cat("Step 3: Creating Precision-Recall curves...\n")

# ============================================================================
# Precision-Recall Curves
# ============================================================================

plot_list_pr <- list()

for (i in seq_along(data_names)) {
  pert <- unname(data_names)[i]

  cur_data <- results_data %>%
    filter(perturbation == pert) %>%
    mutate(
      Precision = ifelse(
        is.nan(Precision) | is.infinite(Precision),
        0,
        Precision
      ),
      Recall = ifelse(is.nan(Recall) | is.infinite(Recall), 0, Recall)
    ) %>%
    filter(!(Recall == 0 & Precision == 0)) %>%
    arrange(desc(Recall))

  cur_data <- cur_data %>%
    group_by(method) %>%
    group_modify(~ complete_curve_data(.x, "PR", "Recall", "Precision")) %>%
    ungroup()

  # Calculate AUPRC for each method
  auprc_data <- cur_data %>%
    group_by(method) %>%
    summarise(
      AUPRC = pracma::trapz(sort(Recall), Precision[order(Recall)]),
      .groups = "drop"
    )

  # Get point at threshold
  pthr_data <- cur_data %>%
    group_by(method) %>%
    filter(abs(threshold - pthr) == min(abs(threshold - pthr))) %>%
    ungroup()

  # Create plot
  p <- ggplot(cur_data, aes(x = Recall, y = Precision, color = method)) +
    geom_line(linewidth = 1) +
    geom_point(
      data = pthr_data,
      aes(x = Recall, y = Precision),
      size = 5,
      shape = 21,
      fill = NA,
      stroke = 1.5,
      color = "black"
    ) +
    scale_color_manual(values = method_colors) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05)) +
    labs(
      title = pert,
      x = "Recall",
      y = "Precision"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.position = "none",
      panel.grid.major = element_line(
        color = "lightgrey",
        linetype = "dashed",
        linewidth = 0.3
      ),
      panel.grid.minor = element_line(
        color = "lightgrey",
        linetype = "dashed",
        linewidth = 0.2
      )
    )

  # Add AUPRC annotations
  for (j in seq_len(nrow(auprc_data))) {
    p <- p +
      annotate(
        "text",
        x = 0.05,
        y = 0.05 + (j - 1) * 0.075,
        label = sprintf("%s: %.2f", auprc_data$method[j], auprc_data$AUPRC[j]),
        color = method_colors[auprc_data$method[j]],
        hjust = 0,
        vjust = 0,
        fontface = "bold",
        size = 4
      )
  }

  plot_list_pr[[i]] <- p
}

# Combine plots
pr_combined <- wrap_plots(plot_list_pr, nrow = 1) +
  plot_annotation(
    caption = sprintf(
      "Black circles represent Precision and Recall at p-value threshold of %.0e",
      pthr
    ),
    theme = theme(
      plot.caption = element_text(hjust = 0.9, face = "italic", color = "gray")
    )
  )

print(pr_combined)

# Save PR plots
ggsave(
  file.path(figure_dir, "PR_curves_perturbations_methods.png"),
  pr_combined,
  width = 16,
  height = 4,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "PR_curves_perturbations_methods.pdf"),
  pr_combined,
  width = 16,
  height = 4
)

cat("  PR curves saved to figures/\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n")
cat("Step 4: Generating summary statistics...\n")

# Calculate AUC and AUPRC for each method and perturbation
summary_stats <- results_data %>%
  group_by(perturbation, method) %>%
  summarise(
    AUC = pracma::trapz(sort(FPR), TPR[order(FPR)]),
    AUPRC = pracma::trapz(sort(Recall), Precision[order(Recall)]),
    .groups = "drop"
  ) %>%
  arrange(perturbation, desc(AUC))

print(summary_stats)

# Save summary
write.csv(
  summary_stats,
  file.path(result_dir, "performance_summary.csv"),
  row.names = FALSE
)

cat("\n")
cat("=", rep("=", 78), "\n", sep = "")
cat("              IDENTIFICATION BENCHMARK COMPLETED                   \n")
cat("=", rep("=", 78), "\n", sep = "")
cat("\n")
cat("Figures saved to:", figure_dir, "\n")
cat("  - ROC_curves_perturbations_methods.png/pdf\n")
cat("  - PR_curves_perturbations_methods.png/pdf\n")
cat("\n")

stop

# ============================================================================
# GROUPING BENCHMARK - Proteoform Detection (COPF and ProteoForge only)
# ============================================================================

cat("\n")
cat("=", rep("=", 78), "\n", sep = "")
cat("           GROUPING BENCHMARK - PROTEOFORM DETECTION              \n")
cat("=", rep("=", 78), "\n", sep = "")
cat("\n")

# Note: PeCorA is excluded as it only flags individual peptides, not groups
grouping_methods <- c("COPF", "ProteoForge", "DCF", "DCFBaseline")
grouping_data_ids <- c("2pep", "050pep", "random") # Exclude 1pep (trivial case)
grouping_data_names <- c(
  "2pep" = "2 Peptides",
  "050pep" = "%50 Peptides",
  "random" = "Random (2 to %50) Peptides"
)

# Helper function to update COPF proteoform grouping
update_copf_grouping <- function(data, pval_thr, score_thr = NULL) {
  # Determine significance
  if (!is.null(score_thr)) {
    data$isSignificant <- !is.na(data$proteoform_score_pval) &
      data$proteoform_score_pval <= pval_thr &
      !is.na(data$proteoform_score) &
      data$proteoform_score >= score_thr
  } else {
    data$isSignificant <- !is.na(data$proteoform_score_pval) &
      data$proteoform_score_pval <= pval_thr
  }

  # Group by protein to determine WithProteoform status
  # A protein has proteoforms if: any peptide is significant AND cluster > 1

  # Aggregate at protein level - efficient vectorized operations
  any_sig <- aggregate(
    isSignificant ~ protein_id,
    data = data,
    FUN = function(x) any(x, na.rm = TRUE)
  )

  n_clusters <- aggregate(
    cluster ~ protein_id,
    data = data,
    FUN = function(x) length(unique(x))
  )

  pert_protein <- aggregate(
    perturbed_protein ~ protein_id,
    data = data,
    FUN = function(x) x[1]
  )

  # Merge protein-level summaries
  protein_summary <- merge(any_sig, n_clusters, by = "protein_id")
  protein_summary <- merge(protein_summary, pert_protein, by = "protein_id")

  # WithProteoform = any significant peptide AND more than 1 cluster
  protein_summary$WithProteoform <- protein_summary$isSignificant &
    protein_summary$cluster > 1

  protein_summary[, c("protein_id", "perturbed_protein", "WithProteoform")]
}

cat("Step 1: Processing grouping data...\n")

# Process grouping results
grouping_results <- data.frame()

for (data_id in grouping_data_ids) {
  for (method in grouping_methods) {
    cat(sprintf("  Processing %s - %s (grouping)...\n", method, data_id))

    # Read data
    file_path <- file.path(
      result_dir,
      paste0(method, "_", data_id, "_result.feather")
    )
    data <- read_feather(file_path)

    # Calculate metrics for each threshold
    metric_data_list <- lapply(thresholds, function(thr) {
      if (method == "COPF") {
        # Update COPF grouping
        result <- update_copf_grouping(data, pval_thr = thr, score_thr = NULL)
        result <- result %>%
          select(protein_id, perturbed_protein, WithProteoform) %>%
          distinct()
      } else if (method == "ProteoForge") {
        # Process ProteoForge grouping
        result <- data %>%
          select(
            protein_id,
            peptide_id,
            perturbed_protein,
            perturbed_peptide,
            adj_pval,
            ClusterID
          ) %>%
          distinct() %>%
          mutate(
            PTM_id = ifelse(adj_pval < thr, ClusterID, -1),
            isPTM = PTM_id > 0
          ) %>%
          group_by(protein_id, PTM_id) %>%
          mutate(isGroup = n() > 1) %>%
          ungroup() %>%
          mutate(dPF = isGroup & isPTM) %>%
          group_by(protein_id) %>%
          summarise(
            perturbed_protein = dplyr::first(perturbed_protein),
            WithProteoform = any(dPF),
            .groups = "drop"
          )
      } else if (method %in% c("DCF", "DCFBaseline")) {
        # Process DCF / DCFBaseline grouping — identical column structure
        result <- data %>%
          select(
            protein_id,
            peptide_id,
            perturbed_protein,
            p_adj,
            dCF
          ) %>%
          distinct() %>%
          mutate(
            isDiscordant = !is.na(p_adj) & p_adj < thr,
            isAlternateForm = !is.na(dCF) & dCF != "dCF-1" & dCF != "dCF0"
          ) %>%
          group_by(protein_id) %>%
          summarise(
            perturbed_protein = dplyr::first(perturbed_protein),
            any_discordant = any(isDiscordant, na.rm = TRUE),
            has_alternate_form = any(
              isDiscordant & isAlternateForm,
              na.rm = TRUE
            ),
            .groups = "drop"
          ) %>%
          mutate(
            WithProteoform = any_discordant & has_alternate_form
          ) %>%
          select(protein_id, perturbed_protein, WithProteoform)
      }

      # Calculate confusion matrix
      true_labels <- result$perturbed_protein
      pred_labels <- result$WithProteoform

      TP <- sum(pred_labels & true_labels, na.rm = TRUE)
      FP <- sum(pred_labels & !true_labels, na.rm = TRUE)
      TN <- sum(!pred_labels & !true_labels, na.rm = TRUE)
      FN <- sum(!pred_labels & true_labels, na.rm = TRUE)

      # Calculate metrics
      TPR <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
      FPR <- ifelse(FP + TN > 0, FP / (FP + TN), 0)
      Precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
      Recall <- TPR

      data.frame(
        threshold = thr,
        TP = TP,
        FP = FP,
        TN = TN,
        FN = FN,
        TPR = TPR,
        FPR = FPR,
        Precision = Precision,
        Recall = Recall
      )
    })

    metric_data <- do.call(rbind, metric_data_list)
    metric_data$method <- method
    metric_data$perturbation <- data_id

    grouping_results <- rbind(grouping_results, metric_data)
  }
}

# Replace perturbation IDs with display names
grouping_results$perturbation <- grouping_data_names[
  grouping_results$perturbation
]
grouping_results$perturbation <- factor(
  grouping_results$perturbation,
  levels = unname(grouping_data_names)
)

cat(sprintf("\nGrouping data processed: %d rows\n", nrow(grouping_results)))

# Save grouping data
grouping_output_file <- file.path(
  result_dir,
  "peptide_grouping_performance_data.feather"
)
write_feather(grouping_results, grouping_output_file)
write.csv(
  grouping_results,
  sub(".feather", ".csv", grouping_output_file),
  row.names = FALSE
)

cat("\n")
cat("Step 2: Creating grouping ROC curves...\n")

# ============================================================================
# Grouping ROC Curves
# ============================================================================

grouping_colors <- c(
  "COPF" = color_palette[1],
  "ProteoForge" = color_palette[3],
  "DCF" = color_palette[4],
  "DCFBaseline" = color_palette[5]
)

plot_list_group_roc <- list()

for (i in seq_along(grouping_data_names)) {
  pert <- unname(grouping_data_names)[i]

  cur_data <- grouping_results %>%
    filter(perturbation == pert) %>%
    group_by(method) %>%
    group_modify(~ complete_curve_data(.x, "ROC", "FPR", "TPR")) %>%
    ungroup() %>%
    arrange(method, FPR)

  # Calculate AUC
  auc_data <- cur_data %>%
    group_by(method) %>%
    summarise(
      AUC = pracma::trapz(sort(FPR), TPR[order(FPR)]),
      .groups = "drop"
    )

  # Get point at threshold
  pthr_data <- cur_data %>%
    filter(abs(threshold - pthr) == min(abs(threshold - pthr)))

  pthr_005_data <- cur_data %>%
    group_by(method) %>%
    filter(abs(threshold - 0.05) == min(abs(threshold - 0.05))) %>%
    ungroup()

  # Create plot
  p <- ggplot(cur_data, aes(x = FPR, y = TPR, color = method)) +
    geom_line(linewidth = 1) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "black"
    ) +
    geom_point(
      data = pthr_data,
      aes(x = FPR, y = TPR),
      size = 5,
      shape = 21,
      fill = NA,
      stroke = 1.5,
      color = "black"
    ) +
    geom_point(
      data = pthr_005_data,
      aes(x = FPR, y = TPR),
      size = 3,
      shape = 21,
      fill = NA,
      stroke = 1.5,
      color = "red"
    ) +
    scale_color_manual(values = grouping_colors) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05)) +
    labs(
      title = pert,
      x = "FPR",
      y = "TPR"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.position = "none",
      panel.grid.major = element_line(
        color = "lightgrey",
        linetype = "dashed",
        linewidth = 0.3
      ),
      panel.grid.minor = element_line(
        color = "lightgrey",
        linetype = "dashed",
        linewidth = 0.2
      )
    )

  # Add AUC annotations
  for (j in seq_len(nrow(auc_data))) {
    p <- p +
      annotate(
        "text",
        x = 0.95,
        y = 0.05 + (j - 1) * 0.1,
        label = sprintf("%s: %.2f", auc_data$method[j], auc_data$AUC[j]),
        color = grouping_colors[auc_data$method[j]],
        hjust = 1,
        vjust = 0,
        fontface = "bold",
        size = 4
      )
  }

  plot_list_group_roc[[i]] <- p
}

# Combine plots
roc_grouping_combined <- wrap_plots(plot_list_group_roc, nrow = 1) +
  plot_annotation(
    caption = sprintf(
      "Grouping: Black circles: p-value threshold %.0e | Red circles: p-value threshold 0.05",
      pthr
    ),
    theme = theme(
      plot.caption = element_text(hjust = 0.9, face = "italic", color = "gray")
    )
  )

print(roc_grouping_combined)

# Save grouping ROC plots
ggsave(
  file.path(figure_dir, "ROC_curves_grouping_methods.png"),
  roc_grouping_combined,
  width = 12,
  height = 4,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "ROC_curves_grouping_methods.pdf"),
  roc_grouping_combined,
  width = 12,
  height = 4
)

cat("  Grouping ROC curves saved to figures/\n")

cat("\n")
cat("Step 3: Creating grouping Precision-Recall curves...\n")

# ============================================================================
# Grouping PR Curves
# ============================================================================

plot_list_group_pr <- list()

for (i in seq_along(grouping_data_names)) {
  pert <- unname(grouping_data_names)[i]

  cur_data <- grouping_results %>%
    filter(perturbation == pert) %>%
    mutate(
      Precision = ifelse(
        is.nan(Precision) | is.infinite(Precision),
        0,
        Precision
      ),
      Recall = ifelse(is.nan(Recall) | is.infinite(Recall), 0, Recall)
    ) %>%
    filter(!(Recall == 0 & Precision == 0)) %>%
    arrange(desc(Recall))

  cur_data <- cur_data %>%
    group_by(method) %>%
    group_modify(~ complete_curve_data(.x, "PR", "Recall", "Precision")) %>%
    ungroup()

  # Calculate AUPRC
  auprc_data <- cur_data %>%
    group_by(method) %>%
    summarise(
      AUPRC = pracma::trapz(sort(Recall), Precision[order(Recall)]),
      .groups = "drop"
    )

  # Get point at threshold
  pthr_data <- cur_data %>%
    group_by(method) %>%
    filter(abs(threshold - pthr) == min(abs(threshold - pthr))) %>%
    ungroup()

  # Create plot
  p <- ggplot(cur_data, aes(x = Recall, y = Precision, color = method)) +
    geom_line(linewidth = 1) +
    geom_point(
      data = pthr_data,
      aes(x = Recall, y = Precision),
      size = 5,
      shape = 21,
      fill = NA,
      stroke = 1.5,
      color = "black"
    ) +
    scale_color_manual(values = grouping_colors) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05)) +
    labs(
      title = pert,
      x = "Recall",
      y = "Precision"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.position = "none",
      panel.grid.major = element_line(
        color = "lightgrey",
        linetype = "dashed",
        linewidth = 0.3
      ),
      panel.grid.minor = element_line(
        color = "lightgrey",
        linetype = "dashed",
        linewidth = 0.2
      )
    )

  # Add AUPRC annotations
  for (j in seq_len(nrow(auprc_data))) {
    p <- p +
      annotate(
        "text",
        x = 0.05,
        y = 0.05 + (j - 1) * 0.1,
        label = sprintf("%s: %.2f", auprc_data$method[j], auprc_data$AUPRC[j]),
        color = grouping_colors[auprc_data$method[j]],
        hjust = 0,
        vjust = 0,
        fontface = "bold",
        size = 4
      )
  }

  plot_list_group_pr[[i]] <- p
}

# Combine plots
pr_grouping_combined <- wrap_plots(plot_list_group_pr, nrow = 1) +
  plot_annotation(
    caption = sprintf(
      "Grouping: Black circles at p-value threshold of %.0e",
      pthr
    ),
    theme = theme(
      plot.caption = element_text(hjust = 0.9, face = "italic", color = "gray")
    )
  )

print(pr_grouping_combined)

# Save grouping PR plots
ggsave(
  file.path(figure_dir, "PR_curves_grouping_methods.png"),
  pr_grouping_combined,
  width = 12,
  height = 4,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "PR_curves_grouping_methods.pdf"),
  pr_grouping_combined,
  width = 12,
  height = 4
)

cat("  Grouping PR curves saved to figures/\n")

# ============================================================================
# Grouping Summary Statistics
# ============================================================================

cat("\n")
cat("Step 4: Generating grouping summary statistics...\n")

# Calculate AUC and AUPRC for grouping
grouping_summary <- grouping_results %>%
  group_by(perturbation, method) %>%
  summarise(
    AUC = pracma::trapz(sort(FPR), TPR[order(FPR)]),
    AUPRC = pracma::trapz(sort(Recall), Precision[order(Recall)]),
    .groups = "drop"
  ) %>%
  arrange(perturbation, desc(AUC))

print(grouping_summary)

# Save grouping summary
write.csv(
  grouping_summary,
  file.path(result_dir, "grouping_performance_summary.csv"),
  row.names = FALSE
)

cat("\n")
cat("Step 5: Creating proteoform count comparison plots...\n")

# ============================================================================
# Proteoform Count Comparisons (DCF vs COPF and ProteoForge)
# ============================================================================

# Function to count proteoforms per protein for each method
count_proteoforms <- function(data, method, pval_thr = 0.001) {
  if (method == "COPF") {
    # COPF: count unique clusters
    # Exclude cluster 1 (canonical) and cluster 100 (singletons/unassigned)
    # This matches DCF's exclusion of dCF0 (canonical) and dCF-1 (singletons)
    result <- data %>%
      filter(
        !is.na(proteoform_score_pval) & proteoform_score_pval <= pval_thr
      ) %>%
      filter(cluster > 1 & cluster != 100) %>%
      group_by(protein_id) %>%
      summarise(
        n_proteoforms = length(unique(cluster)),
        .groups = "drop"
      )
  } else if (method == "ProteoForge") {
    # ProteoForge: count unique ClusterIDs
    # ClusterID 1 is the canonical/major cluster (like COPF cluster 1, DCF's dCF0)
    # Singletons/outliers are handled by the hybrid_outlier_cut clustering method
    # We exclude ClusterID 1 (canonical) to match the other methods
    result <- data %>%
      filter(!is.na(adj_pval) & adj_pval <= pval_thr) %>%
      filter(ClusterID > 1) %>%
      group_by(protein_id) %>%
      summarise(
        n_proteoforms = length(unique(ClusterID)),
        .groups = "drop"
      )
  } else if (method %in% c("DCF", "DCFBaseline")) {
    # DCF / DCFBaseline: count unique dCF labels (excluding dCF0 and dCF-1)
    result <- data %>%
      filter(!is.na(p_adj) & p_adj <= pval_thr) %>%
      filter(!is.na(dCF) & dCF != "dCF0" & dCF != "dCF-1") %>%
      group_by(protein_id) %>%
      summarise(
        n_proteoforms = length(unique(dCF)),
        .groups = "drop"
      )
  }
  result
}

# Create comparison plots for each dataset
comparison_plots <- list()

for (data_id in grouping_data_ids) {
  # Load data for all three methods
  copf_data <- read_feather(file.path(
    result_dir,
    paste0("COPF_", data_id, "_result.feather")
  ))
  pf_data <- read_feather(file.path(
    result_dir,
    paste0("ProteoForge_", data_id, "_result.feather")
  ))
  dcf_data <- read_feather(file.path(
    result_dir,
    paste0("DCF_", data_id, "_result.feather")
  ))
  dcfbaseline_data <- read_feather(file.path(
    result_dir,
    paste0("DCFBaseline_", data_id, "_result.feather")
  ))

  # Count unique proteoform labels per protein (excluding canonical/singletons)
  # COPF: exclude cluster 100 (singletons)
  copf_counts <- copf_data %>%
    filter(cluster != 100) %>%
    group_by(protein_id) %>%
    summarise(n_proteoforms = length(unique(cluster)), .groups = "drop")

  # ProteoForge: exclude ClusterID 1 (canonical)
  pf_counts <- pf_data %>%
    group_by(protein_id) %>%
    summarise(n_proteoforms = length(unique(ClusterID)), .groups = "drop")

  # DCF: exclude dCF-1 (singletons)
  dcf_counts <- dcf_data %>%
    filter(!is.na(dCF) & dCF != "dCF-1") %>%
    group_by(protein_id) %>%
    summarise(n_proteoforms = length(unique(dCF)), .groups = "drop")

  # DCFBaseline: same exclusion logic as DCF
  dcfbaseline_counts <- dcfbaseline_data %>%
    filter(!is.na(dCF) & dCF != "dCF-1") %>%
    group_by(protein_id) %>%
    summarise(n_proteoforms = length(unique(dCF)), .groups = "drop")

  # Merge for comparison (DCF vs COPF)
  copf_comparison <- merge(
    dcf_counts,
    copf_counts,
    by = "protein_id",
    all = TRUE,
    suffixes = c("_dcf", "_copf")
  )
  copf_comparison[is.na(copf_comparison)] <- 0

  # Merge for comparison (DCF vs ProteoForge)
  pf_comparison <- merge(
    dcf_counts,
    pf_counts,
    by = "protein_id",
    all = TRUE,
    suffixes = c("_dcf", "_pf")
  )
  pf_comparison[is.na(pf_comparison)] <- 0

  # Merge for comparison (COPF vs ProteoForge)
  copf_pf_comparison <- merge(
    copf_counts,
    pf_counts,
    by = "protein_id",
    all = TRUE,
    suffixes = c("_copf", "_pf")
  )
  copf_pf_comparison[is.na(copf_pf_comparison)] <- 0

  # Create count matrix for COPF heatmap
  copf_count_table <- table(
    copf_comparison$n_proteoforms_dcf,
    copf_comparison$n_proteoforms_copf
  )
  copf_count_df <- as.data.frame(copf_count_table)
  names(copf_count_df) <- c("DCF", "COPF", "count")
  copf_count_df$DCF <- as.numeric(as.character(copf_count_df$DCF))
  copf_count_df$COPF <- as.numeric(as.character(copf_count_df$COPF))
  copf_count_df <- copf_count_df[copf_count_df$count > 0, ]

  # Create COPF comparison heatmap
  p_copf <- ggplot(
    copf_count_df,
    aes(x = DCF, y = COPF, fill = count)
  ) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = count),
      color = "white",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_gradient(
      low = "#e0f3f8",
      high = "#00A087",
      trans = "log1p",
      breaks = c(1, 5, 10, 50, 100),
      name = "Count"
    ) +
    scale_x_continuous(
      limits = c(0, 10),
      breaks = 0:10,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 10),
      breaks = 0:10,
      expand = c(0, 0)
    ) +
    labs(
      title = paste(grouping_data_names[data_id], "- DCF vs COPF"),
      x = "Number of proteoforms (DCF)",
      y = "Number of proteoforms (COPF)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.title = element_text(size = 10),
      panel.grid.major = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed(ratio = 1)

  # Create count matrix for ProteoForge heatmap
  pf_count_table <- table(
    pf_comparison$n_proteoforms_dcf,
    pf_comparison$n_proteoforms_pf
  )
  pf_count_df <- as.data.frame(pf_count_table)
  names(pf_count_df) <- c("DCF", "ProteoForge", "count")
  pf_count_df$DCF <- as.numeric(as.character(pf_count_df$DCF))
  pf_count_df$ProteoForge <- as.numeric(as.character(pf_count_df$ProteoForge))
  pf_count_df <- pf_count_df[pf_count_df$count > 0, ]

  # Create ProteoForge comparison heatmap
  p_pf <- ggplot(
    pf_count_df,
    aes(x = DCF, y = ProteoForge, fill = count)
  ) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = count),
      color = "white",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_gradient(
      low = "#e0f3f8",
      high = "#00A087",
      trans = "log1p",
      breaks = c(1, 5, 10, 50, 100),
      name = "Count"
    ) +
    scale_x_continuous(
      limits = c(0, 10),
      breaks = 0:10,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 10),
      breaks = 0:10,
      expand = c(0, 0)
    ) +
    labs(
      title = paste(grouping_data_names[data_id], "- DCF vs ProteoForge"),
      x = "Number of proteoforms (DCF)",
      y = "Number of proteoforms (ProteoForge)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.title = element_text(size = 10),
      panel.grid.major = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed(ratio = 1)

  # Create count matrix for COPF vs ProteoForge heatmap
  copf_pf_count_table <- table(
    copf_pf_comparison$n_proteoforms_copf,
    copf_pf_comparison$n_proteoforms_pf
  )
  copf_pf_count_df <- as.data.frame(copf_pf_count_table)
  names(copf_pf_count_df) <- c("COPF", "ProteoForge", "count")
  copf_pf_count_df$COPF <- as.numeric(as.character(copf_pf_count_df$COPF))
  copf_pf_count_df$ProteoForge <- as.numeric(as.character(
    copf_pf_count_df$ProteoForge
  ))
  copf_pf_count_df <- copf_pf_count_df[copf_pf_count_df$count > 0, ]

  # Create COPF vs ProteoForge comparison heatmap
  p_copf_pf <- ggplot(
    copf_pf_count_df,
    aes(x = COPF, y = ProteoForge, fill = count)
  ) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = count),
      color = "white",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_gradient(
      low = "#e0f3f8",
      high = "#00A087",
      trans = "log1p",
      breaks = c(1, 5, 10, 50, 100),
      name = "Count"
    ) +
    scale_x_continuous(
      limits = c(0, 10),
      breaks = 0:10,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 10),
      breaks = 0:10,
      expand = c(0, 0)
    ) +
    labs(
      title = paste(grouping_data_names[data_id], "- COPF vs ProteoForge"),
      x = "Number of proteoforms (COPF)",
      y = "Number of proteoforms (ProteoForge)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.title = element_text(size = 10),
      panel.grid.major = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed(ratio = 1)

  # Merge for comparison (DCFBaseline vs DCF)
  dcfbaseline_dcf_comparison <- merge(
    dcfbaseline_counts,
    dcf_counts,
    by = "protein_id",
    all = TRUE,
    suffixes = c("_dcfbaseline", "_dcf")
  )
  dcfbaseline_dcf_comparison[is.na(dcfbaseline_dcf_comparison)] <- 0

  # Create count matrix for DCFBaseline vs DCF heatmap
  dcfbaseline_dcf_count_table <- table(
    dcfbaseline_dcf_comparison$n_proteoforms_dcfbaseline,
    dcfbaseline_dcf_comparison$n_proteoforms_dcf
  )
  dcfbaseline_dcf_count_df <- as.data.frame(dcfbaseline_dcf_count_table)
  names(dcfbaseline_dcf_count_df) <- c("DCFBaseline", "DCF", "count")
  dcfbaseline_dcf_count_df$DCFBaseline <- as.numeric(
    as.character(dcfbaseline_dcf_count_df$DCFBaseline)
  )
  dcfbaseline_dcf_count_df$DCF <- as.numeric(
    as.character(dcfbaseline_dcf_count_df$DCF)
  )
  dcfbaseline_dcf_count_df <- dcfbaseline_dcf_count_df[
    dcfbaseline_dcf_count_df$count > 0,
  ]

  p_dcfbaseline <- ggplot(
    dcfbaseline_dcf_count_df,
    aes(x = DCFBaseline, y = DCF, fill = count)
  ) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = count),
      color = "white",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_gradient(
      low = "#e0f3f8",
      high = "#00A087",
      trans = "log1p",
      breaks = c(1, 5, 10, 50, 100),
      name = "Count"
    ) +
    scale_x_continuous(
      limits = c(0, 10),
      breaks = 0:10,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 10),
      breaks = 0:10,
      expand = c(0, 0)
    ) +
    labs(
      title = paste(grouping_data_names[data_id], "- DCFBaseline vs DCF"),
      x = "Number of proteoforms (DCFBaseline)",
      y = "Number of proteoforms (DCF)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.title = element_text(size = 10),
      panel.grid.major = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed(ratio = 1)

  comparison_plots[[paste0("copf_", data_id)]] <- p_copf
  comparison_plots[[paste0("pf_", data_id)]] <- p_pf
  comparison_plots[[paste0("copf_pf_", data_id)]] <- p_copf_pf
  comparison_plots[[paste0("dcfbaseline_", data_id)]] <- p_dcfbaseline
}

# Organize plots for 4x3 grid layout
# Row 1: DCF vs COPF (2pep, 050pep, random)
# Row 2: DCF vs ProteoForge (2pep, 050pep, random)
# Row 3: COPF vs ProteoForge (2pep, 050pep, random)
# Row 4: DCFBaseline vs DCF (2pep, 050pep, random)
top_row <- plot_grid(
  comparison_plots[["copf_2pep"]],
  comparison_plots[["copf_050pep"]],
  comparison_plots[["copf_random"]],
  ncol = 3,
  align = "hv",
  axis = "tblr"
)

middle_row <- plot_grid(
  comparison_plots[["pf_2pep"]],
  comparison_plots[["pf_050pep"]],
  comparison_plots[["pf_random"]],
  ncol = 3,
  align = "hv",
  axis = "tblr"
)

bottom_row <- plot_grid(
  comparison_plots[["copf_pf_2pep"]],
  comparison_plots[["copf_pf_050pep"]],
  comparison_plots[["copf_pf_random"]],
  ncol = 3,
  align = "hv",
  axis = "tblr"
)

dcfbaseline_row <- plot_grid(
  comparison_plots[["dcfbaseline_2pep"]],
  comparison_plots[["dcfbaseline_050pep"]],
  comparison_plots[["dcfbaseline_random"]],
  ncol = 3,
  align = "hv",
  axis = "tblr"
)

# Combine rows with title and caption
title <- ggdraw() +
  draw_label(
    "Proteoform Count Comparisons Across Methods",
    fontface = "bold",
    size = 14,
    hjust = 0.5
  )

comparison_combined <- plot_grid(
  title,
  top_row,
  middle_row,
  bottom_row,
  dcfbaseline_row,
  ncol = 1,
  rel_heights = c(0.06, 1, 1, 1, 1)
)

print(comparison_combined)

# Save comparison plots
ggsave(
  file.path(figure_dir, "proteoform_count_comparisons.png"),
  comparison_combined,
  width = 12,
  height = 12,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "proteoform_count_comparisons.pdf"),
  comparison_combined,
  width = 12,
  height = 12
)

cat("  Proteoform count comparison plots saved to figures/\n")

cat("\n")
cat("=", rep("=", 78), "\n", sep = "")
cat("                 GROUPING BENCHMARK COMPLETED                      \n")
cat("=", rep("=", 78), "\n", sep = "")
cat("\n")
cat("Grouping figures saved to:", figure_dir, "\n")
cat("  - ROC_curves_grouping_methods.png/pdf\n")
cat("  - PR_curves_grouping_methods.png/pdf\n")
cat("  - proteoform_count_comparisons.png/pdf\n")
cat("\n")

# Calculate total time
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "secs")

cat("Total execution time:", round(total_time, 2), "seconds\n")
cat("\n")
