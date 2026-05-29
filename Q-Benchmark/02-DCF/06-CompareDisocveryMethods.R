# Libraries
suppressMessages(library(arrow, warn.conflicts = FALSE))
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Change working directory to Q-Benchmark/02-DCF
setwd(here::here("Q-Benchmark/02-DCF"))

# Define output directory for comparison results
output_dir <- "./00-Data/04-PerformanceMetrics"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

################################################################################
# Definitions
################################################################################
# Definition of ground truth discordant peptides in ProteoMaker datasets
#  - A peptide is considered discordant if, after splitting Proteoform_ID by "|",
#    it does not exclusively carry the proteoform ID that occurs in the most
#    peptidoforms (the dominant proteoform). Concordant peptidoforms are those
#    whose every split token equals the dominant ID.

################################################################################
# Generate ground truth table from ProteoMaker datasets
################################################################################
# Datasets
datasets <- list(
  "Low" = "./00-Data/01-Unprepared/ProteoMaker_LowNA.csv",
  "Med" = "./00-Data/01-Unprepared/ProteoMaker_MedNA.csv",
  "High" = "./00-Data/01-Unprepared/ProteoMaker_HighNA.csv"
)

# Create ground truth tables
ground_truth_list <- list()
for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  data <- fread(file_path)

  # Subset to relevant columns
  data_subset <- data[, .(Accession, Peptidoform, Proteoform_ID)]

  # Find the dominant proteoform ID per protein (Accession)
  pf_ids <- data_subset[,
    .(id = unique(unlist(strsplit(Proteoform_ID, "\\|")))),
    by = .(Accession, Peptidoform)
  ]
  dominant_ids <- pf_ids[, .N, by = .(Accession, id)][
    order(Accession, -N)
  ][, .SD[1L], by = Accession][, .(Accession, dominant_id = id)]

  data_subset <- merge(data_subset, dominant_ids, by = "Accession")

  # A peptidoform is discordant if not all of its split Proteoform_IDs equal
  # the dominant proteoform for its protein
  data_subset[,
    isDiscordant := !all(
      unique(unlist(strsplit(Proteoform_ID, "\\|"))) == dominant_id
    ),
    by = .(Accession, Peptidoform)
  ]
  data_subset[, dominant_id := NULL]

  # Store in list
  ground_truth_list[[name]] <- data_subset
}


################################################################################
# Calculate performance metrics for each method for each dataset
# - Metrics: FPR, TPR, Precision, Recall
# - Metrics calculated at different thresholds from 1 to 1e-15
################################################################################

# Global variables
# fmt: skip
thresholds <- seq(0.975, 0.01, by = -0.025)
thresholds <- c(
  thresholds,
  0.009,
  0.008,
  0.007,
  0.006,
  0.005,
  0.004,
  0.003,
  0.002,
  0.001,
  1e-3,
  1e-4,
  1e-5,
  1e-6,
  1e-7,
  1e-8,
  1e-9,
  1e-10,
  1e-11,
  1e-12,
  1e-13,
  1e-14,
  1e-15
)
thresholds <- sort(thresholds, decreasing = FALSE)

########################
## Metrics for PeCorA ##
########################
# Note: All peptides not returned by PeCorA are added back with adj_pval = 1

# Datasets
datasets <- list(
  "LowNA" = "./00-Data/03-Results/PeCorA/LowNA_imputed_PeCorA_results.csv",
  "LowNoNA" = "./00-Data/03-Results/PeCorA/LowNoNA_PeCorA_results.csv",
  "MedNA" = "./00-Data/03-Results/PeCorA/MedNA_imputed_PeCorA_results.csv",
  "MedNoNA" = "./00-Data/03-Results/PeCorA/MedNoNA_PeCorA_results.csv",
  "HighNA" = "./00-Data/03-Results/PeCorA/HighNA_imputed_PeCorA_results.csv",
  "HighNoNA" = "./00-Data/03-Results/PeCorA/HighNoNA_PeCorA_results.csv"
)

PeCorA_vs_groundtruth <- list()

# Combine results with ground truth for each dataset
for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  # Read in results
  results <- fread(file_path)

  # Remove "_all" suffix from Peptide column
  results[, peptide := gsub("_all$", "", peptide)]

  # Get corresponding ground truth table
  ground_truth_name <- gsub("(NA|NoNA)$", "", name)
  ground_truth <- ground_truth_list[[ground_truth_name]]

  # Merge with ground truth by Peptidoform
  merged <- merge(
    ground_truth,
    results[, .(peptide, adj_pval)],
    by.x = "Peptidoform",
    by.y = "peptide",
    all = TRUE
  )

  # Set adj_pval to 1 for NA values (i.e., peptides not returned by PeCorA)
  merged[is.na(adj_pval), adj_pval := 1]

  PeCorA_vs_groundtruth[[name]] <- merged
}

# Initialize data table to store metrics for each dataset at different thresholds
PeCorA_metrics <- data.table(
  Dataset = character(),
  Threshold = numeric(),
  TP = integer(),
  FP = integer(),
  TN = integer(),
  FN = integer(),
  FPR = numeric(),
  TPR = numeric(),
  FDR = numeric(),
  Precision = numeric(),
  Recall = numeric()
)

# Calculate metrics at different adj_pval thresholds
data_metrics <- data.table(
  thresholds = thresholds,
  TP = integer(length(thresholds)),
  FP = integer(length(thresholds)),
  TN = integer(length(thresholds)),
  FN = integer(length(thresholds)),
  FPR = numeric(length(thresholds)),
  TPR = numeric(length(thresholds)),
  FDR = numeric(length(thresholds)),
  Precision = numeric(length(thresholds)),
  Recall = numeric(length(thresholds))
)

for (name in names(PeCorA_vs_groundtruth)) {
  dataset_metrics <- copy(data_metrics)

  merged <- PeCorA_vs_groundtruth[[name]]

  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]

    # Classify peptides as discordant or not based on threshold
    merged[, predicted_discordant := adj_pval < threshold]

    # Calculate TP, FP, TN, FN
    tp_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == TRUE])
    fp_val <- nrow(merged[isDiscordant == FALSE & predicted_discordant == TRUE])
    tn_val <- nrow(merged[
      isDiscordant == FALSE & predicted_discordant == FALSE
    ])
    fn_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == FALSE])

    # Calculate metrics
    fpr_val <- ifelse((fp_val + tn_val) > 0, fp_val / (fp_val + tn_val), NA)
    tpr_val <- ifelse((tp_val + fn_val) > 0, tp_val / (tp_val + fn_val), NA)
    fdr_val <- ifelse((tp_val + fp_val) > 0, fp_val / (tp_val + fp_val), NA)
    precision_val <- ifelse(
      (tp_val + fp_val) > 0,
      tp_val / (tp_val + fp_val),
      NA
    )
    recall_val <- tpr_val

    # Store metrics in data table
    dataset_metrics[
      i,
      `:=`(
        Threshold = threshold,
        TP = tp_val,
        FP = fp_val,
        TN = tn_val,
        FN = fn_val,
        FPR = fpr_val,
        TPR = tpr_val,
        FDR = fdr_val,
        Precision = precision_val,
        Recall = recall_val
      )
    ]
  }

  # Add dataset name to metrics
  dataset_metrics[, Dataset := name]

  # Append to overall metrics table
  PeCorA_metrics <- rbind(PeCorA_metrics, dataset_metrics, fill = TRUE)
}

# Save metrics to CSV
fwrite(
  PeCorA_metrics,
  file.path(output_dir, "PeCorA_Metrics.csv")
)

#############################
## Metrics for ProteoForge ##
#############################

# Datasets
datasets <- list(
  "LowNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNA_result.feather",
  "LowNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNoNA_result.feather",
  "MedNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNA_result.feather",
  "MedNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNoNA_result.feather",
  "HighNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNA_result.feather",
  "HighNoNA" = "./00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNoNA_result.feather"
)

ProteoForge_vs_groundtruth <- list()

# Combine results with ground truth for each dataset
for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  # Read in results
  results <- read_feather(file_path)
  results <- as.data.table(results)

  # Get corresponding ground truth table
  ground_truth_name <- gsub("(NA|NoNA)$", "", name)
  ground_truth <- ground_truth_list[[ground_truth_name]]

  # Convert to long format with feature_id as first row, Sample as column names, and adj_pval as values
  results_wide <- dcast(
    results,
    feature_id ~ Sample,
    value.var = "adj_pval"
  )

  # Remove all but one sample column and rename it to adj_pval
  results_wide <- results_wide[, .(
    feature_id,
    adj_pval = get(names(results_wide)[2])
  )]

  # Remove everything before and including "||"
  results_wide[, feature_id := sub(".*\\|\\|", "", feature_id)]

  # Combine with ground truth by Peptidoform
  merged <- merge(
    ground_truth,
    results_wide[, .(feature_id, adj_pval)],
    by.x = "Peptidoform",
    by.y = "feature_id",
    all = TRUE
  )

  # Set adj_pval to 1 for peptides not returned by ProteoForge
  merged[is.na(adj_pval), adj_pval := 1]

  ProteoForge_vs_groundtruth[[name]] <- merged
}

ProteForge_metrics <- data.table(
  Dataset = character(),
  Threshold = numeric(),
  TP = integer(),
  FP = integer(),
  TN = integer(),
  FN = integer(),
  FPR = numeric(),
  TPR = numeric(),
  FDR = numeric(),
  Precision = numeric(),
  Recall = numeric()
)

data_metrics <- data.table(
  thresholds = thresholds,
  TP = integer(length(thresholds)),
  FP = integer(length(thresholds)),
  TN = integer(length(thresholds)),
  FN = integer(length(thresholds)),
  FPR = numeric(length(thresholds)),
  TPR = numeric(length(thresholds)),
  FDR = numeric(length(thresholds)),
  Precision = numeric(length(thresholds)),
  Recall = numeric(length(thresholds))
)

for (name in names(ProteoForge_vs_groundtruth)) {
  dataset_metrics <- copy(data_metrics)

  merged <- ProteoForge_vs_groundtruth[[name]]

  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]

    # Classify peptides as discordant or not based on threshold
    merged[, predicted_discordant := adj_pval < threshold]

    # Calculate TP, FP, TN, FN
    tp_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == TRUE])
    fp_val <- nrow(merged[isDiscordant == FALSE & predicted_discordant == TRUE])
    tn_val <- nrow(merged[
      isDiscordant == FALSE & predicted_discordant == FALSE
    ])
    fn_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == FALSE])

    # Calculate metrics
    fpr_val <- ifelse((fp_val + tn_val) > 0, fp_val / (fp_val + tn_val), NA)
    tpr_val <- ifelse((tp_val + fn_val) > 0, tp_val / (tp_val + fn_val), NA)
    fdr_val <- ifelse((tp_val + fp_val) > 0, fp_val / (tp_val + fp_val), NA)
    precision_val <- ifelse(
      (tp_val + fp_val) > 0,
      tp_val / (tp_val + fp_val),
      NA
    )
    recall_val <- tpr_val

    # Store metrics in data table
    dataset_metrics[
      i,
      `:=`(
        Threshold = threshold,
        TP = tp_val,
        FP = fp_val,
        TN = tn_val,
        FN = fn_val,
        FPR = fpr_val,
        TPR = tpr_val,
        FDR = fdr_val,
        Precision = precision_val,
        Recall = recall_val
      )
    ]
  }

  # Add dataset name to metrics
  dataset_metrics[, Dataset := name]

  # Append to overall metrics table
  ProteForge_metrics <- rbind(ProteForge_metrics, dataset_metrics, fill = TRUE)
}

# Save metrics to CSV
fwrite(
  ProteForge_metrics,
  file.path(output_dir, "ProteoForge_Metrics.csv")
)

# #####################
# ## Metrics for RPC ##
# #####################
# # Note: peptides where p_adj == NA have p-value set to 1

# Datasets
# datasets <- list(
#   "LowNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerLowNA_results.csv",
#   "LowNoNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerLowNoNA_results.csv",
#   "MedNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerMedNA_results.csv",
#   "MedNoNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerMedNoNA_results.csv",
#   "HighNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerHighNA_results.csv",
#   "HighNoNA" = "./00-Data/03-Results/RPC/RPC_ProteoMakerHighNoNA_results.csv"
# )
#
# RPC_vs_groundtruth <- list()
#
# # Combine results with ground truth for each dataset
# for (name in names(datasets)) {
#   file_path <- datasets[[name]]
#   if (!file.exists(file_path)) {
#     warning(paste("File not found:", file_path))
#     next
#   }
#
#   # Read in results
#   results <- fread(file_path)
#
#   # Subset to relevant columns and rename p_adj to adj_pval
#   results_wide <- results[, .(
#     Peptide_ID,
#     adj_pval = p_adj
#   )]
#
#   # Remove everything before and including the last "_" in Peptide_ID to get Peptidoform
#   results_wide[,
#     Peptidoform := sub(".*_(?=[^_]+$)", "", Peptide_ID, perl = TRUE)
#   ]
#   results_wide[, Peptide_ID := NULL]
#
#   # Get corresponding ground truth table
#   ground_truth_name <- gsub("(NA|NoNA)$", "", name)
#   ground_truth <- ground_truth_list[[ground_truth_name]]
#
#   # Combine with ground truth by Peptidoform
#   merged <- merge(
#     ground_truth,
#     results_wide[, .(Peptidoform, adj_pval)],
#     by = "Peptidoform",
#     all = TRUE
#   )
#
#   # Set adj_pval to 1 for NA values (i.e., peptides not returned by RPC)
#   merged[is.na(adj_pval), adj_pval := 1]
#
#   # Store in list
#   RPC_vs_groundtruth[[name]] <- merged
# }
#
# RPC_metrics <- data.table(
#   Dataset = character(),
#   Threshold = numeric(),
#   TP = integer(),
#   FP = integer(),
#   TN = integer(),
#   FN = integer(),
#   FPR = numeric(),
#   TPR = numeric(),
#   FDR = numeric(),
#   Precision = numeric(),
#   Recall = numeric()
# )
#
# data_metrics <- data.table(
#   thresholds = thresholds,
#   TP = integer(length(thresholds)),
#   FP = integer(length(thresholds)),
#   TN = integer(length(thresholds)),
#   FN = integer(length(thresholds)),
#   FPR = numeric(length(thresholds)),
#   TPR = numeric(length(thresholds)),
#   FDR = numeric(length(thresholds)),
#   Precision = numeric(length(thresholds)),
#   Recall = numeric(length(thresholds))
# )
#
# for (name in names(RPC_vs_groundtruth)) {
#   dataset_metrics <- copy(data_metrics)
#
#   merged <- RPC_vs_groundtruth[[name]]
#
#   for (i in seq_along(thresholds)) {
#     threshold <- thresholds[i]
#
#     # Classify peptides as discordant or not based on threshold
#     merged[, predicted_discordant := adj_pval < threshold]
#
#     # Calculate TP, FP, TN, FN
#     tp_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == TRUE])
#     fp_val <- nrow(merged[isDiscordant == FALSE & predicted_discordant == TRUE])
#     tn_val <- nrow(merged[
#       isDiscordant == FALSE & predicted_discordant == FALSE
#     ])
#     fn_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == FALSE])
#
#     # Calculate metrics
#     fpr_val <- ifelse((fp_val + tn_val) > 0, fp_val / (fp_val + tn_val), NA)
#     tpr_val <- ifelse((tp_val + fn_val) > 0, tp_val / (tp_val + fn_val), NA)
#     fdr_val <- ifelse((tp_val + fp_val) > 0, fp_val / (tp_val + fp_val), NA)
#     precision_val <- ifelse(
#       (tp_val + fp_val) > 0,
#       tp_val / (tp_val + fp_val),
#       NA
#     )
#     recall_val <- tpr_val
#
#     # Store metrics in data table
#     dataset_metrics[
#       i,
#       `:=`(
#         Threshold = threshold,
#         TP = tp_val,
#         FP = fp_val,
#         TN = tn_val,
#         FN = fn_val,
#         FPR = fpr_val,
#         TPR = tpr_val,
#         FDR = fdr_val,
#         Precision = precision_val,
#         Recall = recall_val
#       )
#     ]
#   }
#
#   # Add dataset name to metrics
#   dataset_metrics[, Dataset := name]
#
#   # Append to overall metrics table
#   RPC_metrics <- rbind(RPC_metrics, dataset_metrics, fill = TRUE)
# }
#
# # Save metrics to CSV
# fwrite(
#   RPC_metrics,
#   file.path(output_dir, "RPC_Metrics.csv")
# )

#####################
## Metrics for DCF ##
#####################
# Note: peptides where p_adj == NA have p-value set to 1

# Datasets
datasets <- list(
  "LowNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerLowNA_results.csv",
  "LowNoNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerLowNoNA_results.csv",
  "MedNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerMedNA_results.csv",
  "MedNoNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerMedNoNA_results.csv",
  "HighNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerHighNA_results.csv",
  "HighNoNA" = "./00-Data/03-Results/DCF/DCF_ProteoMakerHighNoNA_results.csv"
)

DCF_vs_groundtruth <- list()

# Combine results with ground truth for each dataset
for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  # Read in results
  results <- fread(file_path)

  # Subset to relevant columns and rename p_adj to adj_pval
  results_wide <- results[, .(
    Peptide_ID,
    adj_pval = p_adj
  )]

  # Remove everything before and including the last "_" in Peptide_ID to get Peptidoform
  results_wide[,
    Peptidoform := sub(".*_(?=[^_]+$)", "", Peptide_ID, perl = TRUE)
  ]
  results_wide[, Peptide_ID := NULL]

  # Get corresponding ground truth table
  ground_truth_name <- gsub("(NA|NoNA)$", "", name)
  ground_truth <- ground_truth_list[[ground_truth_name]]

  # Combine with ground truth by Peptidoform
  merged <- merge(
    ground_truth,
    results_wide[, .(Peptidoform, adj_pval)],
    by = "Peptidoform",
    all = TRUE
  )

  # Set adj_pval to 1 for NA values (i.e., peptides not returned by DCF)
  merged[is.na(adj_pval), adj_pval := 1]

  # Store in list
  DCF_vs_groundtruth[[name]] <- merged
}

DCF_metrics <- data.table(
  Dataset = character(),
  Threshold = numeric(),
  TP = integer(),
  FP = integer(),
  TN = integer(),
  FN = integer(),
  FPR = numeric(),
  TPR = numeric(),
  FDR = numeric(),
  Precision = numeric(),
  Recall = numeric()
)

data_metrics <- data.table(
  thresholds = thresholds,
  TP = integer(length(thresholds)),
  FP = integer(length(thresholds)),
  TN = integer(length(thresholds)),
  FN = integer(length(thresholds)),
  FPR = numeric(length(thresholds)),
  TPR = numeric(length(thresholds)),
  FDR = numeric(length(thresholds)),
  Precision = numeric(length(thresholds)),
  Recall = numeric(length(thresholds))
)

for (name in names(DCF_vs_groundtruth)) {
  dataset_metrics <- copy(data_metrics)

  merged <- DCF_vs_groundtruth[[name]]

  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]

    # Classify peptides as discordant or not based on threshold
    merged[, predicted_discordant := adj_pval < threshold]

    # Calculate TP, FP, TN, FN
    tp_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == TRUE])
    fp_val <- nrow(merged[isDiscordant == FALSE & predicted_discordant == TRUE])
    tn_val <- nrow(merged[
      isDiscordant == FALSE & predicted_discordant == FALSE
    ])
    fn_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == FALSE])

    # Calculate metrics
    fpr_val <- ifelse((fp_val + tn_val) > 0, fp_val / (fp_val + tn_val), NA)
    tpr_val <- ifelse((tp_val + fn_val) > 0, tp_val / (tp_val + fn_val), NA)
    fdr_val <- ifelse((tp_val + fp_val) > 0, fp_val / (tp_val + fp_val), NA)
    precision_val <- ifelse(
      (tp_val + fp_val) > 0,
      tp_val / (tp_val + fp_val),
      NA
    )
    recall_val <- tpr_val

    # Store metrics in data table
    dataset_metrics[
      i,
      `:=`(
        Threshold = threshold,
        TP = tp_val,
        FP = fp_val,
        TN = tn_val,
        FN = fn_val,
        FPR = fpr_val,
        TPR = tpr_val,
        FDR = fdr_val,
        Precision = precision_val,
        Recall = recall_val
      )
    ]
  }

  # Add dataset name to metrics
  dataset_metrics[, Dataset := name]

  # Append to overall metrics table
  DCF_metrics <- rbind(DCF_metrics, dataset_metrics, fill = TRUE)
}

# Save metrics to CSV
fwrite(
  DCF_metrics,
  file.path(output_dir, "DCF_Metrics.csv")
)

##############################
## Metrics for DCFBaseline  ##
##############################
# Note: same output structure as DCF; missing peptides get adj_pval = 1

datasets <- list(
  "LowNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNA_results.csv",
  "LowNoNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNoNA_results.csv",
  "MedNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNA_results.csv",
  "MedNoNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNoNA_results.csv",
  "HighNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNA_results.csv",
  "HighNoNA" = "./00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNoNA_results.csv"
)

DCFBaseline_vs_groundtruth <- list()

for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  results <- fread(file_path)
  results_wide <- results[, .(Peptide_ID, adj_pval = p_adj)]
  results_wide[,
    Peptidoform := sub(".*_(?=[^_]+$)", "", Peptide_ID, perl = TRUE)
  ]
  results_wide[, Peptide_ID := NULL]
  ground_truth_name <- gsub("(NA|NoNA)$", "", name)
  ground_truth <- ground_truth_list[[ground_truth_name]]
  merged <- merge(
    ground_truth,
    results_wide[, .(Peptidoform, adj_pval)],
    by = "Peptidoform",
    all = TRUE
  )
  merged[is.na(adj_pval), adj_pval := 1]
  DCFBaseline_vs_groundtruth[[name]] <- merged
}

DCFBaseline_metrics <- data.table(
  Dataset = character(),
  Threshold = numeric(),
  TP = integer(),
  FP = integer(),
  TN = integer(),
  FN = integer(),
  FPR = numeric(),
  TPR = numeric(),
  FDR = numeric(),
  Precision = numeric(),
  Recall = numeric()
)

data_metrics <- data.table(
  thresholds = thresholds,
  TP = integer(length(thresholds)),
  FP = integer(length(thresholds)),
  TN = integer(length(thresholds)),
  FN = integer(length(thresholds)),
  FPR = numeric(length(thresholds)),
  TPR = numeric(length(thresholds)),
  FDR = numeric(length(thresholds)),
  Precision = numeric(length(thresholds)),
  Recall = numeric(length(thresholds))
)

for (name in names(DCFBaseline_vs_groundtruth)) {
  dataset_metrics <- copy(data_metrics)
  merged <- DCFBaseline_vs_groundtruth[[name]]
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    merged[, predicted_discordant := adj_pval < threshold]
    tp_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == TRUE])
    fp_val <- nrow(merged[isDiscordant == FALSE & predicted_discordant == TRUE])
    tn_val <- nrow(merged[
      isDiscordant == FALSE & predicted_discordant == FALSE
    ])
    fn_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == FALSE])
    fpr_val <- ifelse((fp_val + tn_val) > 0, fp_val / (fp_val + tn_val), NA)
    tpr_val <- ifelse((tp_val + fn_val) > 0, tp_val / (tp_val + fn_val), NA)
    fdr_val <- ifelse((tp_val + fp_val) > 0, fp_val / (tp_val + fp_val), NA)
    precision_val <- ifelse(
      (tp_val + fp_val) > 0,
      tp_val / (tp_val + fp_val),
      NA
    )
    recall_val <- tpr_val
    dataset_metrics[
      i,
      `:=`(
        Threshold = threshold,
        TP = tp_val,
        FP = fp_val,
        TN = tn_val,
        FN = fn_val,
        FPR = fpr_val,
        TPR = tpr_val,
        FDR = fdr_val,
        Precision = precision_val,
        Recall = recall_val
      )
    ]
  }
  dataset_metrics[, Dataset := name]
  DCFBaseline_metrics <- rbind(
    DCFBaseline_metrics,
    dataset_metrics,
    fill = TRUE
  )
}

fwrite(DCFBaseline_metrics, file.path(output_dir, "DCFBaseline_Metrics.csv"))

######################
## Metrics for COPF ##
######################
# Note: peptide_adj_pval = adj_pval for proteoform-2 peptides, 1.0 for canonical.
#       Missing peptides (not returned by COPF) get adj_pval = 1.

datasets <- list(
  "LowNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerLowNA_imputed_results.csv",
  "LowNoNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerLowNoNA_results.csv",
  "MedNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerMedNA_imputed_results.csv",
  "MedNoNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerMedNoNA_results.csv",
  "HighNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerHighNA_imputed_results.csv",
  "HighNoNA" = "./00-Data/03-Results/COPF/COPF_ProteoMakerHighNoNA_results.csv"
)

COPF_vs_groundtruth <- list()

for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  results <- fread(file_path)
  results_wide <- results[,
    .(adj_pval = min(peptide_adj_pval, na.rm = TRUE)),
    by = peptide_id
  ]
  ground_truth_name <- gsub("(NA|NoNA)$", "", name)
  ground_truth <- ground_truth_list[[ground_truth_name]]
  merged <- merge(
    ground_truth,
    results_wide,
    by.x = "Peptidoform",
    by.y = "peptide_id",
    all = TRUE
  )
  merged[is.na(adj_pval), adj_pval := 1]
  COPF_vs_groundtruth[[name]] <- merged
}

COPF_metrics <- data.table(
  Dataset = character(),
  Threshold = numeric(),
  TP = integer(),
  FP = integer(),
  TN = integer(),
  FN = integer(),
  FPR = numeric(),
  TPR = numeric(),
  FDR = numeric(),
  Precision = numeric(),
  Recall = numeric()
)

data_metrics <- data.table(
  thresholds = thresholds,
  TP = integer(length(thresholds)),
  FP = integer(length(thresholds)),
  TN = integer(length(thresholds)),
  FN = integer(length(thresholds)),
  FPR = numeric(length(thresholds)),
  TPR = numeric(length(thresholds)),
  FDR = numeric(length(thresholds)),
  Precision = numeric(length(thresholds)),
  Recall = numeric(length(thresholds))
)

for (name in names(COPF_vs_groundtruth)) {
  dataset_metrics <- copy(data_metrics)
  merged <- COPF_vs_groundtruth[[name]]
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    merged[, predicted_discordant := adj_pval < threshold]
    tp_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == TRUE])
    fp_val <- nrow(merged[isDiscordant == FALSE & predicted_discordant == TRUE])
    tn_val <- nrow(merged[
      isDiscordant == FALSE & predicted_discordant == FALSE
    ])
    fn_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == FALSE])
    fpr_val <- ifelse((fp_val + tn_val) > 0, fp_val / (fp_val + tn_val), NA)
    tpr_val <- ifelse((tp_val + fn_val) > 0, tp_val / (tp_val + fn_val), NA)
    fdr_val <- ifelse((tp_val + fp_val) > 0, fp_val / (tp_val + fp_val), NA)
    precision_val <- ifelse(
      (tp_val + fp_val) > 0,
      tp_val / (tp_val + fp_val),
      NA
    )
    recall_val <- tpr_val
    dataset_metrics[
      i,
      `:=`(
        Threshold = threshold,
        TP = tp_val,
        FP = fp_val,
        TN = tn_val,
        FN = fn_val,
        FPR = fpr_val,
        TPR = tpr_val,
        FDR = fdr_val,
        Precision = precision_val,
        Recall = recall_val
      )
    ]
  }
  dataset_metrics[, Dataset := name]
  COPF_metrics <- rbind(COPF_metrics, dataset_metrics, fill = TRUE)
}

fwrite(COPF_metrics, file.path(output_dir, "COPF_Metrics.csv"))

# ########################
# ## Metrics for RPC2-3 ##
# ########################
# # Note: peptides where p_adj == NA have p-value set to 1

# # Datasets
# datasets <- list(
#   "LowNA" = "./00-Data/03-Results/RPC/RPC2-3_ProteoMakerLowNA_results.csv",
#   "LowNoNA" = "./00-Data/03-Results/RPC/RPC2-3_ProteoMakerLowNoNA_results.csv",
#   "MedNA" = "./00-Data/03-Results/RPC/RPC2-3_ProteoMakerMedNA_results.csv",
#   "MedNoNA" = "./00-Data/03-Results/RPC/RPC2-3_ProteoMakerMedNoNA_results.csv",
#   "HighNA" = "./00-Data/03-Results/RPC/RPC2-3_ProteoMakerHighNA_results.csv",
#   "HighNoNA" = "./00-Data/03-Results/RPC/RPC2-3_ProteoMakerHighNoNA_results.csv"
# )
#
# RPC2_3_vs_groundtruth <- list()
#
# # Combine results with ground truth for each dataset
# for (name in names(datasets)) {
#   file_path <- datasets[[name]]
#   if (!file.exists(file_path)) {
#     warning(paste("File not found:", file_path))
#     next
#   }
#
#   # Read in results
#   results <- fread(file_path)
#
#   # Subset to relevant columns and rename p_adj to adj_pval
#   results_wide <- results[, .(
#     Peptide_ID,
#     adj_pval = p_adj
#   )]
#
#   # Remove everything before and including the last "_" in Peptide_ID to get Peptidoform
#   results_wide[,
#     Peptidoform := sub(".*_(?=[^_]+$)", "", Peptide_ID, perl = TRUE)
#   ]
#   results_wide[, Peptide_ID := NULL]
#
#   # Get corresponding ground truth table
#   ground_truth_name <- gsub("(NA|NoNA)$", "", name)
#   ground_truth <- ground_truth_list[[ground_truth_name]]
#
#   # Combine with ground truth by Peptidoform
#   merged <- merge(
#     ground_truth,
#     results_wide[, .(Peptidoform, adj_pval)],
#     by = "Peptidoform",
#     all = TRUE
#   )
#
#   # Set adj_pval to 1 for NA values (i.e., peptides not returned by RPC2-3)
#   merged[is.na(adj_pval), adj_pval := 1]
#
#   # Store in list
#   RPC2_3_vs_groundtruth[[name]] <- merged
# }
#
# RPC2_3_metrics <- data.table(
#   Dataset = character(),
#   Threshold = numeric(),
#   TP = integer(),
#   FP = integer(),
#   TN = integer(),
#   FN = integer(),
#   FPR = numeric(),
#   TPR = numeric(),
#   FDR = numeric(),
#   Precision = numeric(),
#   Recall = numeric()
# )
#
# data_metrics <- data.table(
#   thresholds = thresholds,
#   TP = integer(length(thresholds)),
#   FP = integer(length(thresholds)),
#   TN = integer(length(thresholds)),
#   FN = integer(length(thresholds)),
#   FPR = numeric(length(thresholds)),
#   TPR = numeric(length(thresholds)),
#   FDR = numeric(length(thresholds)),
#   Precision = numeric(length(thresholds)),
#   Recall = numeric(length(thresholds))
# )
#
# for (name in names(RPC2_3_vs_groundtruth)) {
#   dataset_metrics <- copy(data_metrics)
#
#   merged <- RPC2_3_vs_groundtruth[[name]]
#
#   for (i in seq_along(thresholds)) {
#     threshold <- thresholds[i]
#
#     # Classify peptides as discordant or not based on threshold
#     merged[, predicted_discordant := adj_pval < threshold]
#
#     # Calculate TP, FP, TN, FN
#     tp_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == TRUE])
#     fp_val <- nrow(merged[isDiscordant == FALSE & predicted_discordant == TRUE])
#     tn_val <- nrow(merged[
#       isDiscordant == FALSE & predicted_discordant == FALSE
#     ])
#     fn_val <- nrow(merged[isDiscordant == TRUE & predicted_discordant == FALSE])
#
#     # Calculate metrics
#     fpr_val <- ifelse((fp_val + tn_val) > 0, fp_val / (fp_val + tn_val), NA)
#     tpr_val <- ifelse((tp_val + fn_val) > 0, tp_val / (tp_val + fn_val), NA)
#     fdr_val <- ifelse((tp_val + fp_val) > 0, fp_val / (tp_val + fp_val), NA)
#     precision_val <- ifelse(
#       (tp_val + fp_val) > 0,
#       tp_val / (tp_val + fp_val),
#       NA
#     )
#     recall_val <- tpr_val
#
#     # Store metrics in data table
#     dataset_metrics[
#       i,
#       `:=`(
#         Threshold = threshold,
#         TP = tp_val,
#         FP = fp_val,
#         TN = tn_val,
#         FN = fn_val,
#         FPR = fpr_val,
#         TPR = tpr_val,
#         FDR = fdr_val,
#         Precision = precision_val,
#         Recall = recall_val
#       )
#     ]
#   }
#
#   # Add dataset name to metrics
#   dataset_metrics[, Dataset := name]
#
#   # Append to overall metrics table
#   RPC2_3_metrics <- rbind(RPC2_3_metrics, dataset_metrics, fill = TRUE)
# }
#
# # Save metrics to CSV
# fwrite(
#   RPC2_3_metrics,
#   file.path(output_dir, "RPC2_3_Metrics.csv")
# )

################################################################################
## Combined NA metrics — all methods, NA datasets only                        ##
################################################################################

all_metrics <- rbindlist(
  list(
    PeCorA_metrics[, method := "PeCorA"],
    ProteForge_metrics[, method := "ProteoForge"],
    DCF_metrics[, method := "ComplexoFinder"],
    DCFBaseline_metrics[, method := "ComplexoFinderBaseline"],
    COPF_metrics[, method := "COPF"]
  ),
  use.names = TRUE,
  fill = TRUE
)

# Compute F1 from confusion matrix counts
all_metrics[,
  F1 := ifelse(
    (2 * TP + FP + FN) > 0,
    2 * TP / (2 * TP + FP + FN),
    NA_real_
  )
]

# Simplify Dataset label: LowNA -> Low, MedNA -> Med, HighNA -> High
all_metrics[, NALevel := sub("NA$", "", Dataset)]

fwrite(all_metrics, file.path(output_dir, "Combined_NA_Metrics.csv"))
