# Libraries
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(scico)
library(cowplot)

# Set working directory
setwd(paste0(here::here(), "/Q-Benchmark/03-Clustering"))

# Set output directory for plots
output_dir <- "00-Plots"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Membership threshold mode is expected to be set by 07-RunComparison.R
membership_threshold_mode <- getOption("qbenchmark.membership_threshold_mode")
if (is.null(membership_threshold_mode)) {
  stop(
    "Option 'qbenchmark.membership_threshold_mode' is not set. ",
    "Run this script through 07-RunComparison.R."
  )
}
membership_threshold_mode <- tolower(as.character(membership_threshold_mode))

include_nona <- isTRUE(getOption("qbenchmark.include_nona", TRUE))

filter_nona_datasets <- function(dt) {
  if (isTRUE(include_nona)) {
    return(dt)
  }
  if (!"Dataset" %in% names(dt)) {
    return(dt)
  }
  dt[!grepl("NoNA", Dataset)]
}

resolve_compare_file <- function(base_filename) {
  compare_dir <- "00-Data/05-CompareClustering"
  base_stem <- sub("\\.csv$", "", base_filename)

  candidate <- file.path(
    compare_dir,
    paste0(base_stem, "_", membership_threshold_mode, ".csv")
  )

  if (file.exists(candidate)) {
    return(candidate)
  }

  stop(
    "Could not find comparison input file for ",
    base_filename,
    " in mode '",
    membership_threshold_mode,
    "' (expected: ",
    candidate,
    ")"
  )
}

filter_rpc22_methods <- function(
  data,
  source_label = "input data",
  require_rows = TRUE
) {
  if (!"Method" %in% names(data)) {
    stop("Column 'Method' not found in ", source_label)
  }

  data <- data[!grepl("RPC", Method)]

  if (require_rows && nrow(data) == 0) {
    stop("No DCF-derived methods found in ", source_label)
  }

  data
}

################################################################################
## Plots of ProteoMaker datasets ###############################################
################################################################################

####################
# Helper functions #
####################

# Function to prepare data for plotting
prepare_data_PM_plot <- function(
  file_path,
  subtract_baseline = TRUE,
  baseline_method_prefix = "BlindDCF", # Prefix to identify baseline method
  y_var = "Weighted_Mean_Cluster_Purity" # Column to subtract baseline from
) {
  # Read data
  data <- fread(file_path)
  data <- filter_rpc22_methods(data, source_label = file_path)

  data <- filter_nona_datasets(data)

  # Harmonize dataset labels for matching across methods
  data[, Dataset_Key := gsub(".*_", "", Dataset)]

  if (subtract_baseline) {
    id_cols <- intersect(c("Dataset_Key", "Accession"), names(data))

    if (length(id_cols) == 0) {
      stop("Could not determine ID columns for baseline subtraction.")
    }

    baseline_methods <- c(
      baseline_method_prefix,
      paste0(baseline_method_prefix, "2-2")
    )

    data[,
      Baseline_Method := fifelse(
        grepl("2-2", Method),
        baseline_methods[2],
        baseline_methods[1]
      )
    ]

    baseline_data <- data[
      Method %in% baseline_methods,
      c(id_cols, "Method", y_var),
      with = FALSE
    ]

    setnames(
      baseline_data,
      old = c("Method", y_var),
      new = c("Baseline_Method", "Baseline_Value")
    )

    data <- merge(
      data,
      baseline_data,
      by = c(id_cols, "Baseline_Method"),
      all.x = TRUE,
      sort = FALSE
    )

    data[, (y_var) := get(y_var) - Baseline_Value]

    # Remove baseline methods from plot data
    data <- data[Method != Baseline_Method]

    data[, c("Baseline_Method", "Baseline_Value") := NULL]
  }

  # Keep normalized dataset labels for plotting
  data[, Dataset := Dataset_Key]
  data[, Dataset_Key := NULL]

  # Combine Method and ConstraintThreshold into a single column for plotting
  data[, Method := paste(Method, round(ConstraintThreshold, 2), sep = " ")]

  return(data)
}

# Function to create a base beeswarm plot for PM comparisons
create_PM_base_plot <- function(
  data,
  legend_position = "none",
  y_var,
  y_label,
  title
) {
  data <- data[is.finite(get(y_var))]
  # Create plot
  p <- ggplot(
    data,
    aes(x = Method, y = get(y_var), fill = Dataset)
  ) +
    geom_beeswarm(
      size = 3,
      shape = 21,
      alpha = 0.7,
      method = "compactswarm",
      corral = "gutter",
      corral.width = 0.7
    ) +
    stat_summary(
      data = data,
      mapping = aes(
        x = Method,
        y = get(y_var),
        group = Method,
        color = "Median"
      ),
      fun = median,
      geom = "point",
      inherit.aes = FALSE,
      size = 2.5,
      shape = 18
    ) +
    stat_summary(
      data = data,
      mapping = aes(
        x = Method,
        y = get(y_var),
        group = Method,
        color = "Mean"
      ),
      fun = mean,
      geom = "point",
      inherit.aes = FALSE,
      size = 2.5,
      shape = 18
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = legend_position,
      plot.title = element_text(hjust = 0.5),
      # Remove x-axis title
      axis.title.x = element_blank()
    ) +
    labs(
      title = title,
      y = y_label,
      fill = "Datasets",
      color = "Summary"
    ) +
    scale_fill_scico_d(
      palette = "batlow",
      na.value = "grey80"
    ) +
    scale_color_manual(
      values = c("Mean" = "black", "Median" = "red")
    )

  return(p)
}


###############
# Purity plot #
###############

# Prepare data for purity plot
data <- prepare_data_PM_plot(
  file_path = resolve_compare_file("ProteoMakerClusterPurityEvaluation.csv"),
  subtract_baseline = FALSE,
  y_var = "Weighted_Mean_Cluster_Purity"
)

# Plot
PM_Purity <- create_PM_base_plot(
  data = data,
  y_var = "Weighted_Mean_Cluster_Purity",
  y_label = "WMCP",
  title = "Weighted Mean Cluster Purity (WMCP)"
)

#################
# Fraction plot #
#################

# Prepare data for fraction plot
data <- prepare_data_PM_plot(
  file_path = resolve_compare_file("ProteoMakerFractionEvaluation.csv"),
  subtract_baseline = FALSE,
  y_var = "Top_Fraction"
)

# Create fraction plots
PM_Fraction <- create_PM_base_plot(
  data = data,
  y_var = "Top_Fraction",
  y_label = "Fraction",
  title = "Fraction of Proteoform ID 1 in Top Cluster"
)

#################################
# Sum unique plot (No baseline) #
#################################

# Prepare data for sum unique plot
data <- prepare_data_PM_plot(
  file_path = resolve_compare_file("ProteoMakerSumUniqueProteoforms.csv"),
  subtract_baseline = FALSE
)

# Plot
PM_SumUnique <- create_PM_base_plot(
  data = data,
  y_var = "Sum_Unique_Proteoforms",
  y_label = "Sum",
  title = "Sum of Unique Proteoforms in Top Cluster"
)

######################################
# Total potential plot (No baseline) #
######################################

# Prepare data for total potential plot
data <- prepare_data_PM_plot(
  file_path = resolve_compare_file("ProteoMakerTotalPotentialProteoforms.csv"),
  subtract_baseline = FALSE
)

# Set value column to numeric for log transformation
data[,
  Sum_Potential_Proteoforms := as.numeric(Sum_Potential_Proteoforms)
]

# Log10 transform value column
data[, Sum_Potential_Proteoforms := log10(Sum_Potential_Proteoforms)]

# Plot with log10 scale on y-axis
PM_TotalPotential <- create_PM_base_plot(
  data = data,
  y_var = "Sum_Potential_Proteoforms",
  y_label = "Log10(Total)",
  title = "Total Potential Proteoforms (Log10 Scale)"
)

#####################
# Weighted SSE plot #
#####################

# Prepare data for weighted SSE plot
data <- fread(
  resolve_compare_file("ProteoMakerWeightedSSEEvaluation.csv")
)
data <- filter_rpc22_methods(
  data,
  source_label = "ProteoMakerWeightedSSEEvaluation"
)

# Remove prefix from dataset labels
data[, Dataset := gsub(".*_", "", Dataset)]
data <- filter_nona_datasets(data)

# Combine Method and ConstraintThreshold into a single column for plotting
data[, Method := paste(Method, round(ConstraintThreshold, 2), sep = " ")]

# Plot
PM_WeightedSSE <- create_PM_base_plot(
  data = data,
  y_var = "Weighted_SSE",
  y_label = "Weighted SSE",
  title = "Weighted SSE"
)

# Prepare total SSE data for plotting
# Delete rows with no unique combination of Dataset, Accession, Method, and ConstraintThreshold
total_sse_data <- data
total_sse_data$combined <- paste(
  total_sse_data$Dataset,
  total_sse_data$Accession,
  total_sse_data$Method,
  total_sse_data$ConstraintThreshold
)

total_sse_data <- total_sse_data[!duplicated(total_sse_data$combined), ]

# Plot total SSE
PM_TotalSSE <- create_PM_base_plot(
  data = total_sse_data,
  y_var = "Total_Weighted_SSE",
  y_label = "Total Weighted SSE",
  title = "Total Weighted SSE Per Unique Dataset-Accession-Method-Constraint Combination"
)

##########################################
# Combine plots into one figure and save #
##########################################

combined_plot <- plot_grid(
  PM_Purity,
  PM_Fraction,
  PM_SumUnique,
  PM_TotalPotential,
  PM_WeightedSSE,
  PM_TotalSSE,
  ncol = 2,
  labels = "AUTO"
)

# Add shared legend
legend <- get_legend(
  PM_Purity +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 1))
)

final_PM_plot <- plot_grid(
  combined_plot,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

# Save the final plot
ggsave(
  filename = file.path(
    output_dir,
    paste0("ProteoMakerComparisonPlots_", membership_threshold_mode, ".pdf")
  ),
  plot = final_PM_plot,
  width = 16,
  height = 9,
  dpi = 300
)

################################################################################
## Plots of Protein Complex datasets ###########################################
################################################################################

#################
# Fraction plot #
#################

# Load and prepare data
data <- prepare_data_PM_plot(
  file_path = resolve_compare_file("ComplexFractionEvaluation.csv"),
  subtract_baseline = FALSE,
  y_var = "Top_Fraction"
)

# Plot
Complex_Fraction <- create_PM_base_plot(
  data = data,
  y_var = "Top_Fraction",
  legend_position = "none",
  y_label = "Fraction",
  title = "Fraction of NMPeptide in Top Cluster"
)

#############################
# Complex Weighted SSE plot #
#############################

# Load and prepare data
data <- fread(
  resolve_compare_file("ComplexWeightedSSEEvaluation.csv")
)
data <- filter_rpc22_methods(
  data,
  source_label = "ComplexWeightedSSEEvaluation"
)

# Remove prefix from dataset labels
data[, Dataset := gsub(".*_", "", Dataset)]
data <- filter_nona_datasets(data)

# Combine Method and ConstraintThreshold into a single column for plotting
data[, Method := paste(Method, round(ConstraintThreshold, 2), sep = " ")]

# Plot
Complex_WeightedSSE <- create_PM_base_plot(
  data = data,
  y_var = "Weighted_SSE",
  legend_position = "none",
  y_label = "Weighted SSE",
  title = "Weighted SSE"
)

# Prepare total SSE data for plotting
# Delete rows with no unique combination of Dataset, Accession, Method, and ConstraintThreshold
total_sse_data <- data
total_sse_data$combined <- paste(
  total_sse_data$Dataset,
  total_sse_data$Accession,
  total_sse_data$Method,
  total_sse_data$ConstraintThreshold
)

total_sse_data <- total_sse_data[!duplicated(total_sse_data$combined), ]

# Plot total SSE
Complex_TotalSSE <- create_PM_base_plot(
  data = total_sse_data,
  y_var = "Total_Weighted_SSE",
  legend_position = "none",
  y_label = "Total Weighted SSE",
  title = "Total Weighted SSE Per Unique Dataset-Accession-Method-Constraint Combination"
)

##########################################
# Combine plots into one figure and save #
##########################################

combined_plot <- plot_grid(
  Complex_Fraction,
  Complex_WeightedSSE,
  Complex_TotalSSE,
  ncol = 1,
  labels = "AUTO"
)

# Add shared legend
legend <- get_legend(
  Complex_Fraction +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 3))
)
final_Complex_plot <- plot_grid(
  combined_plot,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.15)
)


ggsave(
  filename = file.path(
    output_dir,
    paste0("ComplexComparisonPlots_", membership_threshold_mode, ".pdf")
  ),
  plot = final_Complex_plot,
  width = 16,
  height = 9,
  dpi = 300
)

################################################################################
## Mean and median comparisons for all metrics #################################
################################################################################

metric_specs <- data.table(
  Metric = c(
    "PM_PurityMM",
    "PM_FractionMM",
    "PM_SumUniqueMM",
    "PM_TotalPotentialMM",
    "PM_WeightedSSEMM",
    "PM_TotalSSEMM",
    "Complex_FractionMM",
    "Complex_WeightedSSEMM",
    "Complex_TotalSSEMM"
  ),
  file_path = c(
    resolve_compare_file("ProteoMakerClusterPurityEvaluation.csv"),
    resolve_compare_file("ProteoMakerFractionEvaluation.csv"),
    resolve_compare_file("ProteoMakerSumUniqueProteoforms.csv"),
    resolve_compare_file("ProteoMakerTotalPotentialProteoforms.csv"),
    resolve_compare_file("ProteoMakerWeightedSSEEvaluation.csv"),
    resolve_compare_file("ProteoMakerWeightedSSEEvaluation.csv"),
    resolve_compare_file("ComplexFractionEvaluation.csv"),
    resolve_compare_file("ComplexWeightedSSEEvaluation.csv"),
    resolve_compare_file("ComplexWeightedSSEEvaluation.csv")
  ),
  y_var = c(
    "Weighted_Mean_Cluster_Purity",
    "Top_Fraction",
    "Sum_Unique_Proteoforms",
    "Sum_Potential_Proteoforms",
    "Weighted_SSE",
    "Total_Weighted_SSE",
    "Top_Fraction",
    "Weighted_SSE",
    "Total_Weighted_SSE"
  ),
  transform = c(
    "identity",
    "identity",
    "identity",
    "log10",
    "identity",
    "identity",
    "identity",
    "identity",
    "identity"
  )
)

summary_results <- rbindlist(
  lapply(seq_len(nrow(metric_specs)), function(i) {
    metric <- metric_specs$Metric[i]
    file_path <- metric_specs$file_path[i]
    y_var <- metric_specs$y_var[i]
    transform <- metric_specs$transform[i]

    data <- prepare_data_PM_plot(
      file_path = file_path,
      subtract_baseline = FALSE,
      y_var = y_var
    )

    if (!y_var %in% names(data)) {
      stop("Column '", y_var, "' not found for metric '", metric, "'.")
    }

    data[, value := as.numeric(get(y_var))]

    if (transform == "log10") {
      data[, value := log10(value)]
    }

    data[
      is.finite(value),
      .(
        Mean = mean(value, na.rm = TRUE),
        Median = median(value, na.rm = TRUE)
      ),
      by = .(Method)
    ][, `:=`(Metric = metric, Transform = transform)][]
  }),
  use.names = TRUE,
  fill = TRUE
)

setcolorder(
  summary_results,
  c("Metric", "Method", "Mean", "Median", "Transform")
)

# Create plots for mean and median comparisons for each metric
for (metric in unique(summary_results$Metric)) {
  metric_data <- summary_results[Metric == metric]
  transform <- metric_data$Transform[1]

  p <- ggplot(
    metric_data,
    aes(x = Method, y = Mean)
  ) +
    geom_point(size = 3, color = "black", shape = 18) +
    geom_point(aes(y = Median), size = 3, color = "red", shape = 18) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank()
    ) +
    labs(
      y = fifelse(transform == "log10", "Log10(Value)", "Value")
    )

  assign(paste0(metric, "_MeanMedianPlot"), p)
}

# Create combined figure with mean/median plots + swarm plots
PM_Fraction_combined <- plot_grid(
  PM_Fraction + theme(legend.position = "none") + labs(title = NULL),
  PM_FractionMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
# Add title to combined plot
PM_Fraction_combined <- plot_grid(
  ggdraw() + draw_label("Fraction of Proteoform ID 1 in Top Cluster"),
  PM_Fraction_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

PM_Purity_combined <- plot_grid(
  PM_Purity + theme(legend.position = "none") + labs(title = NULL),
  PM_PurityMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
# Add title to combined plot
PM_Purity_combined <- plot_grid(
  ggdraw() + draw_label("Weighted Mean Cluster Purity (WMCP)"),
  PM_Purity_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

PM_SumUnique_combined <- plot_grid(
  PM_SumUnique + theme(legend.position = "none") + labs(title = NULL),
  PM_SumUniqueMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
# Add title to combined plot
PM_SumUnique_combined <- plot_grid(
  ggdraw() + draw_label("Sum of Unique Proteoforms in Top Cluster"),
  PM_SumUnique_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

PM_TotalPotential_combined <- plot_grid(
  PM_TotalPotential +
    theme(legend.position = "none") +
    labs(title = NULL),
  PM_TotalPotentialMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
# Add title to combined plot
PM_TotalPotential_combined <- plot_grid(
  ggdraw() + draw_label("Total Potential Proteoforms (Log10 Scale)"),
  PM_TotalPotential_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

PM_WeightedSSE_combined <- plot_grid(
  PM_WeightedSSE + theme(legend.position = "none") + labs(title = NULL),
  PM_WeightedSSEMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
# Add title to combined plot
PM_WeightedSSE_combined <- plot_grid(
  ggdraw() + draw_label("Weighted SSE"),
  PM_WeightedSSE_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

PM_TotalSSE_combined <- plot_grid(
  PM_TotalSSE + theme(legend.position = "none") + labs(title = NULL),
  PM_TotalSSEMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
# Add title to combined plot
PM_TotalSSE_combined <- plot_grid(
  ggdraw() +
    draw_label(
      "Total Weighted SSE Per Unique Dataset-Accession-Method-Constraint Combination"
    ),
  PM_TotalSSE_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

PM_combined <- plot_grid(
  PM_Purity_combined,
  PM_Fraction_combined,
  PM_SumUnique_combined,
  PM_TotalPotential_combined,
  PM_WeightedSSE_combined,
  PM_TotalSSE_combined,
  ncol = 2,
  labels = "AUTO"
)

# Add shared legend at the bottom
legend <- get_legend(
  PM_Purity +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 1))
)
final_PM_combined <- plot_grid(
  PM_combined,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

ggsave(
  filename = file.path(
    output_dir,
    paste0(
      "ProteoMakerComparisonPlots_CombinedWithMeanMedian_",
      membership_threshold_mode,
      ".pdf"
    )
  ),
  plot = final_PM_combined,
  width = 16,
  height = 9,
  dpi = 300
)

Complex_Fraction_combined <- plot_grid(
  Complex_Fraction +
    theme(legend.position = "none") +
    labs(title = NULL),
  Complex_FractionMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
Complex_Fraction_combined <- plot_grid(
  ggdraw() + draw_label("Fraction of NMPeptide in Top Cluster"),
  Complex_Fraction_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

Complex_WeightedSSE_combined <- plot_grid(
  Complex_WeightedSSE +
    theme(legend.position = "none") +
    labs(title = NULL),
  Complex_WeightedSSEMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
Complex_WeightedSSE_combined <- plot_grid(
  ggdraw() + draw_label("Weighted SSE"),
  Complex_WeightedSSE_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

Complex_TotalSSE_combined <- plot_grid(
  Complex_TotalSSE +
    theme(legend.position = "none") +
    labs(title = NULL),
  Complex_TotalSSEMM_MeanMedianPlot,
  nrow = 1,
  rel_widths = c(1, 0.5)
)
Complex_TotalSSE_combined <- plot_grid(
  ggdraw() +
    draw_label(
      "Total Weighted SSE Per Unique Dataset-Accession-Method-Constraint Combination"
    ),
  Complex_TotalSSE_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

Complex_combined <- plot_grid(
  Complex_Fraction_combined,
  Complex_WeightedSSE_combined,
  Complex_TotalSSE_combined,
  ncol = 1,
  labels = "AUTO"
)

# Add shared legend at the bottom
legend <- get_legend(
  Complex_Fraction +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 3))
)
final_Complex_combined <- plot_grid(
  Complex_combined,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)
final_Complex_combined
ggsave(
  filename = file.path(
    output_dir,
    paste0(
      "ComplexComparisonPlots_CombinedWithMeanMedian_",
      membership_threshold_mode,
      ".pdf"
    )
  ),
  plot = final_Complex_combined,
  width = 16,
  height = 9,
  dpi = 300
)

################################################################################
## Compare position only potential proteoforms #################################
################################################################################

plot_data <- prepare_data_PM_plot(
  file_path = resolve_compare_file(
    "ProteoMakerPositionPotentialProteoforms.csv"
  ),
  subtract_baseline = FALSE,
  y_var = "Sum_Position_Proteoforms"
)

# Remove features in dataset HighNA or HighNoNA and accession Q12150
plot_data <- plot_data[
  !(Dataset %in% c("HighNA", "HighNoNA") & Accession == "Q12150"),
]

plot_data$Sum_Position_Proteoforms <- log10(plot_data$Sum_Position_Proteoforms)

PM_position <- create_PM_base_plot(
  data = plot_data,
  y_var = "Sum_Position_Proteoforms",
  y_label = "Log10(Sum)",
  title = "Sum of Position-Only Potential Proteoforms (Log10 Scale)"
)

plot_data <- prepare_data_PM_plot(
  file_path = resolve_compare_file("ComplexPositionPotentialProteoforms.csv"),
  subtract_baseline = FALSE,
  y_var = "Sum_Position_Proteoforms"
)

Complex_position <- create_PM_base_plot(
  data = plot_data,
  y_var = "Sum_Position_Proteoforms",
  y_label = "Sum",
  title = "Sum of Position-Only Potential Proteoforms"
)

position_Combined <- plot_grid(
  PM_position + theme(legend.position = "none") + labs(title = NULL),
  Complex_position + theme(legend.position = "none") + labs(title = NULL),
  ncol = 1,
  labels = "AUTO"
)

# Add shared legend at the bottom and title to combined plot
legend <- get_legend(
  PM_position +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 1))
)
final_position_combined <- plot_grid(
  ggdraw() + draw_label("Sum of Position-Only Potential Proteoforms "),
  position_Combined,
  legend,
  ncol = 1,
  rel_heights = c(0.05, 1, 0.05)
)

ggsave(
  filename = file.path(
    output_dir,
    paste0(
      "PositionPotentialProteoforms_Comparison_",
      membership_threshold_mode,
      ".pdf"
    )
  ),
  plot = final_position_combined,
  width = 16,
  height = 9,
  dpi = 300
)
