# Libraries
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(ggplot2, warn.conflicts = FALSE))
suppressMessages(library(scico, warn.conflicts = FALSE))
suppressMessages(library(cowplot, warn.conflicts = FALSE))

# Change working directory to Q-Benchmark/02-DCF
setwd(here::here("Q-Benchmark/02-DCF"))


# Define output directory for comparison results
output_dir <- "./00-Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load metrics for all methods
PeCora_metrics <- fread("./00-Data/04-PerformanceMetrics/PeCora_Metrics.csv")
ProteoForge_metrics <- fread(
  "./00-Data/04-PerformanceMetrics/ProteoForge_Metrics.csv"
)
DCF_metrics <- fread("./00-Data/04-PerformanceMetrics/DCF_Metrics.csv")
DCFBaseline_metrics <- fread(
  "./00-Data/04-PerformanceMetrics/DCFBaseline_Metrics.csv"
)
RPC_metrics <- fread("./00-Data/04-PerformanceMetrics/RPC_Metrics.csv")
COPF_metrics <- fread("./00-Data/04-PerformanceMetrics/COPF_Metrics.csv")

# Combine all metrics into a single data table for plotting
# COPF is only included for NoNA datasets
all_metrics <- rbind(
  PeCora_metrics[, Method := "PeCorA"],
  ProteoForge_metrics[, Method := "ProteoForge"],
  DCF_metrics[, Method := "DCF"],
  DCFBaseline_metrics[, Method := "DCF Baseline"],
  RPC_metrics[, Method := "RPC"],
  COPF_metrics[Dataset %in% c("LowNoNA", "MedNoNA", "HighNoNA")][, Method := "COPF"],
  fill = TRUE
)

# Set color palette for plotting
color_palette <- scico(palette = "managua", n = 6)

# Create color mapping for methods
method_colors <- setNames(
  color_palette,
  c(
    "DCF",
    "DCF Baseline",
    "RPC",
    "PeCorA",
    "ProteoForge",
    "COPF"
  )
)

.calc_auc_trapz <- function(fpr, tpr) {
  ord <- order(fpr)
  fpr <- fpr[ord]
  tpr <- tpr[ord]
  sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2, na.rm = TRUE)
}

unique_thresholds <- sort(unique(all_metrics$Threshold), decreasing = FALSE)
unique_thresholds <- unique_thresholds[is.finite(unique_thresholds)]
x_breaks <- unique_thresholds[seq(1, length(unique_thresholds), by = 3)]
x_label_breaks <- c(0.001, 0.05, 0.25, 0.5, 0.75)

n_thr <- length(unique_thresholds)
n1 <- min(9, n_thr)
n2 <- min(9, max(n_thr - n1, 0))
n3 <- max(n_thr - n1 - n2, 0)

axis_positions <- numeric(n_thr)
if (n1 > 0) {
  axis_positions[1:n1] <- seq(0, 1 / 3, length.out = n1)
}
if (n2 > 0) {
  axis_positions[(n1 + 1):(n1 + n2)] <- seq(1 / 3, 2 / 3, length.out = n2)
}
if (n3 > 0) {
  axis_positions[(n1 + n2 + 1):n_thr] <- seq(2 / 3, 1, length.out = n3)
}

threshold_axis_trans <- scales::trans_new(
  name = "threshold_axis_trans",
  transform = function(x) {
    approx(
      x = unique_thresholds,
      y = axis_positions,
      xout = x,
      rule = 2,
      ties = "ordered"
    )$y
  },
  inverse = function(y) {
    approx(
      x = axis_positions,
      y = unique_thresholds,
      xout = y,
      rule = 2,
      ties = "ordered"
    )$y
  }
)

################################
## Plotting function for TPR ##
################################
plot_TPR <- function(metrics_dt, dataset) {
  if (!dataset %in% metrics_dt$Dataset) {
    stop("Dataset not found in metrics data table: ", dataset)
  }

  plot_data <- metrics_dt[Dataset == dataset]

  method_order <- names(method_colors)
  plot_data[, Method := factor(Method, levels = method_order)]
  plot_data <- plot_data[order(Method, Threshold)]

  ggplot(plot_data) +
    geom_line(
      aes(
        x = Threshold,
        y = TPR,
        color = Method,
        group = Method
      ),
      linewidth = 1.2
    ) +
    geom_point(
      aes(
        x = Threshold,
        y = TPR,
        fill = Method,
        shape = Method
      ),
      size = 2.6,
      color = "black",
      stroke = 0.8
    ) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    scale_shape_manual(
      values = c(
        "DCF" = 23,
        "PeCorA" = 22,
        "ProteoForge" = 24,
        "RPC" = 21,
        "DCF Baseline" = 25,
        "COPF" = 8
      )
    ) +
    scale_x_continuous(
      trans = threshold_axis_trans,
      breaks = x_label_breaks,
      labels = c("0.001", "0.05", "0.25", "0.5", "0.75"),
      limits = range(unique_thresholds),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = paste(dataset),
      x = "Adjusted p-value",
      y = "TPR"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    ) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey50")
}

################################
## Plotting function for FPR ##
################################
plot_FPR <- function(metrics_dt, dataset) {
  if (!dataset %in% metrics_dt$Dataset) {
    stop("Dataset not found in metrics data table: ", dataset)
  }

  plot_data <- metrics_dt[Dataset == dataset]

  method_order <- names(method_colors)
  plot_data[, Method := factor(Method, levels = method_order)]
  plot_data <- plot_data[order(Method, Threshold)]

  ggplot(plot_data) +
    geom_line(
      aes(
        x = Threshold,
        y = FPR,
        color = Method,
        group = Method
      ),
      linewidth = 1.2
    ) +
    geom_point(
      aes(
        x = Threshold,
        y = FPR,
        fill = Method,
        shape = Method
      ),
      size = 2.6,
      color = "black",
      stroke = 0.8
    ) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    scale_shape_manual(
      values = c(
        "DCF" = 23,
        "PeCorA" = 22,
        "ProteoForge" = 24,
        "RPC" = 21,
        "DCF Baseline" = 25,
        "COPF" = 8
      )
    ) +
    scale_x_continuous(
      trans = threshold_axis_trans,
      breaks = x_label_breaks,
      labels = c("0.001", "0.05", "0.25", "0.5", "0.75"),
      limits = range(unique_thresholds),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = paste(dataset),
      x = "Adjusted p-value",
      y = "FPR"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    ) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey50")
}

###################################
## Plotting function for Recall ##
###################################
plot_Recall <- function(metrics_dt, dataset) {
  if (!dataset %in% metrics_dt$Dataset) {
    stop("Dataset not found in metrics data table: ", dataset)
  }

  plot_data <- metrics_dt[Dataset == dataset]

  if (!"Recall" %in% names(plot_data)) {
    stop("metrics_dt must contain 'Recall' column")
  }

  method_order <- names(method_colors)
  plot_data[, Method := factor(Method, levels = method_order)]
  plot_data <- plot_data[order(Method, Threshold)]

  ggplot(plot_data) +
    geom_line(
      aes(
        x = Threshold,
        y = Recall,
        color = Method,
        group = Method
      ),
      linewidth = 1.2
    ) +
    geom_point(
      aes(
        x = Threshold,
        y = Recall,
        fill = Method,
        shape = Method
      ),
      size = 2.6,
      color = "black",
      stroke = 0.8
    ) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    scale_shape_manual(
      values = c(
        "DCF" = 23,
        "PeCorA" = 22,
        "ProteoForge" = 24,
        "RPC" = 21,
        "DCF Baseline" = 25,
        "COPF" = 8
      )
    ) +
    scale_x_continuous(
      trans = threshold_axis_trans,
      breaks = x_label_breaks,
      labels = c("0.001", "0.05", "0.25", "0.5", "0.75"),
      limits = range(unique_thresholds),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = paste(dataset),
      x = "Adjusted p-value",
      y = "Recall"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    ) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey50")
}

####################################
## Plotting function for Precision ##
####################################
plot_Precision <- function(metrics_dt, dataset) {
  if (!dataset %in% metrics_dt$Dataset) {
    stop("Dataset not found in metrics data table: ", dataset)
  }

  plot_data <- metrics_dt[Dataset == dataset]

  if (!"Precision" %in% names(plot_data)) {
    stop("metrics_dt must contain 'Precision' column")
  }

  method_order <- names(method_colors)
  plot_data[, Method := factor(Method, levels = method_order)]
  plot_data <- plot_data[order(Method, Threshold)]

  ggplot(plot_data) +
    geom_line(
      aes(
        x = Threshold,
        y = Precision,
        color = Method,
        group = Method
      ),
      linewidth = 1.2
    ) +
    geom_point(
      aes(
        x = Threshold,
        y = Precision,
        fill = Method,
        shape = Method
      ),
      size = 2.6,
      color = "black",
      stroke = 0.8
    ) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    scale_shape_manual(
      values = c(
        "DCF" = 23,
        "PeCorA" = 22,
        "ProteoForge" = 24,
        "RPC" = 21,
        "DCF Baseline" = 25,
        "COPF" = 8
      )
    ) +
    scale_x_continuous(
      trans = threshold_axis_trans,
      breaks = x_label_breaks,
      labels = c("0.001", "0.05", "0.25", "0.5", "0.75"),
      limits = range(unique_thresholds),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = paste(dataset),
      x = "Adjusted p-value",
      y = "Precision"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    ) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey50")
}

###############################
## Plotting function for FDR ##
###############################
plot_FDR <- function(metrics_dt, dataset) {
  if (!dataset %in% metrics_dt$Dataset) {
    stop("Dataset not found in metrics data table: ", dataset)
  }

  plot_data <- metrics_dt[Dataset == dataset]

  if (!"FDR" %in% names(plot_data)) {
    stop("metrics_dt must contain 'FDR' column")
  }

  method_order <- names(method_colors)
  plot_data[, Method := factor(Method, levels = method_order)]
  plot_data <- plot_data[order(Method, Threshold)]

  ggplot(plot_data) +
    geom_line(
      aes(
        x = Threshold,
        y = FDR,
        color = Method,
        group = Method
      ),
      linewidth = 1.2
    ) +
    geom_point(
      aes(
        x = Threshold,
        y = FDR,
        fill = Method,
        shape = Method
      ),
      size = 2.6,
      color = "black",
      stroke = 0.8
    ) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    scale_shape_manual(
      values = c(
        "DCF" = 23,
        "PeCorA" = 22,
        "ProteoForge" = 24,
        "RPC" = 21,
        "DCF Baseline" = 25,
        "COPF" = 8
      )
    ) +
    scale_x_continuous(
      trans = threshold_axis_trans,
      breaks = x_label_breaks,
      labels = c("0.001", "0.05", "0.25", "0.5", "0.75"),
      limits = range(unique_thresholds),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = paste(dataset),
      x = "Adjusted p-value",
      y = "FDR"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    ) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey50")
}

#####################################
## Creating plots for all datasets ##
#####################################

# TPR plots
tpr1 <- plot_TPR(all_metrics, "LowNA")
tpr2 <- plot_TPR(all_metrics, "MedNA")
tpr3 <- plot_TPR(all_metrics, "HighNA")
tpr4 <- plot_TPR(all_metrics, "LowNoNA")
tpr5 <- plot_TPR(all_metrics, "MedNoNA")
tpr6 <- plot_TPR(all_metrics, "HighNoNA")

# FPR plots
fpr1 <- plot_FPR(all_metrics, "LowNA")
fpr2 <- plot_FPR(all_metrics, "MedNA")
fpr3 <- plot_FPR(all_metrics, "HighNA")
fpr4 <- plot_FPR(all_metrics, "LowNoNA")
fpr5 <- plot_FPR(all_metrics, "MedNoNA")
fpr6 <- plot_FPR(all_metrics, "HighNoNA")

# Recall plots
recall1 <- plot_Recall(all_metrics, "LowNA")
recall2 <- plot_Recall(all_metrics, "MedNA")
recall3 <- plot_Recall(all_metrics, "HighNA")
recall4 <- plot_Recall(all_metrics, "LowNoNA")
recall5 <- plot_Recall(all_metrics, "MedNoNA")
recall6 <- plot_Recall(all_metrics, "HighNoNA")

# Precision plots
prec1 <- plot_Precision(all_metrics, "LowNA")
prec2 <- plot_Precision(all_metrics, "MedNA")
prec3 <- plot_Precision(all_metrics, "HighNA")
prec4 <- plot_Precision(all_metrics, "LowNoNA")
prec5 <- plot_Precision(all_metrics, "MedNoNA")
prec6 <- plot_Precision(all_metrics, "HighNoNA")

# FDR plots
fdr1 <- plot_FDR(all_metrics, "LowNA")
fdr2 <- plot_FDR(all_metrics, "MedNA")
fdr3 <- plot_FDR(all_metrics, "HighNA")
fdr4 <- plot_FDR(all_metrics, "LowNoNA")
fdr5 <- plot_FDR(all_metrics, "MedNoNA")
fdr6 <- plot_FDR(all_metrics, "HighNoNA")

# Combined TPR plots
combined_plot_TPR <- plot_grid(
  tpr1,
  tpr2,
  tpr3,
  # tpr4,
  # tpr5,
  # tpr6,
  ncol = 3
)

# Helper to build a legend from a given metrics subset
.make_legend <- function(metrics_dt) {
  shape_values <- c(
    "DCF" = 23, "PeCorA" = 22, "ProteoForge" = 24,
    "RPC" = 21, "DCF Baseline" = 25, "COPF" = 8
  )
  get_legend(
    ggplot(metrics_dt) +
      geom_line(
        aes(x = Threshold, y = TPR, color = Method, group = Method),
        linewidth = 1.2
      ) +
      geom_point(
        aes(x = Threshold, y = TPR, fill = Method, shape = Method),
        size = 2.6, color = "black", stroke = 0.8
      ) +
      scale_color_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      scale_shape_manual(values = shape_values) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "lines")
      )
  )
}

# NA plots legend (no COPF); NoNA plots legend (includes COPF)
legend_plot <- .make_legend(all_metrics[Method != "COPF"])
legend_plot_noNA <- .make_legend(all_metrics)

combined_plot_TPR_with_legend <- plot_grid(
  combined_plot_TPR,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

# Combined FPR plots
combined_plot_FPR <- plot_grid(
  fpr1,
  fpr2,
  fpr3,
  # fpr4,
  # fpr5,
  # fpr6,
  ncol = 3
)

combined_plot_FPR_with_legend <- plot_grid(
  combined_plot_FPR,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.05)
)


# Combined Recall plots
combined_plot_Recall <- plot_grid(
  recall1,
  recall2,
  recall3,
  # recall4,
  # recall5,
  # recall6,
  ncol = 3
)

combined_plot_Recall_with_legend <- plot_grid(
  combined_plot_Recall,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

# Combined Precision plots
combined_plot_Precision <- plot_grid(
  prec1,
  prec2,
  prec3,
  # prec4,
  # prec5,
  # prec6,
  ncol = 3
)

combined_plot_Precision_with_legend <- plot_grid(
  combined_plot_Precision,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

# Combined FDR plots
combined_plot_FDR <- plot_grid(
  fdr1,
  fdr2,
  fdr3,
  # fdr4,
  # fdr5,
  # fdr6,
  ncol = 3
)

combined_plot_FDR_with_legend <- plot_grid(
  combined_plot_FDR,
  legend_plot,
  ncol = 1,
  scale = 0.95,
  rel_heights = c(1, 0.1),
  labels = c("E)", ""),
  label_size = 20,
  label_fontface = "bold"
)

all_combined1 <- plot_grid(
  combined_plot_TPR,
  combined_plot_FPR,
  combined_plot_Recall,
  combined_plot_Precision,
  scale = 0.95,
  ncol = 2,
  labels = c("A)", "B)", "C)", "D)"),
  label_size = 20,
  label_fontface = "bold"
)

all_combined2 <- plot_grid(
  all_combined1,
  combined_plot_FDR_with_legend,
  ncol = 1,
  rel_heights = c(1, 0.5),
  rel_widths = c(1, 1),
  labels = c("", "E)"),
  label_size = 20,
  label_fontface = "bold"
)

# Make equivalent plots with LowNoNA, MedNoNA, HighNoNA datasets
combined_plot_TPR_noNA <- plot_grid(
  tpr4,
  tpr5,
  tpr6,
  ncol = 3
)
combined_plot_FPR_noNA <- plot_grid(
  fpr4,
  fpr5,
  fpr6,
  ncol = 3
)
combined_plot_Recall_noNA <- plot_grid(
  recall4,
  recall5,
  recall6,
  ncol = 3
)
combined_plot_Precision_noNA <- plot_grid(
  prec4,
  prec5,
  prec6,
  ncol = 3
)
combined_plot_FDR_noNA <- plot_grid(
  fdr4,
  fdr5,
  fdr6,
  ncol = 3
)

combined_plot_FDR_noNA_with_legend <- plot_grid(
  combined_plot_FDR_noNA,
  legend_plot_noNA,
  scale = 0.95,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

combined1_noNA <- plot_grid(
  combined_plot_TPR_noNA,
  combined_plot_FPR_noNA,
  combined_plot_Recall_noNA,
  combined_plot_Precision_noNA,
  scale = 0.95,
  ncol = 2,
  labels = c("A)", "B)", "C)", "D)"),
  label_size = 20,
  label_fontface = "bold"
)

combined2_noNA <- plot_grid(
  combined1_noNA,
  combined_plot_FDR_noNA_with_legend,
  ncol = 1,
  rel_heights = c(1, 0.5),
  rel_widths = c(1, 1),
  labels = c("", "E)"),
  label_size = 20,
  label_fontface = "bold"
)

# Save the combined plots
ggsave(
  filename = file.path(output_dir, "06-CompareDiscoveryMethods-NA.pdf"),
  plot = all_combined2,
  width = 16,
  height = 9,
  dpi = 300
)
ggsave(
  filename = file.path(output_dir, "06-CompareDiscoveryMethods-NoNA.pdf"),
  plot = combined2_noNA,
  width = 16,
  height = 9,
  dpi = 300
)
