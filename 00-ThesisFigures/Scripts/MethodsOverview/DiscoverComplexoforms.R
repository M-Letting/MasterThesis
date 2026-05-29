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

######################################
## Discover complexoforms base idea ##
######################################

data <- readRDS(
  "ComplexoFinder/10-Results/B-WICH_chromatin_remodelling_complex/complexoform_results.rds"
)

data <- data$dCF$peptide_results

intensity_cols <- grep("H9", names(data), value = TRUE)

data <- data[, c("Identifier", intensity_cols), with = FALSE]

# Calculate average intensity per condition (ie D17, D19, D21 etc from H9_B1_D17, H9_B1_D19, H9_B1_D21)
intensity_cols <- grep("H9", names(data), value = TRUE)
data_long <- melt(
  data = data,
  id.vars = c("Identifier"),
  measure.vars = intensity_cols,
  variable.name = "sample",
  value.name = "intensity"
)

data_long[, condition := sub("H9_B[1-3]_D*", "\\1", sample)]


# Find median per condition
median_long <- data_long[,
  .(median_intensity = median(intensity, na.rm = TRUE)),
  by = c("condition")
]

median_long$condition <- as.numeric(median_long$condition)

# Smooth median profile with loess
loess_fit <- loess(median_intensity ~ condition, data = median_long, span = 0.4)
median_long$loess <- predict(loess_fit)


# --- choose two example peptides ---
example_peptides <- "BAZ1B_Q9UIG0_Q9UIG0 [114-126]_NA_RmCys"

# mean per condition for those peptides
example_mean <- data_long[
  Identifier %in% example_peptides,
  .(mean_intensity = mean(intensity, na.rm = TRUE)),
  by = .(Identifier, condition)
]

median_long$pep1 <- example_mean$mean_intensity

example_peptides <- "SF3B1_O75533_O75533 [319-333]_O75533 1xPhospho [T326(100)]_Phospho"

example_mean <- data_long[
  Identifier %in% example_peptides,
  .(mean_intensity = mean(intensity, na.rm = TRUE)),
  by = .(Identifier, condition)
]

median_long$pep2 <- example_mean$mean_intensity

# Smooth example peptides with loess
loess_fit_pep1 <- loess(pep1 ~ condition, data = median_long, span = 0.4)
median_long$pep1_loess <- predict(loess_fit_pep1)
loess_fit_pep2 <- loess(pep2 ~ condition, data = median_long, span = 0.4)
median_long$pep2_loess <- predict(loess_fit_pep2)

p <- ggplot(median_long) +
  geom_line(aes(x = condition, y = loess), linewidth = 1, color = "#3d3846ff") +
  geom_point(aes(x = condition, y = loess), size = 2, color = "#3d3846ff") +
  labs(
    title = "Testing for Discordant Peptide Profiles",
    x = "Day",
    y = "Intensity"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(20, 200, by = 20)) +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_line(
    aes(x = condition, y = pep1_loess),
    linewidth = 1,
    color = "#006d2c"
  ) +
  geom_point(aes(x = condition, y = pep1_loess), size = 2, color = "#006d2c") +
  geom_line(
    aes(x = condition, y = pep2_loess),
    linewidth = 1,
    color = "#be4e16ff"
  ) +
  geom_point(aes(x = condition, y = pep2_loess), size = 2, color = "#be4e16ff")

# Add legend saying "Median profile", "Cysteine modified peptide", "Phospho peptide" with corresponding colors
legend <- cowplot::get_legend(
  ggplot(median_long) +
    geom_line(
      aes(x = condition, y = loess, color = "Median profile"),
      linewidth = 1
    ) +
    geom_line(
      aes(x = condition, y = pep1_loess, color = "Cysteine modified peptide"),
      linewidth = 1
    ) +
    geom_line(
      aes(x = condition, y = pep2_loess, color = "Phospho peptide"),
      linewidth = 1
    ) +
    scale_color_manual(
      values = c(
        "Median profile" = "#3d3846ff",
        "Cysteine modified peptide" = "#006d2c",
        "Phospho peptide" = "#be4e16ff"
      )
    ) +
    labs(color = "Profile type") +
    theme(legend.title = element_blank(), legend.position = "bottom")
)
p <- cowplot::plot_grid(
  p + theme(legend.position = "none"),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.15)
)

p
# Save as pdf
ggsave(
  filename = file.path(output_dir, "DiscoverComplexoforms_PeptideProfiles.pdf"),
  plot = p,
  width = 6,
  height = 3,
  dpi = 300
)
