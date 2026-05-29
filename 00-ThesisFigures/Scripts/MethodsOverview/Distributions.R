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

##################################################
## Filtering, transformation, and normalization ##
##################################################

# Boxplot before filtering, transformation, and normalization
data <- readRDS("ComplexoFinder/00-Data/H9_data_NoLogNoNorm.rds")

# Remove NMprotein from data
data <- data[!grepl("NMprotein", names(data))]

data <- rbindlist(data)

data$Accession <- NULL
data$`Gene name` <- NULL
data$`Position in master protein` <- NULL
data$`Modifications in master protein` <- NULL

data_long <- melt(
  data = data,
  id.vars = c("datatype"),
  variable.name = "sample",
  value.name = "intensity"
)

# Average per condition (ie D17, D19, D21 etc from H9_B1_D17, H9_B1_D19, H9_B1_D21) and datatype
data_long[, condition := sub("H9_B[1-3]_(D\\d+)", "\\1", sample)]

data_long$datatype <- factor(
  data_long$datatype,
  levels = c("NMpeptide", "LysAc", "Phospho", "Deglyco", "FreeCys", "RmCys")
)

data_plot <- data_long

data_plot$intensity <- log2(data_plot$intensity)

data_plot$typeCond <- paste(data_plot$datatype, data_plot$condition, sep = "_")

# all with the same color scheme, but different line types for conditions
data_plot$color <- datatype_colors[data_plot$datatype]

# Plot as density curves (one curve per condition per datatype)
density_noNorm <- ggplot(
  data_plot,
  aes(
    x = intensity,
    group = interaction(datatype, condition, drop = TRUE),
    color = datatype
  )
) +
  geom_density(fill = NA, linewidth = 0.4) +
  scale_color_manual(values = scales::alpha(datatype_colors, 0.5)) +
  scale_x_continuous(limits = c(3, 10)) +
  labs(
    title = "Intensities",
    x = "Log2(Intensity)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave2(
  filename = "00-ThesisFigures/Figures/Log2IntensityDistributions_NoNorm.pdf",
  plot = density_noNorm,
  width = 6,
  height = 4,
  dpi = 300
)

# Boxplot after filtering, transformation, and normalization
data <- readRDS("ComplexoFinder/00-Data/H9_data.rds")

# Remove NMprotein from data
data <- data[!grepl("NMprotein", names(data))]

data <- rbindlist(data)

data$Accession <- NULL
data$`Gene name` <- NULL
data$`Position in master protein` <- NULL
data$`Modifications in master protein` <- NULL

data_long <- melt(
  data = data,
  id.vars = c("datatype"),
  variable.name = "sample",
  value.name = "intensity"
)

# Average per condition (ie D17, D19, D21 etc from H9_B1_D17, H9_B1_D19, H9_B1_D21) and datatype
data_long[, condition := sub("H9_B[1-3]_(D\\d+)", "\\1", sample)]

data_long$datatype <- factor(
  data_long$datatype,
  levels = c("NMpeptide", "LysAc", "Phospho", "Deglyco", "FreeCys", "RmCys")
)

data_plot <- data_long

data_plot$typeCond <- paste(data_plot$datatype, data_plot$condition, sep = "_")

# all with the same color scheme, but different line types for conditions
data_plot$color <- datatype_colors[data_plot$datatype]

# Plot as density curves (one curve per condition per datatype)
density_norm <- ggplot(
  data_plot,
  aes(
    x = intensity,
    group = interaction(datatype, condition, drop = TRUE),
    color = datatype
  )
) +
  geom_density(fill = NA, linewidth = 0.4) +
  scale_color_manual(values = scales::alpha(datatype_colors, 0.5)) +
  scale_x_continuous(limits = c(-3.5, 3.5)) +
  labs(
    title = "Median Centered Intensities",
    x = "Log2(Intensity)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# add legend only for datatypes using cowplot only for the legend
legend <- get_legend(
  ggplot(data_plot, aes(x = intensity, y = c(1:10303272), color = datatype)) +
    geom_line(linewidth = 0.4) +
    scale_color_manual(values = scales::alpha(datatype_colors)) +
    labs(color = "Data Type") +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      nrow = 1
    )
)

ggsave2(
  filename = "00-ThesisFigures/Figures/Log2IntensityDistributions_Norm_Legend.pdf",
  plot = legend,
  width = 6,
  height = 2,
  dpi = 300
)

ggsave2(
  filename = "00-ThesisFigures/Figures/Log2IntensityDistributions_Norm.pdf",
  plot = density_norm,
  width = 6,
  height = 4,
  dpi = 300
)

density_norm

length(data_plot$intensity)
