# Libraries
library(data.table)
library(gt)

# Set working directory
setwd(here::here())

# Define input and output paths
input <- "ComplexoFinder/10-Results/Eukaryotic_translation_initiation_factor_3_complex/cluster_summary.csv"
output_dir <- "00-ThesisFigures/Tables"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load data
data <- fread(input)

data <- data[Datatype == "Phospho"]
data <- data[`membership of cluster 6` > 0.5]

# Extract peptide range from Position (strip leading accession prefix and brackets)
data[, PeptidePosition := sub("^[^ ]+ \\[(.*)\\]$", "\\1", Position)]

# Parse Modification column to extract phosphosite count and residue locations
# Example format: "Q14152 1xPhospho [S929(100)]"
data[, NumPhosphosites := as.integer(
  sub(".*?(\\d+)xPhospho.*", "\\1", Modification)
)]

data[, PhosphositeLocation := gsub(
  "\\([^)]*\\)", "",                          # remove confidence scores in ()
  regmatches(Modification, regexpr("\\[.*\\]", Modification))
)]
data[, PhosphositeLocation := gsub("[\\[\\]]", "", PhosphositeLocation)]
data[, PhosphositeLocation := trimws(PhosphositeLocation)]

# Build final table
final_table <- data[, .(
  Gene        = `Gene name`,
  Position    = PeptidePosition,
  `# Phosphosites` = NumPhosphosites,
  `Phosphosite Location` = PhosphositeLocation
)]

setorder(final_table, Gene, Position)

# Create gt table matching thesis style
gt_table <- gt(final_table) |>
  cols_label(
    Gene                   = "Gene",
    Position               = "Peptide Position",
    `# Phosphosites`       = "# Phosphosites",
    `Phosphosite Location` = "Phosphosite Location"
  ) |>
  cols_align(align = "left",  columns = c(Gene, Position, `Phosphosite Location`)) |>
  cols_align(align = "right", columns = `# Phosphosites`) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_borders(sides = "left", color = "#cccccc", weight = px(1)),
    locations = list(
      cells_body(columns = 2:4),
      cells_column_labels(columns = 2:4)
    )
  ) |>
  tab_style(
    style     = cell_fill(color = "#f5f5f5"),
    locations = cells_body(rows = seq(1, nrow(final_table), by = 2))
  ) |>
  tab_style(
    style     = cell_fill(color = "#e0e0e0"),
    locations = cells_column_labels()
  ) |>
  tab_options(
    table.border.top.style           = "hidden",
    table.border.bottom.style        = "hidden",
    column_labels.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black",
    table.font.names                  = "Arial"
  )

gtsave(gt_table, file.path(output_dir, "EIF3PhosphoSummary.png"))
