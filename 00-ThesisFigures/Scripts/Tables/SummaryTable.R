# libraries
library(data.table)
library(gt)

# Set working directory
setwd(here::here())

# Set output directory for figures
output_dir <- "00-ThesisFigures/Tables"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define input
NMpeptide <- fread("ComplexoFinder/00-Data/NonModified_Peptide.csv")
LysAC <- fread("ComplexoFinder/00-Data/AcetylatedLysines_Peptide.csv")
Deglyco <- fread("ComplexoFinder/00-Data/Deglycosylated_Peptide.csv")
FreeCys <- fread("ComplexoFinder/00-Data/FreeCysteines_Peptide.csv")
Phospho <- fread("ComplexoFinder/00-Data/Phospho_Peptide.csv")
RmCys <- fread("ComplexoFinder/00-Data/ReversiblyModifiedCysteines_Peptide.csv")

# Create H9 and IMR90 specific datasets
metadata_cols <- c(
  "Accession",
  "Gene name",
  "Position in master protein",
  "Modifications in master protein",
  "datatype"
)
H9_cols <- grep("H9", names(NMpeptide), value = TRUE)
IMR90_cols <- grep("IMR90", names(NMpeptide), value = TRUE)

# H9 dataset with metadata
H9_NMpeptide <- NMpeptide[, c(metadata_cols, H9_cols), with = FALSE]
H9_LysAC <- LysAC[, c(metadata_cols, H9_cols), with = FALSE]
H9_Deglyco <- Deglyco[, c(metadata_cols, H9_cols), with = FALSE]
H9_FreeCys <- FreeCys[, c(metadata_cols, H9_cols), with = FALSE]
H9_Phospho <- Phospho[, c(metadata_cols, H9_cols), with = FALSE]
H9_RmCys <- RmCys[, c(metadata_cols, H9_cols), with = FALSE]

# Remove rows with NA in all intensity columns
H9_NMpeptide <- H9_NMpeptide[
  rowSums(is.na(H9_NMpeptide[, H9_cols, with = FALSE])) < length(H9_cols)
]
H9_LysAC <- H9_LysAC[
  rowSums(is.na(H9_LysAC[, H9_cols, with = FALSE])) < length(H9_cols)
]
H9_Deglyco <- H9_Deglyco[
  rowSums(is.na(H9_Deglyco[, H9_cols, with = FALSE])) < length(H9_cols)
]
H9_FreeCys <- H9_FreeCys[
  rowSums(is.na(H9_FreeCys[, H9_cols, with = FALSE])) < length(H9_cols)
]
H9_RmCys <- H9_RmCys[
  rowSums(is.na(H9_RmCys[, H9_cols, with = FALSE])) < length(H9_cols)
]
H9_Phospho <- H9_Phospho[
  rowSums(is.na(H9_Phospho[, H9_cols, with = FALSE])) < length(H9_cols)
]

# IMR90 dataset with metadata
IMR90_NMpeptide <- NMpeptide[, c(metadata_cols, IMR90_cols), with = FALSE]
IMR90_LysAC <- LysAC[, c(metadata_cols, IMR90_cols), with = FALSE]
IMR90_Deglyco <- Deglyco[, c(metadata_cols, IMR90_cols), with = FALSE]
IMR90_FreeCys <- FreeCys[, c(metadata_cols, IMR90_cols), with = FALSE]
IMR90_Phospho <- Phospho[, c(metadata_cols, IMR90_cols), with = FALSE]
IMR90_RmCys <- RmCys[, c(metadata_cols, IMR90_cols), with = FALSE]

# Remove rows with NA in all intensity columns
IMR90_NMpeptide <- IMR90_NMpeptide[
  rowSums(is.na(IMR90_NMpeptide[, IMR90_cols, with = FALSE])) <
    length(IMR90_cols)
]
IMR90_LysAC <- IMR90_LysAC[
  rowSums(is.na(IMR90_LysAC[, IMR90_cols, with = FALSE])) < length(IMR90_cols)
]
IMR90_Deglyco <- IMR90_Deglyco[
  rowSums(is.na(IMR90_Deglyco[, IMR90_cols, with = FALSE])) < length(IMR90_cols)
]
IMR90_FreeCys <- IMR90_FreeCys[
  rowSums(is.na(IMR90_FreeCys[, IMR90_cols, with = FALSE])) < length(IMR90_cols)
]
IMR90_RmCys <- IMR90_RmCys[
  rowSums(is.na(IMR90_RmCys[, IMR90_cols, with = FALSE])) < length(IMR90_cols)
]
IMR90_Phospho <- IMR90_Phospho[
  rowSums(is.na(IMR90_Phospho[, IMR90_cols, with = FALSE])) < length(IMR90_cols)
]

# Create summary table
summary_table <- data.table(
  Datatype = c(
    "NMpeptide",
    "Phospho",
    "LysAc",
    "Deglyco",
    "FreeCys",
    "RmCys"
  ),
  Combined = c(
    nrow(NMpeptide),
    nrow(Phospho),
    nrow(LysAC),
    nrow(Deglyco),
    nrow(FreeCys),
    nrow(RmCys)
  ),
  H9 = c(
    nrow(H9_NMpeptide),
    nrow(H9_Phospho),
    nrow(H9_LysAC),
    nrow(H9_Deglyco),
    nrow(H9_FreeCys),
    nrow(H9_RmCys)
  ),
  IMR90 = c(
    nrow(IMR90_NMpeptide),
    nrow(IMR90_Phospho),
    nrow(IMR90_LysAC),
    nrow(IMR90_Deglyco),
    nrow(IMR90_FreeCys),
    nrow(IMR90_RmCys)
  )
)

names(summary_table) <- c(
  "Datatype",
  "# Features Combined",
  "# Features H9",
  "# Features IMR90"
)

# Create and save publication-quality table with gt
gt_table <- gt(summary_table) |>
  cols_label(
    Datatype = "Data Type",
    `# Features H9` = "# Features H9",
    `# Features IMR90` = "# Features IMR90",
    `# Features Combined` = "# Features Total"
  ) |>
  cols_align(align = "right", columns = 2:4) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#cccccc", weight = px(1)),
    locations = list(
      cells_body(columns = 2:4),
      cells_column_labels(columns = 2:4)
    )
  ) |>
  tab_style(
    style = cell_fill(color = "#f5f5f5"),
    locations = cells_body(rows = seq(1, nrow(summary_table), by = 2))
  ) |>
  tab_style(
    style = cell_fill(color = "#e0e0e0"),
    locations = cells_column_labels()
  ) |>
  tab_options(
    table.border.top.style = "hidden",
    table.border.bottom.style = "hidden",
    column_labels.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black",
    table.font.names = "Arial"
  )

gtsave(gt_table, file.path(output_dir, "SummaryTable.pdf"))

# Filter tables to only include peptides with > 40% non-NA values in combined, H9, and IMR90 datasets
filter_threshold <- 0.4

# Filter tables based on the threshold
H9_NMpeptide <- H9_NMpeptide[
  rowMeans(!is.na(H9_NMpeptide[, H9_cols, with = FALSE])) > filter_threshold
]
H9_LysAC <- H9_LysAC[
  rowMeans(!is.na(H9_LysAC[, H9_cols, with = FALSE])) > filter_threshold
]
H9_Deglyco <- H9_Deglyco[
  rowMeans(!is.na(H9_Deglyco[, H9_cols, with = FALSE])) > filter_threshold
]
H9_FreeCys <- H9_FreeCys[
  rowMeans(!is.na(H9_FreeCys[, H9_cols, with = FALSE])) > filter_threshold
]
H9_Phospho <- H9_Phospho[
  rowMeans(!is.na(H9_Phospho[, H9_cols, with = FALSE])) > filter_threshold
]
H9_RmCys <- H9_RmCys[
  rowMeans(!is.na(H9_RmCys[, H9_cols, with = FALSE])) > filter_threshold
]

IMR90_NMpeptide <- IMR90_NMpeptide[
  rowMeans(!is.na(IMR90_NMpeptide[, IMR90_cols, with = FALSE])) >
    filter_threshold
]
IMR90_LysAC <- IMR90_LysAC[
  rowMeans(!is.na(IMR90_LysAC[, IMR90_cols, with = FALSE])) > filter_threshold
]
IMR90_Deglyco <- IMR90_Deglyco[
  rowMeans(!is.na(IMR90_Deglyco[, IMR90_cols, with = FALSE])) > filter_threshold
]
IMR90_FreeCys <- IMR90_FreeCys[
  rowMeans(!is.na(IMR90_FreeCys[, IMR90_cols, with = FALSE])) > filter_threshold
]
IMR90_Phospho <- IMR90_Phospho[
  rowMeans(!is.na(IMR90_Phospho[, IMR90_cols, with = FALSE])) > filter_threshold
]
IMR90_RmCys <- IMR90_RmCys[
  rowMeans(!is.na(IMR90_RmCys[, IMR90_cols, with = FALSE])) > filter_threshold
]

# Update summary table with filtered counts
summary_filtered <- data.table(
  Datatype = c(
    "NMpeptide",
    "Phospho",
    "LysAc",
    "Deglyco",
    "FreeCys",
    "RmCys"
  ),
  H9 = c(
    nrow(H9_NMpeptide),
    nrow(H9_Phospho),
    nrow(H9_LysAC),
    nrow(H9_Deglyco),
    nrow(H9_FreeCys),
    nrow(H9_RmCys)
  ),
  IMR90 = c(
    nrow(IMR90_NMpeptide),
    nrow(IMR90_Phospho),
    nrow(IMR90_LysAC),
    nrow(IMR90_Deglyco),
    nrow(IMR90_FreeCys),
    nrow(IMR90_RmCys)
  )
)

names(summary_filtered) <- c(
  "Datatype",
  "# Features H9 (filtered)",
  "# Features IMR90 (filtered)"
)

# Merge with original summary table
summary_final <- merge(
  summary_table,
  summary_filtered,
  by = "Datatype",
  sort = FALSE
)

# Add total row
summary_final <- rbind(
  summary_final,
  data.table(
    Datatype = "Total",
    `# Features Combined` = sum(summary_final$`# Features Combined`),
    `# Features H9` = sum(summary_final$`# Features H9`),
    `# Features IMR90` = sum(summary_final$`# Features IMR90`),
    `# Features H9 (filtered)` = sum(summary_final$`# Features H9 (filtered)`),
    `# Features IMR90 (filtered)` = sum(
      summary_final$`# Features IMR90 (filtered)`
    )
  )
)

# Create and save publication-quality table with gt
gt_table_final <- gt(summary_final) |>
  cols_label(
    Datatype = "Data Type",
    `# Features H9` = "# Features H9",
    `# Features IMR90` = "# Features IMR90",
    `# Features Combined` = "# Features Total",
    `# Features H9 (filtered)` = "# Features H9 (filtered)",
    `# Features IMR90 (filtered)` = "# Features IMR90 (filtered)"
  ) |>
  cols_align(align = "right", columns = 2:6) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#cccccc", weight = px(1)),
    locations = list(
      cells_body(columns = 2:6),
      cells_column_labels(columns = 2:6)
    )
  ) |>
  tab_style(
    style = cell_fill(color = "#f5f5f5"),
    locations = cells_body(rows = seq(1, nrow(summary_final), by = 2))
  ) |>
  tab_style(
    style = cell_fill(color = "#e0e0e0"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_borders(sides = "top", color = "black", weight = px(2)),
    locations = cells_body(rows = nrow(summary_final))
  ) |>
  tab_options(
    table.border.top.style = "hidden",
    table.border.bottom.style = "hidden",
    column_labels.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black",
    table.font.names = "Arial"
  )

gtsave(gt_table_final, file.path(output_dir, "SummaryTable_Filtered.png"))

# Add total row to no-filtered table
summary_table <- rbind(
  summary_table,
  data.table(
    Datatype = "Total",
    `# Features H9` = sum(summary_table$`# Features H9`),
    `# Features IMR90` = sum(summary_table$`# Features IMR90`),
    `# Features Combined` = sum(summary_table$`# Features Combined`)
  )
)

# Create table without the filtered columns
gt_table_no_filtered <- gt(summary_table) |>
  cols_label(
    Datatype = "Data Type",
    `# Features H9` = "# Features H9",
    `# Features IMR90` = "# Features IMR90",
    `# Features Combined` = "# Features Total"
  ) |>
  cols_align(align = "right", columns = 2:4) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#cccccc", weight = px(1)),
    locations = list(
      cells_body(columns = 2:4),
      cells_column_labels(columns = 2:4)
    )
  ) |>
  tab_style(
    style = cell_fill(color = "#f5f5f5"),
    locations = cells_body(rows = seq(1, nrow(summary_table), by = 2))
  ) |>
  tab_style(
    style = cell_fill(color = "#e0e0e0"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_borders(sides = "top", color = "black", weight = px(2)),
    locations = cells_body(rows = nrow(summary_table))
  ) |>
  tab_options(
    table.border.top.style = "hidden",
    table.border.bottom.style = "hidden",
    column_labels.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black",
    table.font.names = "Arial"
  )
gtsave(
  gt_table_no_filtered,
  file.path(output_dir, "SummaryTable_NoFiltered.png")
)
