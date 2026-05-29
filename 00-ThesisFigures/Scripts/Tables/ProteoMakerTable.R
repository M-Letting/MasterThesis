# Libraries
library(data.table)
library(gt)

# Set working directory
setwd(here::here())

# Output directory
output_dir <- "00-ThesisFigures/Tables"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Input directory
input_dir <- "Q-Benchmark/02-DCF/00-Data/01-Unprepared"

# Dataset definitions: name, PropModPerProt, PercDetectedVal
datasets <- list(
  list(name = "LowNoNA", prop_mod = 1, perc_detected = 1.0),
  list(name = "LowNA", prop_mod = 1, perc_detected = 0.5),
  list(name = "MedNoNA", prop_mod = 5, perc_detected = 1.0),
  list(name = "MedNA", prop_mod = 5, perc_detected = 0.5),
  list(name = "HighNoNA", prop_mod = 20, perc_detected = 1.0),
  list(name = "HighNA", prop_mod = 20, perc_detected = 0.5)
)

# Build summary rows
summary_rows <- lapply(datasets, function(d) {
  path <- file.path(input_dir, paste0("ProteoMaker_", d$name, ".csv"))
  df <- fread(path)

  # Intensity columns follow the pattern C_<cond>_R_<rep>
  intensity_cols <- grep("^C_[0-9]+_R_[0-9]+$", names(df), value = TRUE)

  n_peptides <- nrow(df)
  pct_missing <- round(
    mean(is.na(df[, ..intensity_cols])) * 100,
    digits = 1
  )

  data.table(
    Dataset = d$name,
    PropModPerProt = d$prop_mod,
    PercDetectedVal = d$perc_detected,
    `# Peptides` = n_peptides,
    `% Missing` = pct_missing
  )
})

summary_table <- rbindlist(summary_rows)

# Format table with gt
gt_table <- gt(summary_table) |>
  cols_label(
    Dataset = "Dataset",
    PropModPerProt = "Proteoforms Multiplier",
    PercDetectedVal = "Detection Rate",
    `# Peptides` = "# Peptides",
    `% Missing` = "% Missing"
  ) |>
  fmt_number(columns = "# Peptides", use_seps = TRUE, decimals = 0) |>
  fmt_number(columns = "% Missing", decimals = 1) |>
  fmt_number(columns = "PercDetectedVal", decimals = 1) |>
  cols_align(align = "right", columns = 2:5) |>
  cols_align(align = "left", columns = 1) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#cccccc", weight = px(1)),
    locations = list(
      cells_body(columns = 2:5),
      cells_column_labels(columns = 2:5)
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

gtsave(gt_table, file.path(output_dir, "ProteoMakerTable.png"))
