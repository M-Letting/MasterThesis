library(data.table)
library(gt)
library(scico)
library(arrow)

setwd(here::here())

# ── Constants ─────────────────────────────────────────────────────────────────

TARGET_THRESHOLD <- 0.01

method_levels <- c(
  "ComplexoFinder",
  "ComplexoFinderBaseline",
  "ProteoForge",
  "PeCorA",
  "COPF"
)
method_labels <- c(
  "ComplexoFinder" = "ComplexoFinder",
  "ComplexoFinderBaseline" = "ComplexoFinder Baseline",
  "ProteoForge" = "ProteoForge",
  "PeCorA" = "PeCorA",
  "COPF" = "COPF"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

rename_methods <- function(dt) {
  dt[method == "DCF", method := "ComplexoFinder"]
  dt[method == "DCFBaseline", method := "ComplexoFinderBaseline"]
  dt
}

# For each method × condition, pick the row closest to TARGET_THRESHOLD,
# then return the mean of `metric` across conditions.
mean_at_threshold <- function(dt, threshold_col, metric_col, group_cols) {
  at_thr <- dt[,
    .SD[which.min(abs(get(threshold_col) - TARGET_THRESHOLD))],
    by = group_cols
  ]
  at_thr[, .(val = mean(get(metric_col), na.rm = TRUE)), by = method]
}

# Build a named vector over method_levels; NA where method is absent.
to_row <- function(means_dt, fmt = "%.2f") {
  out <- setNames(rep(NA_real_, length(method_levels)), method_levels)
  for (i in seq_len(nrow(means_dt))) {
    m <- means_dt$method[i]
    if (m %in% method_levels) out[m] <- means_dt$val[i]
  }
  ifelse(is.na(out), "—", sprintf(fmt, out))
}

# ── 1. SWATH-MS Detection ─────────────────────────────────────────────────────

swath_id <- fread(
  "Q-Benchmark/01-FindDCF/data/results/peptide_identification_performance_data.csv"
)
rename_methods(swath_id)

swath_id_tpr <- mean_at_threshold(
  swath_id,
  "threshold",
  "TPR",
  c("method", "perturbation")
)
swath_id_fpr <- mean_at_threshold(
  swath_id,
  "threshold",
  "FPR",
  c("method", "perturbation")
)

# ── 2. ProteoMaker Detection ──────────────────────────────────────────────────

pm_id <- fread(
  "Q-Benchmark/02-DCF/00-Data/04-PerformanceMetrics/Combined_NA_Metrics.csv"
)
rename_methods(pm_id)

# Per-dataset helper: value at target threshold per method for one specific dataset
pm_det_val <- function(dataset_name, metric) {
  sub <- pm_id[Dataset == dataset_name]
  at_thr <- sub[,
    .SD[which.min(abs(Threshold - TARGET_THRESHOLD))],
    by = method
  ]
  at_thr[, .(method, val = get(metric))]
}

# ── 3. SWATH-MS Grouping ──────────────────────────────────────────────────────

swath_grp <- fread(
  "Q-Benchmark/01-FindDCF/data/results/peptide_grouping_performance_data.csv"
)
rename_methods(swath_grp)

swath_grp_tpr <- mean_at_threshold(
  swath_grp,
  "threshold",
  "TPR",
  c("method", "perturbation")
)
swath_grp_fpr <- mean_at_threshold(
  swath_grp,
  "threshold",
  "FPR",
  c("method", "perturbation")
)

# ── 4. ProteoMaker Grouping — MAD ─────────────────────────────────────────────
# MAD = mean absolute deviation of (found proteoform count − ProteoMaker count)
# Only ComplexoFinder, ComplexoFinderBaseline, and ProteoForge ran this benchmark.

base_dir <- "Q-Benchmark/02-DCF"
ref_files <- list(
  LowNoNA = file.path(
    base_dir,
    "00-Data/02-Prepared/ProteoMakerLowNoNA_wide_processed.csv"
  ),
  MedNoNA = file.path(
    base_dir,
    "00-Data/02-Prepared/ProteoMakerMedNoNA_wide_processed.csv"
  ),
  HighNoNA = file.path(
    base_dir,
    "00-Data/02-Prepared/ProteoMakerHighNoNA_wide_processed.csv"
  )
)
dcf_files <- list(
  LowNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerLowNA_results.csv"
  ),
  MedNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerMedNA_results.csv"
  ),
  HighNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerHighNA_results.csv"
  ),
  LowNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerLowNoNA_results.csv"
  ),
  MedNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerMedNoNA_results.csv"
  ),
  HighNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCF/DCF_ProteoMakerHighNoNA_results.csv"
  )
)
dcfbl_files <- list(
  LowNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNA_results.csv"
  ),
  MedNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNA_results.csv"
  ),
  HighNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNA_results.csv"
  ),
  LowNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerLowNoNA_results.csv"
  ),
  MedNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerMedNoNA_results.csv"
  ),
  HighNoNA = file.path(
    base_dir,
    "00-Data/03-Results/DCFBaseline/DCFBaseline_ProteoMakerHighNoNA_results.csv"
  )
)
pf_files <- list(
  LowNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNA_result.feather"
  ),
  MedNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNA_result.feather"
  ),
  HighNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNA_result.feather"
  ),
  LowNoNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerLowNoNA_result.feather"
  ),
  MedNoNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerMedNoNA_result.feather"
  ),
  HighNoNA = file.path(
    base_dir,
    "00-Data/03-Results/ProteoForge/ProteoForge_ProteoMakerHighNoNA_result.feather"
  )
)

# Reference key: for NA datasets use the corresponding NoNA file
ref_key <- function(name) {
  if (grepl("NoNA", name)) name else sub("NA$", "NoNA", name)
}

count_ref_proteoforms <- function(ref_data, prot) {
  ids <- as.character(ref_data[Accession == prot, Proteoform_ID])
  ids <- unlist(strsplit(ids, "\\|"), use.names = FALSE)
  ids <- ids[ids != "" & !is.na(ids)]
  length(unique(ids))
}

# Load found-vs-reference counts for one method across all datasets
load_dcf_counts <- function(file_list, method_name) {
  rows <- rbindlist(lapply(names(file_list), function(nm) {
    path <- file_list[[nm]]
    if (!file.exists(path)) {
      return(NULL)
    }
    res <- fread(path)
    grps <- res[!is.na(dCF) & dCF != "dCF-1", .N, by = .(Accession, dCF)]
    found <- grps[N > 1, .(Found = .N), by = Accession]
    found <- found[!grepl("\\|", Accession)]
    ref <- fread(ref_files[[ref_key(nm)]])
    found[, Ref := sapply(Accession, count_ref_proteoforms, ref_data = ref)]
    found[, .(Dataset = nm, Method = method_name, diff = Found - Ref)]
  }))
  rows
}

load_pf_counts <- function(file_list) {
  rows <- rbindlist(lapply(names(file_list), function(nm) {
    path <- file_list[[nm]]
    if (!file.exists(path)) {
      return(NULL)
    }
    pf_res <- as.data.table(read_feather(path))[Sample == "C1_R1"]
    grps <- pf_res[, .(n = .N), by = .(protein_id, ClusterID)]
    found <- grps[n > 1, .(Found = .N), by = protein_id]
    found <- found[!grepl("\\|", protein_id)]
    ref <- fread(ref_files[[ref_key(nm)]])
    found[, Ref := sapply(protein_id, count_ref_proteoforms, ref_data = ref)]
    found[, .(Dataset = nm, Method = "ProteoForge", diff = Found - Ref)]
  }))
  rows
}

compute_mad <- function(diff_vec) {
  mean(abs(diff_vec - mean(diff_vec, na.rm = TRUE)), na.rm = TRUE)
}

mapping_counts <- rbindlist(list(
  load_dcf_counts(dcf_files, "ComplexoFinder"),
  load_dcf_counts(dcfbl_files, "ComplexoFinderBaseline"),
  load_pf_counts(pf_files)
))

# Per-dataset MAD helper
pm_grp_mad_ds <- function(dataset_name) {
  out <- mapping_counts[
    Dataset == dataset_name,
    .(val = compute_mad(diff)),
    by = Method
  ]
  setnames(out, "Method", "method")
  out
}

# ── Assemble summary table ────────────────────────────────────────────────────

make_row <- function(benchmark, metric, means_dt) {
  vals <- to_row(means_dt)
  data.table(
    Benchmark = benchmark,
    Metric = metric,
    ComplexoFinder = vals["ComplexoFinder"],
    ComplexoFinderBaseline = vals["ComplexoFinderBaseline"],
    ProteoForge = vals["ProteoForge"],
    PeCorA = vals["PeCorA"],
    COPF = vals["COPF"]
  )
}

table_data <- rbindlist(list(
  # SWATH-MS Detection (mean across perturbations)
  make_row("SWATH-MS - Detection", "TPR (p = 0.01)", swath_id_tpr),
  make_row("SWATH-MS - Detection", "FPR (p = 0.01)", swath_id_fpr),
  # SWATH-MS Grouping (mean across perturbations)
  make_row("SWATH-MS - Grouping", "TPR (p = 0.01)", swath_grp_tpr),
  make_row("SWATH-MS - Grouping", "FPR (p = 0.01)", swath_grp_fpr),
  # ProteoMaker Detection — Low
  make_row(
    "ProteoMaker Detection - Low",
    "TPR (50% NA)",
    pm_det_val("LowNA", "TPR")
  ),
  make_row(
    "ProteoMaker Detection - Low",
    "FDR (50% NA)",
    pm_det_val("LowNA", "FDR")
  ),
  make_row(
    "ProteoMaker Detection - Low",
    "TPR (No NA)",
    pm_det_val("LowNoNA", "TPR")
  ),
  make_row(
    "ProteoMaker Detection - Low",
    "FDR (No NA)",
    pm_det_val("LowNoNA", "FDR")
  ),
  # ProteoMaker Detection — Medium
  make_row(
    "ProteoMaker Detection - Medium",
    "TPR (50% NA)",
    pm_det_val("MedNA", "TPR")
  ),
  make_row(
    "ProteoMaker Detection - Medium",
    "FDR (50% NA)",
    pm_det_val("MedNA", "FDR")
  ),
  make_row(
    "ProteoMaker Detection - Medium",
    "TPR (No NA)",
    pm_det_val("MedNoNA", "TPR")
  ),
  make_row(
    "ProteoMaker Detection - Medium",
    "FDR (No NA)",
    pm_det_val("MedNoNA", "FDR")
  ),
  # ProteoMaker Detection — High
  make_row(
    "ProteoMaker Detection - High",
    "TPR (50% NA)",
    pm_det_val("HighNA", "TPR")
  ),
  make_row(
    "ProteoMaker Detection - High",
    "FDR (50% NA)",
    pm_det_val("HighNA", "FDR")
  ),
  make_row(
    "ProteoMaker Detection - High",
    "TPR (No NA)",
    pm_det_val("HighNoNA", "TPR")
  ),
  make_row(
    "ProteoMaker Detection - High",
    "FDR (No NA)",
    pm_det_val("HighNoNA", "FDR")
  ),
  # ProteoMaker Grouping — Low
  make_row(
    "ProteoMaker Grouping - Low",
    "MAD (50% NA)",
    pm_grp_mad_ds("LowNA")
  ),
  make_row(
    "ProteoMaker Grouping - Low",
    "MAD (No NA)",
    pm_grp_mad_ds("LowNoNA")
  ),
  # ProteoMaker Grouping — Medium
  make_row(
    "ProteoMaker Grouping - Medium",
    "MAD (50% NA)",
    pm_grp_mad_ds("MedNA")
  ),
  make_row(
    "ProteoMaker Grouping - Medium",
    "MAD (No NA)",
    pm_grp_mad_ds("MedNoNA")
  ),
  # ProteoMaker Grouping — High
  make_row(
    "ProteoMaker Grouping - High",
    "MAD (50% NA)",
    pm_grp_mad_ds("HighNA")
  ),
  make_row(
    "ProteoMaker Grouping - High",
    "MAD (No NA)",
    pm_grp_mad_ds("HighNoNA")
  )
))

# ── gt table ──────────────────────────────────────────────────────────────────

gt_tbl <- table_data |>
  gt(groupname_col = "Benchmark") |>
  cols_label(
    Metric = "Metric",
    ComplexoFinder = "ComplexoFinder",
    ComplexoFinderBaseline = "ComplexoFinder Baseline",
    ProteoForge = "ProteoForge",
    PeCorA = "PeCorA",
    COPF = "COPF"
  ) |>
  tab_header(
    title = md("**Benchmark Performance Summary**"),
    subtitle = "Values at threshold p = 0.01"
  ) |>
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_row_groups()
  ) |>
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = list(cell_fill(color = "#f0f4f8"), cell_text(weight = "bold")),
    locations = cells_column_labels(columns = ComplexoFinder)
  ) |>
  tab_style(
    style = list(cell_fill(color = "#f0f4f8")),
    locations = cells_body(columns = ComplexoFinder)
  ) |>
  cols_width(
    ComplexoFinder ~ px(140),
    ComplexoFinderBaseline ~ px(140),
    ProteoForge ~ px(140),
    PeCorA ~ px(140),
    COPF ~ px(140)
  ) |>
  cols_align(align = "center", columns = -Metric) |>
  tab_options(
    table.font.size = px(13),
    row_group.font.weight = "bold",
    heading.align = "left"
  )

print(gt_tbl)

# ── Save ──────────────────────────────────────────────────────────────────────

gtsave(gt_tbl, "00-ThesisFigures/Figures/Tables/BenchmarkSummaryTable.html")
gtsave(gt_tbl, "00-ThesisFigures/Figures/Tables/BenchmarkSummaryTable.pdf")
gtsave(gt_tbl, "00-ThesisFigures/Figures/Tables/BenchmarkSummaryTable.png")

cat(
  "Saved: 00-ThesisFigures/Figures/Tables/BenchmarkSummaryTable.html/.pdf/.png\n"
)

fwrite(table_data, "00-ThesisFigures/Figures/Tables/BenchmarkSummaryTable.csv")

# ── Qualitative summary table ─────────────────────────────────────────────────

summary_data <- data.table(
  Scenario = c(
    "Simple perturbation detection",
    "Complex proteoform detection",
    "Proteoform count estimation",
    "Increasing proteoform complexity"
  ),
  `Best-performing method` = c(
    "ProteoForge",
    "ComplexoFinder",
    "ComplexoFinder",
    "ComplexoFinder"
  ),
  Reason = c(
    "Better specificity",
    "Better recovery of underlying proteoform structures",
    "Lowest MAD",
    "Relative performance improves as complexity increases"
  )
)

gt_summary <- summary_data |>
  gt() |>
  tab_header(title = md("**Benchmark — Qualitative Summary**")) |>
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = `Best-performing method`,
      rows = `Best-performing method` == "ComplexoFinder"
    )
  ) |>
  cols_width(
    Scenario ~ px(220),
    `Best-performing method` ~ px(200),
    Reason ~ px(300)
  ) |>
  cols_align(align = "left") |>
  tab_options(
    table.font.size = px(13),
    heading.align = "left"
  )

print(gt_summary)

gtsave(
  gt_summary,
  "00-ThesisFigures/Figures/Tables/BenchmarkQualitativeSummary.html"
)
gtsave(
  gt_summary,
  "00-ThesisFigures/Figures/Tables/BenchmarkQualitativeSummary.pdf"
)
gtsave(
  gt_summary,
  "00-ThesisFigures/Figures/Tables/BenchmarkQualitativeSummary.png"
)

cat(
  "Saved: 00-ThesisFigures/Figures/Tables/BenchmarkQualitativeSummary.html/.pdf/.png\n"
)
