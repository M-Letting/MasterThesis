# Libraries
library(data.table)

# Set directories
setwd(here::here())

output_dir <- "Q-Benchmark/02-DCF/00-Data/02-Prepared"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

source("ComplexoFinder/02-Identification/IdentifyComplexes.R")

################################################################################
## 1. Load pre-processed data                                                 ##
################################################################################

H9_LogNorm <- readRDS("ComplexoFinder/00-Data/H9_data.rds")
H9_Raw <- readRDS("ComplexoFinder/00-Data/H9_data_NoLogNoNorm.rds")

################################################################################
## 2. Load lookup tables                                                      ##
################################################################################

EBI_LT <- readRDS("ComplexoFinder/00-LookupTables/ebi_cp_lookup.rds")

################################################################################
## 3. Identify and prune complexes                                            ##
################################################################################

build_universe <- function(data_list) {
  accessions <- unlist(lapply(data_list, function(dt) {
    unlist(lapply(dt[["Accession"]], function(x) {
      if (is.na(x) || x == "") return(character(0))
      trimws(strsplit(x, ";")[[1]])
    }))
  }))
  accessions <- sub("-1$", "", trimws(accessions))
  unique(accessions[nzchar(accessions) & !is.na(accessions)])
}

find_and_prune <- function(data_list) {
  universe <- build_universe(data_list)

  cx <- find_complexes(
    data_list = data_list,
    lookup_table = EBI_LT,
    id_col_data = "Accession",
    id_col_lookup = "uniprot_id"
  )
  cx_pruned <- prune_complexes_stat(
    all_complexes = cx$complexes,
    lookup_table = EBI_LT,
    id_col_data = "Accession",
    id_col_lookup = "uniprot_id",
    universe = universe,
    min_found = 2,
    min_total = 2,
    alpha = 0.05,
    adjust_method = "BH"
  )
  for (complex in names(cx_pruned$keep)) {
    cx_pruned$keep[[complex]] <- resolve_ambiguous_peptides(
      cx_pruned$keep[[complex]],
      lookup_table = EBI_LT
    )
  }
  cx_pruned$keep
}

cat("Finding complexes in log-normalized data...\n")
LogNorm_complexes <- find_and_prune(H9_LogNorm)

cat("Finding complexes in raw data...\n")
Raw_complexes <- find_and_prune(H9_Raw)

################################################################################
## 4. Define complexes to extract                                             ##
################################################################################

# fmt:skip
complexes <- list(
  "B-Wich"                           = "B-WICH chromatin remodelling complex",
  "CRD-mediated"                     = "CRD-mediated mRNA stability complex",
  "Multiaminoacyl-tRNA"              = "Multiaminoacyl-tRNA synthetase complex",
  "Major-Spliceosomal-B"             = "Major Spliceosomal B complex",
  "26S-proteasome"                   = "26S proteasome complex",
  "60S-cytosolic-large-ribosomal-subunit" = "60S cytosolic large ribosomal subunit",
  "Nuclear-pore-complex"             = "Nuclear pore complex",
  "Eukaryotic-translation-initiation-factor-3-complex" = "Eukaryotic translation initiation factor 3 complex",
  "Dynein-1-complex-variant-1"       = "Dynein-1 complex, variant 1",
  "Intraflagellar-transport-complex-B" = "Intraflagellar transport complex B",
  "BSSADCRC"                         = "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  "Dynactin-complex"                 = "Dynactin complex",
  "Laminin-213"                      = "Laminin-213 complex"
)

################################################################################
## 5. Format and save                                                         ##
################################################################################

format_complex <- function(dt) {
  dt <- copy(dt)
  intensity_cols <- grep("^H9_B\\d+_D\\d+$", names(dt), value = TRUE)

  dt <- dt[datatype != "NM"]

  dt[,
    Peptide := paste(
      `Gene name`,
      Accession,
      `Position in master protein`,
      datatype,
      sep = "_"
    )
  ]

  dt <- dt[, c("complex_id", "Peptide", intensity_cols), with = FALSE]

  setnames(dt, intensity_cols, sub("^H9_", "", intensity_cols))
  intensity_cols <- sub("^H9_", "", intensity_cols)

  dt <- dt[!grepl("NOT FOUND", Peptide, fixed = TRUE)]

  Y <- as.matrix(dt[, ..intensity_cols])
  min_non_na <- ceiling(ncol(Y) * 0.4)
  cond_labels <- sub("^B\\d+_(D\\d+)$", "\\1", intensity_cols)
  unique_conds <- unique(cond_labels)
  cond_counts <- vapply(unique_conds, function(cond) {
    cols <- which(cond_labels == cond)
    rowSums(!is.na(Y[, cols, drop = FALSE]))
  }, numeric(nrow(Y)))

  keep <- (rowSums(!is.na(Y)) >= min_non_na) &
    (rowSums(cond_counts >= 2L) >= 2L)
  dt[keep]
}

for (short_name in names(complexes)) {
  full_name <- complexes[[short_name]]
  cat("Processing:", full_name, "\n")

  missing <- character(0)
  if (!full_name %in% names(LogNorm_complexes)) {
    missing <- c(missing, "LogNorm")
  }
  if (!full_name %in% names(Raw_complexes)) {
    missing <- c(missing, "Raw")
  }
  if (length(missing) > 0) {
    warning(full_name, " not found in: ", paste(missing, collapse = ", "))
    next
  }

  lognorm_dt <- format_complex(LogNorm_complexes[[full_name]])
  raw_dt <- format_complex(Raw_complexes[[full_name]])

  fwrite(
    lognorm_dt,
    file.path(output_dir, paste0(short_name, "_WideLogNorm.csv"))
  )
  fwrite(
    raw_dt,
    file.path(output_dir, paste0(short_name, "_WideNoLogNoNorm.csv"))
  )

  cat(
    "  Saved",
    short_name,
    "| LogNorm:",
    nrow(lognorm_dt),
    "rows",
    "| NoLogNoNorm:",
    nrow(raw_dt),
    "rows\n"
  )
}
