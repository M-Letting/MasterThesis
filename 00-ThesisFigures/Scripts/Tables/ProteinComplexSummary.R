# Libraries
library(data.table)
library(gt)

# Set working directory to folder
setwd(here::here("ComplexoFinder"))

# Check if data files exist, if not run data loading and handling script
if (
  !file.exists("00-Data/H9_data.rds") ||
    !file.exists("00-Data/IMR90_data.rds")
) {
  source("01-DataProcessing/ProcessData.R")
}

# Check if required lookup tables exist, if not run script to create them
if (
  !file.exists("00-LookupTables/ebi_cp_lookup.rds") ||
    !file.exists("00-LookupTables/human_proteome.rds")
) {
  source("00-LookupTables/01-GetLookupTables.R")
}

# Check if output directory exists, if not create it
if (!dir.exists("10-Results")) {
  dir.create("10-Results")
}

# Source functions
source("02-Identification/AssignComplexes.R")
source("02-Identification/IdentifyComplexes.R")
source("03-Constraints/Constraints.R")
source("04-Clustering/VsClust.R")
source("04-Clustering/VsClustWrappers.R")
source("04-Clustering/ClusterSummary.R")
source("05-Quantification/Quantification.R")
source("09-Visualization/ExpressionMembershipPlot.R")
source("09-Visualization/SequencePlot.R")
source("09-Visualization/QuantificationPlot.R")

################################################################################
## Options #####################################################################
################################################################################

# Choose which dataset to analyze
dataset <- "H9" # Options: "H9", "IMR90"

# Verbose output
verbose <- TRUE # TRUE for detailed output, FALSE for minimal output

# Include NM (protein level abundance) in analysis
include_NM <- FALSE # TRUE to include protein-level data, FALSE to exclude

# Include constrainsts in clustering
include_constraints <- TRUE # TRUE to include constraints, FALSE to exclude

# Specify complexes to analyze - for all options, see EBI complex lookup table
complexes <- c(
  "B-WICH chromatin remodelling complex",
  "CRD-mediated mRNA stability complex",
  "Multiaminoacyl-tRNA synthetase complex",
  "Major Spliceosomal B complex",
  "26S proteasome complex",
  "60S cytosolic large ribosomal subunit",
  "Nuclear pore complex",
  "Eukaryotic translation initiation factor 3 complex",
  "Dynein-1 complex, variant 1",
  "Intraflagellar transport complex B",
  "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  "Neuron-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  "Calcineurin-Calmodulin-AKAP5 complex, gamma-R1 variant",
  "SNARE complex STX4-SNAP29-SEC22b",
  "Cortical microtubule stabilization complex, KANK1 variant",
  "LIFT actin modulation complex",
  "Dynein-1 complex, variant 4",
  "WASH complex, variant WASHC1/WASHC2C",
  "Cav1.2 voltage-gated calcium channel complex, CACNA2D1-CACNB3 variant",
  "Dynactin complex",
  "Laminin-213 complex"
)

################################################################################
## Run ComplexoFinder ##########################################################
################################################################################

#############
# Load data #
#############

if (dataset == "H9") {
  data_list <- readRDS("00-Data/H9_data.rds")
} else if (dataset == "IMR90") {
  data_list <- readRDS("00-Data/IMR90_data.rds")
} else {
  stop("Invalid dataset choice. Please choose 'H9' or 'IMR90'.")
}

# Define intensity column patterns for each dataset (after loading data)
if (dataset == "H9") {
  intensity_cols <- grep(
    "^H9_B[0-9]+_D[0-9]+$",
    names(data_list[[1]]),
    value = TRUE
  )
} else if (dataset == "IMR90") {
  intensity_cols <- grep("^IMR90", names(data_list[[1]]), value = TRUE)
} else {
  stop("Invalid dataset option. Choose 'H9' or 'IMR90'.")
}

n_rep <- 3
n_cond <- length(intensity_cols) / n_rep
if (n_cond != as.integer(n_cond)) {
  stop("Cannot infer n_cond: length(intensity_cols) is not divisible by n_rep")
}
n_cond <- as.integer(n_cond)

cond_regex <- if (dataset == "H9") {
  "^H9_B[0-9]+_(D[0-9]+)$"
} else {
  "^IMR90_B[0-9]+_(D[0-9]+)$"
}

###############################################################
# Find complexes using EBI annotated list and prune complexes #
###############################################################

# Lookup table for EBI complexes
EBI_LT <- readRDS("00-LookupTables/ebi_cp_lookup.rds")

# Find complexes using EBI annotated list
complexes_EBI <- find_complexes(
  data_list = data_list,
  lookup_table = EBI_LT,
  id_col_data = "Accession",
  id_col_lookup = "uniprot_id"
)

# Build universe from all proteins identified across datasets (canonical accessions)
universe <- unique(sub(
  "-1$",
  "",
  trimws(unlist(lapply(data_list, function(dt) {
    unlist(lapply(dt[["Accession"]], function(x) {
      if (is.na(x) || x == "") {
        return(character(0))
      }
      trimws(strsplit(x, ";")[[1]])
    }))
  })))
))
universe <- universe[nzchar(universe) & !is.na(universe)]

# Prune complexes using statistical test
complexes_EBI <- prune_complexes_stat(
  all_complexes = complexes_EBI$complexes,
  lookup_table = EBI_LT,
  id_col_data = "Accession",
  id_col_lookup = "uniprot_id",
  universe = universe,
  min_found = 3,
  min_total = 3,
  alpha = 0.05,
  adjust_method = "BH"
)

# Cleanup
rm(universe)

###############################
# Remove NM data if requested #
###############################
if (!include_NM) {
  for (complex in names(complexes_EBI$keep)) {
    data <- complexes_EBI$keep[[complex]]
    if (!is.null(data)) {
      data <- data[datatype != "NM"]
      complexes_EBI$keep[[complex]] <- data
    }
  }
}

#####################################################
# Check for duplicate complexes in complexes vector #
#####################################################
if (length(complexes) != length(unique(complexes))) {
  stop(
    "Duplicate complex names found in 'complexes' vector. Please ensure all complex names are unique."
  )
}

##########################################################################
# Check that all specified complexes are present in the pruned complexes #
##########################################################################
missing_complexes <- setdiff(complexes, names(complexes_EBI$keep))
if (length(missing_complexes) > 0) {
  warning(
    paste(
      "The following complexes were specified but not found in the pruned complexes:",
      paste(missing_complexes, collapse = ", ")
    )
  )
}


##########################################################################
# Resolve ambiguous peptides: if a peptide maps to multiple proteins,    #
# keep it only when exactly one of those proteins is a complex member,   #
# updating the accession to that member. Remove otherwise.               #
##########################################################################
for (complex in names(complexes_EBI$keep)) {
  data <- complexes_EBI$keep[[complex]]

  data_clean <- resolve_ambiguous_peptides(
    data,
    accession_col = "Accession",
    position_col = "Position in master protein",
    lookup_table = EBI_LT,
    id_col_lookup = "uniprot_id"
  )

  complexes_EBI$keep[[complex]] <- data_clean
}


##########################################################################
# Create summary table for each complex, including name,  number of      #
# proteins, features per view, and total peptides (NMpeptide + all PTMs) #
##########################################################################
summary_table <- data.table(
  Complex = character(),
  NumProteins = integer(),
  NumAcLys = integer(),
  NumDeglyco = integer(),
  NumFreeCys = integer(),
  NumRmCys = integer(),
  NumPhospho = integer(),
  NumPeptides = integer(),
  TotalFeatures = integer()
)

for (complex in complexes) {
  data <- complexes_EBI$keep[[complex]]
  if (is.null(data)) {
    next
  }

  num_proteins <- length(unique(unlist(strsplit(data$Accession, ";"))))
  num_ac_lys <- nrow(data[datatype == "LysAC"])
  num_deglyco <- nrow(data[datatype == "Deglyco"])
  num_free_cys <- nrow(data[datatype == "FreeCys"])
  num_rm_cys <- nrow(data[datatype == "RmCys"])
  num_phospho <- nrow(data[datatype == "Phospho"])
  num_peptides <- nrow(data[datatype == "NMpeptide"])
  total_features <- nrow(data)

  summary_table <- rbind(
    summary_table,
    data.table(
      Complex = complex,
      NumProteins = num_proteins,
      NumAcLys = num_ac_lys,
      NumDeglyco = num_deglyco,
      NumFreeCys = num_free_cys,
      NumRmCys = num_rm_cys,
      NumPhospho = num_phospho,
      NumPeptides = num_peptides,
      TotalFeatures = total_features
    )
  )
}

# Make nice table with gt
gt_table <- gt(summary_table) |>
  cols_label(
    Complex = "Complex",
    NumProteins = "# Proteins",
    NumAcLys = "# LysAc",
    NumDeglyco = "# Deglyco",
    NumFreeCys = "# FreeCys",
    NumRmCys = "# RmCys",
    NumPhospho = "# Phospho",
    NumPeptides = "# NM Peptides",
    TotalFeatures = "Total Features"
  ) |>
  fmt_number(columns = 2:9, use_seps = TRUE, decimals = 0) |>
  cols_align(align = "right", columns = 2:9) |>
  cols_align(align = "left", columns = 1) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#cccccc", weight = px(1)),
    locations = list(
      cells_body(columns = 2:9),
      cells_column_labels(columns = 2:9)
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
setwd(here::here())
gtsave(
  gt_table,
  "00-ThesisFigures/Tables/ProteinComplexSummary.png",
  vwidth = 2900,
  vheight = 2100
)

