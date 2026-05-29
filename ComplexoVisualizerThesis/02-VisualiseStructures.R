# Libraries
library(scico)

# Define ChimeraX path
chimeraX_path <- # Insert path to ChimeraX executable here

# Set working directory to project root
setwd(here::here())

# Output directories
visualization_output_dir <- "ComplexoVisualiserThesis/00-Data/Results"
assembled_complexes_dir  <- "ComplexoVisualiserThesis/00-Data/AssembledComplexes"
af_cache_dir             <- "ComplexoVisualiserThesis/00-Data/AfCache"

for (d in c(visualization_output_dir, assembled_complexes_dir, af_cache_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Source functions
source("ComplexoVisualizer/MakeCXC.R")

# Shared color palettes (used by MakeCXC.R helpers)
datatype_palette <- function() {
  c(
    "NM"        = "#00441b",
    "NMpeptide" = "#02d054ff",
    "RmCys"     = "#006d2c",
    "Deglyco"   = "#99d38d",
    "Phospho"   = "#be4e16",
    "FreeCys"   = "#ff8000",
    "LysAc"     = "#fec389"
  )
}

protein_chain_palette <- function(n) {
  scico(n, palette = "batlow", end = 0.9)
}

################################################################################
## Complex definitions #########################################################
################################################################################

summaries_dir <- "ComplexoVisualiserThesis/00-Data/Summaries"

Proteasome_26S <- list(
  name            = "26S proteasome complex",
  cluster_summary = file.path(summaries_dir, "26S-proteasome_ComplexoFinderSummary.csv"),
  structure_path  = "ComplexoVisualizer/Structures/5l4g.cif",
  #fmt: skip
  chain_map = c(
    "A" = "P60900", "N" = "P60900",  # PSMA6
    "B" = "P25787", "O" = "P25787",  # PSMA2
    "C" = "P25789", "P" = "P25789",  # PSMA4
    "D" = "O14818", "Q" = "O14818",  # PSMA7
    "E" = "P28066", "R" = "P28066",  # PSMA5
    "F" = "P25786", "S" = "P25786",  # PSMA1
    "G" = "P25788", "T" = "P25788",  # PSMA3
    "1" = "P20618", "U" = "P20618",  # PSMB1
    "2" = "P49721", "V" = "P49721",  # PSMB2
    "3" = "P49720", "W" = "P49720",  # PSMB3
    "4" = "P28070", "X" = "P28070",  # PSMB4
    "5" = "P28074", "Y" = "P28074",  # PSMB5
    "6" = "P28072", "Z" = "P28072",  # PSMB6
    "7" = "Q99436", "8" = "Q99436",  # PSMB7
    "H" = "P35998",                  # PSMC2
    "I" = "P62191",                  # PSMC1
    "J" = "P62195",                  # PSMC5
    "K" = "P43686",                  # PSMC4
    "L" = "P62333",                  # PSMC6
    "M" = "P17980"                   # PSMC3
  )
)

Dynein_1_var_1 <- list(
  name            = "Dynein-1 complex, variant 1",
  cluster_summary = file.path(summaries_dir, "Dynein-1-complex-variant-1_ComplexoFinderSummary.csv"),
  structure_path  = "ComplexoVisualizer/Structures/9e13.cif",
  #fmt: skip
  chain_map = c(
    "A" = "Q14204", "B" = "Q14204",  # DYNC1H1
    "C" = "Q13409", "D" = "Q13409",  # DYNC1I2
    "E" = "O43237", "F" = "O43237"   # DYNC1LI2
  )
)

Brain_SWI_SNF <- list(
  name            = "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  cluster_summary = file.path(summaries_dir, "Brain-SWI-SNF_ComplexoFinderSummary.csv"),
  structure_path  = "",
  chain_map       = c(),
  subunit_af_accessions = c(
    "ACTB"    = "P60709",
    "ACTL6B"  = "O94805",
    "ARID1A"  = "O14497",
    "SMARCA2" = "P51531",
    "SMARCB1" = "Q12824",
    "SMARCC2" = "Q8TAQ2",
    "SMARCD2" = "Q92925",
    "SMARCE1" = "Q969G3"
  ),
  reference_pdb = "6LTJ",
  #fmt: skip
  reference_chain_map = c(
    "P60709" = "K",  # ACTB    (exact)
    "O14497" = "L",  # ARID1A  (exact)
    "Q12824" = "M",  # SMARCB1 (exact)
    "Q8TAQ2" = "N",  # SMARCC2 (exact)
    "Q969G3" = "Q",  # SMARCE1 (exact)
    "P51531" = "I",  # SMARCA2 -> SMARCA4 chain
    "O94805" = "J",  # ACTL6B  -> ACTL6A  chain
    "Q92925" = "P"   # SMARCD2 -> SMARCD1 chain
  )
)

eIF3 <- list(
  name            = "Eukaryotic translation initiation factor 3 complex",
  cluster_summary = file.path(summaries_dir, "Eukaryotic-translation-initiation-factor-3-complex_ComplexoFinderSummary.csv"),
  structure_path  = "",
  chain_map       = c(),
  subunit_af_accessions = c(
    "EIF3A" = "Q14152",
    "EIF3B" = "P55884",
    "EIF3C" = "Q99613",
    "EIF3D" = "O15371",
    "EIF3E" = "P60228",
    "EIF3F" = "O00303",
    "EIF3G" = "O75821",
    "EIF3H" = "O15372",
    "EIF3I" = "Q13347",
    "EIF3J" = "O75822",
    "EIF3K" = "Q9UBQ5",
    "EIF3L" = "Q9Y262",
    "EIF3M" = "Q7L2H7"
  ),
  reference_pdb = "6ZMW",
  #fmt: skip
  reference_chain_map = c(
    "Q14152" = "u",  # EIF3A
    "P55884" = "1",  # EIF3B
    "Q99613" = "y",  # EIF3C
    "O15371" = "x",  # EIF3D
    "P60228" = "v",  # EIF3E
    "O00303" = "4",  # EIF3F
    "O75821" = "o",  # EIF3G
    "O15372" = "8",  # EIF3H
    "Q13347" = "2",  # EIF3I
    "O75822" = "z",  # EIF3J
    "Q9UBQ5" = "3",  # EIF3K
    "Q9Y262" = "5",  # EIF3L
    "Q7L2H7" = "6"   # EIF3M
  )
)

################################################################################
## Run #########################################################################
################################################################################

# Complexes with a pre-assembled structure file
structure_complexes <- list(
  Proteasome_26S,
  Dynein_1_var_1
)

for (cx in structure_complexes) {
  folder_name <- gsub("[^A-Za-z0-9]+", "_", cx$name)
  cxc_dir  <- file.path(visualization_output_dir, folder_name, "cxcFiles")
  dist_dir <- file.path(visualization_output_dir, folder_name, "distanceFiles")

  make_cxc(
    complex_obj          = cx,
    output_dir           = cxc_dir,
    membership_threshold = 0.5,
    chain_lighten_factor = 0.6,
    save_image           = FALSE
  )

  make_cxc_cluster_combined(
    complex_obj          = cx,
    output_dir           = cxc_dir,
    membership_threshold = 0.5,
    save_image           = FALSE
  )

  compute_ptm_distances(
    complex_obj          = cx,
    output_dir           = dist_dir,
    membership_threshold = 0.5
  )

  compute_ptm_distances_combined(
    complex_obj = cx,
    output_dir  = dist_dir
  )
}

# Complexes built from individual AlphaFold subunit models
af_subunit_complexes <- list(
  Brain_SWI_SNF,
  eIF3
)

for (cx in af_subunit_complexes) {
  folder_name <- gsub("[^A-Za-z0-9]+", "_", cx$name)
  cxc_dir  <- file.path(visualization_output_dir, folder_name, "cxcFiles")
  dist_dir <- file.path(visualization_output_dir, folder_name, "distanceFiles")

  make_cxc_af_subunits(
    complex_obj          = cx,
    output_dir           = cxc_dir,
    membership_threshold = 0.5,
    chain_lighten_factor = 0.6,
    save_image           = FALSE
  )

  make_cxc_af_subunits_cluster_combined(
    complex_obj          = cx,
    output_dir           = cxc_dir,
    membership_threshold = 0.5,
    save_image           = FALSE
  )

  # Per-cluster distances must run first to populate the AF coordinate cache
  compute_ptm_distances_af_subunits(
    complex_obj          = cx,
    output_dir           = dist_dir,
    membership_threshold = 0.5,
    cache_dir            = assembled_complexes_dir,
    af_cache_dir         = af_cache_dir
  )

  compute_ptm_distances_af_subunits_combined(
    complex_obj  = cx,
    output_dir   = dist_dir,
    cache_dir    = assembled_complexes_dir,
    af_cache_dir = af_cache_dir
  )
}
