# Libraries
library(scico)

# Define ChimeraX path
chimeraX_path <- # Insert path to ChimeraX executable here

# Set working directory to project root
setwd(here::here())

# Define output directory for visualizations
visualization_output_dir <- "ComplexoVisualizer/Results"

if (!dir.exists(visualization_output_dir)) {
  dir.create(visualization_output_dir, recursive = TRUE)
}

# Directory for assembled complexes (AF subunits superimposed via matchmaker)
assembled_complexes_dir <- "ComplexoVisualizer/AssembledComplexes"

if (!dir.exists(assembled_complexes_dir)) {
  dir.create(assembled_complexes_dir, recursive = TRUE)
}

# Directory for raw fetched AlphaFold subunit structures
af_cache_dir <- "ComplexoVisualizer/AfCache"

if (!dir.exists(af_cache_dir)) {
  dir.create(af_cache_dir, recursive = TRUE)
}

# Source functions
source("ComplexoVisualizer/MakeCXC.R")

# Define color palette for types of datatype
datatype_palette <- function() {
  c(
    "NM" = "#00441b",
    "NMpeptide" = "#02d054ff",
    "RmCys" = "#006d2c",
    "Deglyco" = "#99d38d",
    "Phospho" = "#be4e16",
    "FreeCys" = "#ff8000",
    "LysAc" = "#fec389"
  )
}

# Define protein chain color palette
protein_chain_palette <- function(n) {
  scico(n, palette = "batlow", end = 0.9)
}

# Define paths for each complex.
# chain_map: named character vector, names = chain IDs in the structure file,
# values = UniProt accessions as they appear in the cluster_summary Accession column.
Cav1.2 <- list(
  name = "Cav1.2 voltage-gated calcium channel complex, CACNA2D1-CACNB2 variant",
  cluster_summary = "ComplexoFinder/10-Results/Cav1.2_voltage-gated_calcium_channel_complex,_CACNA2D1-CACNB2_variant/cluster_summary.csv",
  structure_path = "ComplexoVisualizer/Structures/fold_cav12_voltage_gaet_calcium_channel_cacna2d1_cacnb2_model_0.cif",
  # AlphaFold3 multimer: chains assigned by sequence length (A=660aa, B=1103aa, C=2221aa)
  #fmt: skip
  chain_map = c(
    "A" = "Q08289", # CACNB2   (beta-2 subunit)
    "B" = "P54289", # CACNA2D1 (alpha-2/delta-1 subunit)
    "C" = "Q13936" # CACNA1C  (Cav1.2 pore-forming alpha-1 subunit)
  )
)
SNARE <- list(
  name = "SNARE complex STX4-SNAP29-SEC22b",
  cluster_summary = "ComplexoFinder/10-Results/SNARE_complex_STX4-SNAP29-SEC22b/cluster_summary.csv",
  structure_path = "ComplexoVisualizer/Structures/fold_snare_complex_stx4_snap29_sec22b_model_0.cif",
  # AlphaFold3 multimer: chains assigned by sequence length (A=297aa, B=258aa, C=215aa)
  #fmt: skip
  chain_map = c(
    "A" = "Q12846", # STX4   (syntaxin-4)
    "B" = "O95721", # SNAP29 (synaptosomal-associated protein 29)
    "C" = "O75396" # SEC22B (vesicle-trafficking protein SEC22b)
  )
)
Proteasome_26S <- list(
  name = "26S proteasome complex",
  cluster_summary = "ComplexoFinder/10-Results/26S_proteasome_complex/cluster_summary.csv",
  # PDB 5L4G: one 20S core particle + one 19S RP base (AAA+ ring, PSMC1-6).
  # PSMD1-14 and ADRM1 from the cluster summary are NOT present in this structure;
  # peptides from those subunits will be silently skipped by make_cxc.
  structure_path = "ComplexoVisualizer/Structures/5l4g.cif",
  #fmt: skip
  chain_map = c(
    # 20S alpha ring – copy 1 / copy 2 (order matches entity order in struct_ref)
    "A" = "P60900", "N" = "P60900",  # PSMA6
    "B" = "P25787", "O" = "P25787",  # PSMA2
    "C" = "P25789", "P" = "P25789",  # PSMA4
    "D" = "O14818", "Q" = "O14818",  # PSMA7
    "E" = "P28066", "R" = "P28066",  # PSMA5
    "F" = "P25786", "S" = "P25786",  # PSMA1
    "G" = "P25788", "T" = "P25788",  # PSMA3
    # 20S beta ring – copy 1 / copy 2
    "1" = "P20618", "U" = "P20618",  # PSMB1
    "2" = "P49721", "V" = "P49721",  # PSMB2
    "3" = "P49720", "W" = "P49720",  # PSMB3
    "4" = "P28070", "X" = "P28070",  # PSMB4
    "5" = "P28074", "Y" = "P28074",  # PSMB5
    "6" = "P28072", "Z" = "P28072",  # PSMB6
    "7" = "Q99436", "8" = "Q99436",  # PSMB7
    # 19S regulatory particle base (AAA+ ring) – single copy
    "H" = "P35998",  # PSMC2
    "I" = "P62191",  # PSMC1
    "J" = "P62195",  # PSMC5
    "K" = "P43686",  # PSMC4
    "L" = "P62333",  # PSMC6
    "M" = "P17980"   # PSMC3
  )
)
Calcineurin_Calmodulin_AKAP5 <- list(
  name = "Calcineurin-Calmodulin-AKAP5 complex, gamma-R1 variant",
  cluster_summary = "ComplexoFinder/10-Results/Calcineurin-Calmodulin-AKAP5_complex,_gamma-R1_variant/cluster_summary.csv",
  structure_path = "ComplexoVisualizer/Structures/fold_calcineurin_calmodulin_akap5_model_0.cif",
  # AlphaFold3 multimer: chains assigned by sequence length (A=170aa, B=512aa, C=427aa, D=149aa)
  #fmt: skip
  chain_map = c(
    "A" = "P63098", # PPP3R1 (calcineurin regulatory subunit B, type 1)
    "B" = "P48454", # PPP3CC (calcineurin catalytic subunit gamma)
    "C" = "P24588", # AKAP5  (A-kinase anchor protein 5)
    "D" = "P0DP23" # CALM1  (calmodulin-1)
  )
)
Dynein_1_var_1 <- list(
  name = "Dynein-1 complex, variant 4",
  cluster_summary = "ComplexoFinder/10-Results/Dynein-1_complex,_variant_4/cluster_summary.csv",
  # PDB 9DGP: motor domain of DYNC1H1 only (pig homolog, single chain A).
  # DYNC1I1, DYNC1LI1, DYNLL1, DYNLRB1, DYNLT1 are absent from this structure;
  # peptides from those subunits will be silently skipped by make_cxc.
  structure_path = "ComplexoVisualizer/Structures/9DGP.cif",
  chain_map = c(
    "A" = "Q14204" # DYNC1H1 (dynein cytoplasmic 1 heavy chain 1)
  )
)

# No combined structure exists for the two complexes below.
# make_cxc_af_subunits() fetches individual AlphaFold models per subunit and,
# when reference_pdb is provided, uses matchmaker to dock them into the
# reference coordinate frame before coloring.

SWI_SNF <- list(
  name = "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  cluster_summary = "ComplexoFinder/10-Results/Brain-specific_SWI_SNF_ATP-dependent_chromatin_remodeling_complex,_ARID1A-SMARCA2_variant/cluster_summary.csv",
  structure_path = "",
  chain_map = c(),
  subunit_af_accessions = c(
    "ACTB" = "P60709",
    "ACTL6B" = "O94805",
    "ARID1A" = "O14497",
    "SMARCA2" = "P51531",
    "SMARCB1" = "Q12824",
    "SMARCC2" = "Q8TAQ2",
    "SMARCD2" = "Q92925",
    "SMARCE1" = "Q969G3"
  ),
  # PDB 6LTJ: cryo-EM of the human canonical BAF complex.
  # 5 exact matches + 3 paralogs in the same structural position:
  #   SMARCA4 (chain I) stands in for SMARCA2 (~75% identity)
  #   ACTL6A  (chain J) stands in for ACTL6B  (~80% identity)
  #   SMARCD1 (chain P) stands in for SMARCD2 (~68% identity)
  # matchmaker uses sequence alignment, so all 8 subunits dock correctly.
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
  name = "Eukaryotic translation initiation factor 3 complex",
  cluster_summary = "ComplexoFinder/10-Results/Eukaryotic_translation_initiation_factor_3_complex/cluster_summary.csv",
  structure_path = "",
  chain_map = c(),
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
  # PDB 6ZMW: cryo-EM of human eIF3 on a 48S translation initiation complex.
  # All 13 subunits are present; every AlphaFold model can be positioned.
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

Neuron_SWI_SNF <- list(
  name = "Neuron-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  cluster_summary = "ComplexoFinder/10-Results/Neuron-specific_SWI_SNF_ATP-dependent_chromatin_remodeling_complex,_ARID1A-SMARCA2_variant/cluster_summary.csv",
  structure_path = "",
  chain_map = c(),
  subunit_af_accessions = c(
    "ACTB" = "P60709",
    "ACTL6B" = "O94805",
    "ARID1A" = "O14497",
    "DPF1" = "Q92782",
    "DPF3" = "Q92784",
    "SMARCA2" = "P51531",
    "SMARCB1" = "Q12824",
    "SMARCC1" = "Q92922",
    "SMARCC2" = "Q8TAQ2",
    "SMARCD1" = "Q96GM5",
    "SMARCD3" = "Q6STE5",
    "SMARCE1" = "Q969G3"
  ),
  # PDB 9A0K: cryo-EM of human nBAF complex.
  # 7 exact matches + 3 paralogs in equivalent structural positions:
  #   SMARCA4 (chain 4) stands in for SMARCA2 (~75% identity)
  #   ACTL6A  (chain 6) stands in for ACTL6B  (~80% identity)
  #   DPF2    (chain P) stands in for DPF1 and DPF3 (~55% identity)
  # SMARCD3 has no reference chain and will be left spatially unpositioned.
  reference_pdb = "9A0K",
  #fmt: skip
  reference_chain_map = c(
    "P60709" = "B",  # ACTB    (exact)
    "O14497" = "A",  # ARID1A  (exact)
    "Q12824" = "b",  # SMARCB1 (exact)
    "Q92922" = "1",  # SMARCC1 (exact)
    "Q8TAQ2" = "2",  # SMARCC2 (exact)
    "Q96GM5" = "D",  # SMARCD1 (exact)
    "Q969G3" = "E",  # SMARCE1 (exact)
    "P51531" = "4",  # SMARCA2 -> SMARCA4 chain
    "O94805" = "6",  # ACTL6B  -> ACTL6A  chain
    "Q92782" = "P",  # DPF1    -> DPF2    chain
    "Q92784" = "P"   # DPF3    -> DPF2    chain
  )
)

################################################################################
## Run #########################################################################
################################################################################

# Complexes with a combined structure file -> make_cxc()
structure_complexes <- list(
  Cav1.2,
  SNARE,
  Proteasome_26S,
  Calcineurin_Calmodulin_AKAP5,
  Dynein_1_var_1
)

for (cx in structure_complexes) {
  folder_name <- gsub("[^A-Za-z0-9]+", "_", cx$name)
  cxc_dir <- file.path(visualization_output_dir, folder_name, "cxcFiles")
  dist_dir <- file.path(visualization_output_dir, folder_name, "distanceFiles")

  make_cxc(
    complex_obj = cx,
    output_dir = cxc_dir,
    membership_threshold = 0.5,
    chain_lighten_factor = 0.6,
    save_image = FALSE
  )

  make_cxc_cluster_combined(
    complex_obj = cx,
    output_dir = cxc_dir,
    membership_threshold = 0.5,
    save_image = FALSE
  )

  compute_ptm_distances(
    complex_obj = cx,
    output_dir = dist_dir,
    membership_threshold = 0.5
  )
}

# Complexes without a combined structure -> make_cxc_af_subunits()
# Subunits are fetched live from AlphaFold when the .cxc script runs in ChimeraX.
af_subunit_complexes <- list(
  SWI_SNF,
  eIF3,
  Neuron_SWI_SNF
)

for (cx in af_subunit_complexes) {
  folder_name <- gsub("[^A-Za-z0-9]+", "_", cx$name)
  cxc_dir <- file.path(visualization_output_dir, folder_name, "cxcFiles")
  dist_dir <- file.path(visualization_output_dir, folder_name, "distanceFiles")

  make_cxc_af_subunits(
    complex_obj = cx,
    output_dir = cxc_dir,
    membership_threshold = "adaptive",
    chain_lighten_factor = 0.6,
    save_image = FALSE
  )

  make_cxc_af_subunits_cluster_combined(
    complex_obj = cx,
    output_dir = cxc_dir,
    membership_threshold = "adaptive",
    save_image = FALSE
  )

  compute_ptm_distances_af_subunits(
    complex_obj = cx,
    output_dir = dist_dir,
    membership_threshold = "adaptive",
    cache_dir = assembled_complexes_dir,
    af_cache_dir = af_cache_dir
  )
}
