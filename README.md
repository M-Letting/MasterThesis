# Master's Thesis — Elucidation of the behavior of protein complexes and their PTMs in developing organoids using new statistical methods


**Author:** Malthe Letting

This repository contains the code developed for my Master's thesis. The central contribution is **ComplexoFinder**, a pipeline for identifying and characterising protein complexoforms from bottom-up LC-MS/MS data, together with benchmarks validating each methodological step and tools for visualising the results.

---

## Repository Structure

### `ComplexoFinder/`
The main analysis pipeline. Starting point is `RunComplexoFinder.R`.

| Folder | Contents |
|---|---|
| `00-Data/` | Input data loading and storage |
| `00-LookupTables/` | Reference tables (EBI complex portal, human proteome) |
| `01-DataProcessing/` | Pre-processing of data |
| `02-Identification/` | Complex identification and assignment logic |
| `03-Constraints/` | Sequence-derived constraint generation |
| `04-Clustering/` | VSClust clustering with constraints |
| `05-Quantification/` | Quantification of complex abundance |
| `09-Visualization/` | In-pipeline plotting helpers |
| `10-Results/` | Pipeline output (generated at runtime) |

### `Q-Benchmark/`
Systematic benchmarks evaluating each major step of ComplexoFinder.

| Folder | Type of Benchmark |
|---|---|
| `01-FindDCF/` | Benchmarks based on SWATH-MS data (similar to benchmark used by COPF and ProteoForge) |
| `02-DCF/` | Benchmarks using ProteoMaker in-silico generated datasets |
| `03-Clustering/` | Constrained clustering vs. standard clustering |
| `04-Structures/` | Structural assessment of identified complexoforms |

### `ComplexoVisualiserThesis/`
Functions and more used for creating protein complex structures and highlighting modified peptides based on results of ComplexoFinder. Generates the structures and distance plots used in final version of the thesis.

### `ComplexoVisualizer/`
Functions and more used for creating protein complex structures and highlighting modified peptides based on results of ComplexoFinder and ProteoForge (used for "Q-Benchmark/04-Structures").

### `00-ThesisFigures/`
Scripts (`Scripts/`) and outputs (`Figures/`, `Tables/`) for most figures, plots, and tables used during my thesis.

---

## Getting Started

1. Clone the repository.
2. Install R dependencies by running the install script from the repository root:

```r
source("install_packages.R")
```

This installs all required CRAN, Bioconductor, and GitHub packages. See [install_packages.R](install_packages.R) for the full list.

3. Install Python dependencies (required for `Q-Benchmark/` scripts that use ProteoForge):

```bash
pip install -r requirements.txt
```

The ProteoForge package itself is downloaded automatically from GitHub the first time the relevant scripts are run. See [requirements.txt](requirements.txt) for details.

4. Run the main pipeline:

```r
source("ComplexoFinder/RunComplexoFinder.R")
```

Benchmarks in `Q-Benchmark/` can be run independently; each subfolder contains numbered scripts that should be executed in order. These can also be sourced through the "10-Run" scripts in each subfolder.

### External data

Large data files are excluded from this repository. Before running `Q-Benchmark/01-FindDCF/`, download the required SWATH-MS input file from:

> **[https://doi.org/10.5281/zenodo.17795845](https://doi.org/10.5281/zenodo.17795845)** (originally from the [ProteoForge GitHub](https://github.com/LangeLab/ProteoForge_Analysis))

Place the downloaded file at:

```
Q-Benchmark/01-FindDCF/data/site02_global_q_0.01_applied_to_local_global.txt
```

The cerebral organoid (CBO) proteomics data used as input in `ComplexoFinder/` and `Q-Benchmark/` is not included in this repository, as the associated article is currently under review. The data will be added once the article is published or permission has been granted by the authors.

---

