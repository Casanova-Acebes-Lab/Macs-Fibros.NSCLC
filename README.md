# Macs-Fibros.NSCLC

Reproducible single-cell RNA-seq analysis of the lung tumour microenvironment developed in the Casanova-Acebes Lab.

## Overview

This repository contains the R workflow used to reconstruct, annotate and analyze the single-cell RNA-seq dataset used in this study.

The experimental design includes three mouse backgrounds:

- `WT`
- `CD169DTR`
- `CCR2KO`

with naive and tumour-bearing conditions.

The analysis focuses on stromal and immune populations of the lung tumour microenvironment, with particular emphasis on fibroblast heterogeneity and macrophage-fibroblast interactions. The workflow includes global quality control and clustering, compartment-specific subclustering and annotation, manuscript figure generation, trajectory analysis and external validation.

## Repository structure

```text
Macs-Fibros.NSCLC/
├── R/                         Analysis and figure-generation scripts
├── data/
│   ├── raw/                   Primary sample objects (not tracked by Git)
│   ├── barcodes.doublets/     Doublet barcode lists
│   ├── HVG/                   Historical highly variable genes
│   ├── signatures/            Cell-annotation signatures
│   ├── Figure4K/              Figure-specific inputs
│   ├── Figure6/               Figure-specific inputs
│   └── SuppFigure2/           External-validation documentation
├── env/                       Conda/R environment snapshots
├── results/                   Generated objects, tables and plots
├── config.tsv
└── README.md
```

All scripts should be run from the repository root.

## Software environment

The final analysis environment uses R 4.4.3 and contains the packages required for both the main Seurat workflow and the Monocle3 trajectory analysis.

Create the environment with:

```bash
conda env create -f env/environment.yml
```

Activate it with:

```bash
conda activate Macs-Fibros.NSCLC
```

Additional package and environment snapshots are stored in `env/`, including:

```text
env/R-packages.tsv
env/conda-packages.tsv
env/conda-linux-64.lock.txt
env/conda-history.txt
```

## Configuration

Paths are defined in `config.tsv`.

The default configuration is:

```text
key      value
rsdir    results
datadir  data
outdir   results/plots
```

Analysis objects are stored under `results/objects/`, plots under `results/plots/`, and source-data tables under `results/tables/`.

## Input data

### Primary study data

Large sample-level R objects are not stored in this GitHub repository.

The expected filenames and directory structure are documented in:

```text
data/raw/README.md
```

The primary data will be deposited in an appropriate public repository. Accession information will be added when available.

### Doublet barcodes

Previously identified doublet barcodes required during QC are versioned in:

```text
data/barcodes.doublets/
```

The repository contains one small barcode object for each of the 12 samples.

### Historical highly variable genes

The historical set of 2,000 highly variable genes used to reproduce the original global analysis is stored in:

```text
data/HVG/fibroblasts_historical_HVGs.rds
```

PCA and UMAP coordinates are recalculated by the workflow and are not stored as input files.

### Annotation signatures

Curated signatures used during cell-type annotation are stored in:

```text
data/signatures/
```

### External validation datasets

Previously published datasets used for Supplementary Figure 2D-F are not redistributed in this repository.

Their provenance, access information and expected local file structure are documented in:

```text
data/SuppFigure2/README.md
```

## Running the main analysis

Run the analysis from the repository root.

### 1. Quality control and global clustering

```bash
Rscript R/01_QC.R
```

This script loads the 12 sample-level objects, removes previously identified doublets, applies cell-level QC, performs normalization and scaling, uses the historical 2,000-gene HVG set, and recalculates PCA, nearest-neighbour graph, clustering and UMAP.

The reconstructed QC workflow was validated against the historical analysis and reproduced 51,195 cells and the original 26 global clusters.

Main output:

```text
results/objects/fibroblasts.QC.rds
```

### 2. Macrophage, monocyte and dendritic-cell subclustering

```bash
Rscript R/02_Subcluster.Macros.R
```

This script performs dedicated subclustering and annotation of the macrophage, monocyte and dendritic-cell compartments.

Main output:

```text
results/objects/fibroblasts.Macrophage_DC.annotated.rds
```

### 3. Fibroblast and mesothelial subclustering

```bash
Rscript R/03_Subcluster.Fibros.R
```

This step reconstructs the fibroblast/mesothelial compartment and defines the stromal populations used throughout the study:

```text
Fibroblasts Universal
Fibroblasts Activated 1
Fibroblasts Activated 2
Mesothelial
```

Main output:

```text
results/objects/fibroblasts.Fibros.annotated.rds
```

### 4. CD45-negative annotation

```bash
Rscript R/04_CD45neg.R
```

This script combines the fibroblast/mesothelial annotation with the remaining CD45-negative populations and annotates endothelial, epithelial, mesenchymal and tumour-associated populations.

Main output:

```text
results/objects/fibroblasts.CD45NP.annotated.rds
```

### 5. Cancer-cell analysis

```bash
Rscript R/05_CancerCells.R
```

This is a dedicated analysis of tumour epithelial cells. It uses the QC object together with the fibroblast/mesothelial annotation.

### 6. T-cell analysis

```bash
Rscript R/06_Tcells.R
```

This script performs selection, reclustering and annotation of the T-cell compartment.

Main output:

```text
results/objects/fibroblasts.CD45.Tcells.annotated.rds
```

### 7. Final annotation merge

```bash
Rscript R/07_Join.Clustering.R
```

This script combines the compartment-specific annotations into the final annotated single-cell object.

Main downstream object:

```text
results/objects/fibroblasts.Not.Integrated.annotated.rds
```

This object is used by most downstream figure-generation scripts.

## Figure generation

After the main annotation workflow has been completed, manuscript figures can be generated with the corresponding scripts:

```text
R/08_Figure1.R
R/09_Figure2.R
R/10_Figure3.R
R/11_Figure4.R
R/12_Figure4_DFK.R
R/13_Figure6.R
R/14_SuppFigure2.R
R/15_SuppFigure4.R
R/16_SuppFigure6.R
```

In general, plots are written under:

```text
results/plots/
```

and figure source-data tables under:

```text
results/tables/
```

## Figure 2 trajectory analysis

Figure 2 uses Monocle3 and requires explicit selection of trajectory root nodes.

First run:

```bash
Rscript R/09_Figure2.R
```

To inspect the principal graph and node labels:

```bash
Rscript R/09_Figure2_select_root_static.R
```

The root nodes selected for the final analysis are:

```text
Y_3
Y_17
```

Save the selected roots with:

```bash
Rscript R/09_Figure2_select_root_static.R Y_3,Y_17
```

Then rerun Figure 2:

```bash
Rscript R/09_Figure2.R
```

Fibroblast and mesothelial cells are retained during graph construction, root-node selection and pseudotime calculation. Mesothelial cells are excluded only from the final fibroblast-specific trajectory outputs.

## Supplementary Figure 2 external validation

Supplementary Figure 2D-F requires external single-cell datasets that are not included in this repository.

See:

```text
data/SuppFigure2/README.md
```

for the original studies, data-access information and expected local paths.

Supplementary Figure 2A-C use the primary study dataset and do not require these external objects.

## Data availability

Large primary study data are not distributed through GitHub.

The primary dataset will be deposited in an appropriate public repository, and accession information will be added when available.

Previously published datasets used for external validation are not redistributed here. Their provenance and access information are documented in `data/SuppFigure2/README.md`.

## Citation

Citation information for the associated manuscript will be added upon publication.
