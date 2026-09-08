# External datasets for Supplementary Figure 2

Supplementary Figure 2D-F uses previously published single-cell RNA-seq
datasets for external validation of the fibroblast transcriptional
signatures defined in this study.

The external datasets used in Supplementary Figure 2 are not
redistributed in this repository. Publicly available datasets should be
obtained from their original repositories, whereas the processed Cardoso
et al. object used here was provided directly by the authors.

`R/14_SuppFigure2.R` expects the files to be arranged as follows:

```text
data/SuppFigure2/Fibros.External/
├── Cardoso/
│   └── Total_Integration_fibroblasts.rds
├── Shannon/
│   ├── Mouse_SS_Fibro.RDS
│   ├── Mouse_PS_Fibro.RDS
│   └── Human_PS_Fibro.RDS
└── Sheppard/
    └── GSM6428697_seurat_v4_obj.rds
```

## Buechler et al. fibroblast atlas

Supplementary Figure 2D uses the fibroblast atlases described in:

Buechler MB, Pradhan RN, Krishnamurty AT, et al.
Cross-tissue organization of the fibroblast lineage.
Nature. 2021;593:575-579.
DOI: 10.1038/s41586-021-03549-5

The study integrates fibroblasts from multiple mouse and human tissues
and disease states. The processed fibroblast resources were made
available by the authors through the fibroXplorer resource.

The analysis expects the following processed objects:

```text
data/SuppFigure2/Fibros.External/Shannon/Mouse_SS_Fibro.RDS
data/SuppFigure2/Fibros.External/Shannon/Mouse_PS_Fibro.RDS
data/SuppFigure2/Fibros.External/Shannon/Human_PS_Fibro.RDS
```

These correspond to the mouse steady-state, mouse perturbed-state and
human perturbed-state fibroblast atlases, respectively.

## Cardoso et al.

Supplementary Figure 2E uses fibroblast data associated with:

Cardoso EC, Lee H, England FJ, et al.
Early fibrotic niches establish tumour-permissive microenvironments.
Nature. 2026.

The processed fibroblast object used in this analysis was provided
directly by the authors and is not redistributed in this repository.

The analysis expects the object at:

```text
data/SuppFigure2/Fibros.External/Cardoso/Total_Integration_fibroblasts.rds
```

Access to this processed object should be requested from the authors of
the original study unless a public version becomes available in the
future.

## Tsukui et al.

Supplementary Figure 2F uses data from:

Tsukui T, Wolters PJ, Sheppard D.
Alveolar fibroblast lineage orchestrates lung inflammation and fibrosis.
Nature. 2024.
DOI: 10.1038/s41586-024-07660-1

The single-cell RNA-seq dataset generated in that study is publicly
available under:

GSE210341

The processed all-sample dataset used here corresponds to:

GSM6428697

The analysis expects the processed Seurat object at:

```text
data/SuppFigure2/Fibros.External/Sheppard/GSM6428697_seurat_v4_obj.rds
```

## Running the external validation

Once the external objects have been placed in the paths above, run from
the repository root:

```bash
Rscript R/14_SuppFigure2.R
```

These external datasets are required only for Supplementary Figure 2D-F.
They are not required for the primary single-cell analysis or for
Supplementary Figure 2A-C.
