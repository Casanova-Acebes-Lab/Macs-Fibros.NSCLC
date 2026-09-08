# Primary input data

The primary sample-level input data used in this analysis are not stored
in this GitHub repository because of their size.

These data will be deposited in an appropriate public repository.
The accession number and download information will be added here once
the dataset becomes publicly available.

## Expected directory structure

`R/01_QC.R` expects the sample-level R objects to be organized as:

```text
data/raw/
├── CCR2KO/
│   ├── m1.Bis.rds
│   ├── m2.Bis.rds
│   ├── m3.Bis.rds
│   └── m4.Bis.rds
├── WT/
│   ├── m10.rds
│   ├── m12.rds
│   ├── m13.rds
│   └── m14.rds
└── CD169DTR/
    ├── m15.rds
    ├── m16.rds
    ├── m17.rds
    └── m18.rds

```

Small doublet-barcode lists required during QC are versioned separately
under `data/barcodes.doublets/`.

After downloading the primary data and placing the files in the
directories above, run the analysis from the repository root.

Public data accession: **pending**
