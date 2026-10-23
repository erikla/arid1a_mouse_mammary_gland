# ARID1A multiome (RNA) analysis with Seurat

R code for the RNA modality of 10x Genomics Multiome (snRNA + snATAC) data from
mouse mammary gland: Arid1a WT vs KO, with and without estradiol (E2). ATAC
processing is done separately with ArchR; ArchR is used here only to keep RNA
cells that also pass ATAC QC.

## Contents

| File | Description |
|---|---|
| `arid1a_seurat.R` | Configuration, `mergeMultiomeSeurat()` (build/QC/normalize/cluster), `main()` (downstream analysis), analysis functions and plotting helpers |
| `venn_helpers.R` | `plotVenn2()`, `plotVenn3()` |

`mergeMultiomeSeurat()` and `main()` are exploratory, section-by-section
analyses (see the `#### Section ####` headers). They are meant to be run
interactively, not top to bottom.

## Configuration

Paths are set through environment variables (e.g. in `~/.Renviron`); defaults are
relative to `ARID1A_BASE_DIR`:

| Variable | Purpose | Default |
|---|---|---|
| `ARID1A_BASE_DIR` | project root (`robj/`, `plot/` created here) | `.` |
| `ARID1A_RAW_DIR` | Cell Ranger ARC output, one directory per sample | `<base>/cellranger` |
| `ARID1A_ARCHR_DIR` | ArchR project (ATAC QC cell list) | `<base>/ArchR` |
| `ARID1A_DATA_DIR` | input tables (marker spreadsheets, scanpy exports) | `<base>/data` |
| `ARID1A_SHEETS_DIR` | output spreadsheets | `<base>/spreadsheets` |
| `ARID1A_ANNOT_DIR` | ribosomal genes, cell-cycle genes, GMT files | `<base>/annotation` |
| `ARID1A_LIB_DIR` | directory containing `venn_helpers.R` | `.` |
| `GPROFILER_HALLMARK_ID`, `GPROFILER_G5_ID` | gProfiler custom-organism IDs for `hotspot_go()` (create with `gprofiler2::upload_GMT_file()`) | unset |

## Dependencies

R >= 4.2 with Seurat (v4/v5-style `@assays` access), SeuratObject,
SingleCellExperiment, glmGamPoi, sctransform, ggplot2, gridExtra, cowplot,
patchwork, data.table, magrittr, ComplexHeatmap, circlize, pals, viridis,
ggrastr, ggrepel, miloR, presto, gprofiler2, writexl, readxl, clustree, MASS,
plotly, Hmisc, dplyr.

`venn_helpers.R` contains the Venn diagram helpers used in `main()` (needs VennDiagram).

## Notes

* Sample MMO2_CRE_E2_Estrogen_multiome is excluded from the analysis.
* The raw sequencing data are not included.

## License

GPL-3.0-or-later (see `LICENSE`). Parts of the plotting code are adapted from Seurat (GPL-3).
