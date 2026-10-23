# ARID1A multiome (ATAC) analysis with ArchR

R code for the ATAC modality of 10x Genomics Multiome (snRNA + snATAC) data from
mouse mammary gland: Arid1a WT vs KO, with and without estradiol (E2). The RNA
modality is analysed with Seurat (separate project); its object is imported here
for RNA-ATAC integration.

## Contents

| File | Description |
|---|---|
| `arid1a_multiome_archR.R` | Configuration, `BuildProj()` (create ArchR project), `main()` (downstream analysis), analysis functions, plotting helpers and patched ArchR internals |

`BuildProj()` and `main()` are exploratory, section-by-section analyses (see the
`#### Section ####` headers). They are meant to be run interactively, not top to
bottom.

## Configuration

Paths are set through environment variables (e.g. in `~/.Renviron`); defaults are
relative to `ARID1A_BASE_DIR`:

| Variable | Purpose | Default |
|---|---|---|
| `ARID1A_BASE_DIR` | project root (`archR/`, `MultiomeSeurat/`, `spreadsheets/`) | `.` |
| `ARID1A_RAW_DIR` | Cell Ranger ARC output, one directory per sample | `<base>/cellranger` |
| `ARID1A_ANNOT_DIR` | GTF, Vierstra motifs, ribosomal / mitochondrial gene lists | `<base>/annotation` |
| `ARID1A_SEURAT_RDS` | Seurat object from the RNA analysis | `<base>/MultiomeSeurat/robj/arid1a_sctransform_multiome_regress.rds` |
| `MACS2_PATH` | macs2 executable | `macs2` on the PATH |

## Dependencies

R >= 4.2, [ArchR](https://github.com/GreenleafLab/ArchR) (with
`ArchR::installExtraPackages()`), macs2, and: BSgenome.Mmusculus.UCSC.mm10,
ComplexHeatmap, circlize, Seurat, SingleCellExperiment, TFBSTools, chromVARmotifs,
ggseqlogo, seqLogo, Polychrome, pals, ggrepel, shadowtext, ggrastr, ggh4x,
data.table, magrittr, rtracklayer, biomaRt, org.Hs.eg.db, org.Mm.eg.db,
Orthology.eg.db, igraph, doMC, doParallel, foreach, mgcv, XICOR, pheatmap, plotly,
magick, metap.

## Notes

* Sample MMO2 is excluded from the analysis.
* The raw sequencing data are not included.
* Some functions are modified copies of ArchR code (MIT license).

## License

GPL-3.0-or-later (see `LICENSE`).
