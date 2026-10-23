# ARID1A loss in mouse mammary gland: single-nucleus multiome analysis

Code for the analysis of 10x Genomics Multiome (snRNA-seq + snATAC-seq) data from
mouse mammary gland: Arid1a WT vs KO, with and without estradiol (E2).

| Directory | Modality | Tool |
|---|---|---|
| [`seurat/`](seurat) | RNA | Seurat |
| [`archr/`](archr) | ATAC | ArchR |

The RNA analysis is run first; the ArchR project imports the resulting Seurat
object for RNA-ATAC integration. See the README in each directory for
configuration and dependencies.

## License

GPL-3.0-or-later, see [`LICENSE`](LICENSE).
