#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
#################################################################################
# ARID1A multiome (ATAC) analysis with ArchR
#
# Processes the ATAC modality of 10x Genomics Multiome (RNA + ATAC) samples from
# mouse mammary gland: Arid1a WT vs KO, with and without estradiol (E2).
# The RNA modality is processed separately with Seurat (see the MultiomeSeurat
# project); the Seurat object is imported here for RNA-ATAC integration.
#
# Contents
#   1. Configuration (paths, thresholds, marker gene lists, palettes)
#   2. BuildProj()            - create the ArchR project from Cell Ranger output
#   3. main()                 - downstream analysis: clustering, peaks, markers,
#                               motif / chromVAR analysis, integration with RNA
#   4. Analysis functions     - correlateAccessibility*(), NonRedundantMotifs(), ...
#   5. Plotting helpers       - UmapPlotEnhance(), DoHeatmapEnhance(), ...
#   6. Patched ArchR functions (addDeviationsMatrix, filterDoubletsFix, ...)
#
# main() and BuildProj() are exploratory notebooks-in-a-function. They are not
# meant to be run top to bottom: run them section by section (each
# "#### Section ####" header) in an interactive session.
#
# Requires: R >= 4.2, ArchR (GreenleafLab/ArchR), macs2 on the PATH, and the
# packages loaded below / listed in the README.
################################################################################

#### Libraries ####
library(ArchR)
library(ComplexHeatmap)
library(gridExtra)
library(ggrepel)
library(shadowtext)
library(BSgenome.Mmusculus.UCSC.mm10)
library(data.table)
library(magrittr)

options(bitmapType = 'cairo')
addArchRThreads(threads = 8)
addArchRGenome("mm10")

#### Configuration ####
# All locations are configurable. Set the environment variables below (e.g. in
# ~/.Renviron) or edit the defaults. Relative defaults are resolved against the
# directory the script is run from.
#   ARID1A_BASE_DIR   project root (contains archR/, MultiomeSeurat/, spreadsheets/)
#   ARID1A_RAW_DIR    Cell Ranger ARC output; one sub-directory per sample
#   ARID1A_ANNOT_DIR  annotation files (GTF, Vierstra motifs, ribosomal/mito genes)
#   ARID1A_SEURAT_RDS Seurat object produced by the MultiomeSeurat analysis
#   MACS2_PATH        macs2 executable (default: macs2 on the PATH)
getDir <- function(var, default) {
  normalizePath(Sys.getenv(var, unset = default), mustWork = FALSE)
}
base.dir         <- getDir('ARID1A_BASE_DIR', '.')
archr.dir        <- file.path(base.dir, 'archR')
data.dir         <- file.path(archr.dir, 'data')
plot.dir         <- file.path(archr.dir, 'Plots')
robj.dir         <- file.path(archr.dir, 'robj')
spreadsheets.dir <- file.path(base.dir, 'spreadsheets')
raw.dir          <- getDir('ARID1A_RAW_DIR', file.path(base.dir, 'cellranger'))
annot.dir        <- getDir('ARID1A_ANNOT_DIR', file.path(base.dir, 'annotation'))
multiRNA.file    <- getDir('ARID1A_SEURAT_RDS',
                           file.path(base.dir, 'MultiomeSeurat', 'robj',
                                     'arid1a_sctransform_multiome_regress.rds'))
pathToMacs2      <- Sys.getenv('MACS2_PATH', unset = Sys.which('macs2'))
for (d in c(data.dir, plot.dir, robj.dir, spreadsheets.dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Minimum number of ATAC fragments per cell (cell filtering and QC plots)
nfragMin <- 1000

# Cell Ranger ARC mm10 gene annotation
gtf.file <- file.path(annot.dir, 'genes.gtf.gz')

markers_list <- list(
  #https://www.nature.com/articles/s41467-020-15548-7#Sec2 (Supplemental figure 7)
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4331521/
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5696379/ - Visavader
  'Mammary Luminal progenitor'=c('Anpep','Brca1','Cd14','Cd47','Elf5','Egfr','Ltf','Mfge8','Prom1','Tnfrsf11a'),
  'Mammary Luminal mature'=c('Cd24a','Epcam',"Erbb2","Esr1","Foxa1","Gata3",'Krt8','Krt18',"Krt19","Muc1","Pgr","Trp53"),
  'Mammary basal'=c('Anxa8','Itga6','Krt4','Krt5','Krt14','Krt15','Krt17','Mki67','Trp63','Trp73'),
  ## zoom Andrew Ewald
  'mammary epithelium 2-week'=c('Sparc','Tgfbr3','Lgals1','Sorbs1','Timp3','Scn7a','Vim','Cadm1','Prdm16'), #Rxrg
  'mammary epithelium mature luminal'=c('Areg','Cited1','Cldn4','Epcam','Krt8','Ly6d','Prlr','Ptn','Stc2','Krt18'),
  'mammary epithelium progenitor'=c('Aldh1a3','Cd14','Cd81','Csn1s1','Csn3','Elf5','Jund','Lalba','Muc15','Tspan8','Wnt5b'),
  'mammary epithelium basal'=c('Fstl1','Palm','Cd55','Adamts2','Lrp1','Col5a1','Tpm2','Myl9','Myh11','Cnn1','Krt14','Acta2',
                               'Krt5','Cxcl14','Krt17','Snai2','Col18a1','Vcan','Nrg1','Lama1','Wif1','Bmp7','Tgfb1i1',
                               'Cacna2d1','Tenm2','Synpo2'), #'Serpinb11', 'Mrvi1'
  'AR dependent'=c('Ar','Fkbp5','Hpgd','Fgfr2','Lgr4','Egf','Egfr'),
  'Cell Cycle'=c("Pcna","Mki67","Top2a"),
  'Bach markers'=c('Cited1','Prlr','Esr1','Areg','Foxa1','Ly6a','Aldh1a3','Kit',
                   'Cd14','Lypd3','Rspo1','Atp6v1c2','Fabp3','Thrsp','Wap','Glycam1','Olah','1500015O10Rik',
                   'Col7a1','Moxd1','Mia','Emid1','Pdpn','Col9a2','Fbln2','Igfbp3','Fst','Li17b','Oxtr','Krt15','Igfbp6','Igfbp2',
                   'Tns1','Gng11','Procr','Igfbp7','Nrip2','Notch3'),
  'KO' = c("Arid1a")
)

markers_list[['Mammary Luminal mature']] <- sort(unique(c(markers_list[['Mammary Luminal mature']],
                                                          markers_list[['mammary epithelium mature luminal']])))
markers_list[['Mammary Luminal progenitor']] <- sort(unique(c(markers_list[['Mammary Luminal progenitor']],
                                                              markers_list[['mammary epithelium progenitor']])))
markers_list[['Mammary basal']] <- sort(unique(c(markers_list[['Mammary basal']],
                                                 markers_list[['mammary epithelium basal']])))

# To find whether genes exist
# gs <- getFeatures(ArchRProj = proj, useMatrix = 'GeneScoreMatrix')
# for(i in markers_list){ lapply(i, function(x) { n <- length(grep(x, gs, ignore.case = T, value = T)); if(n==0) print(x) } );  }
# grep("Foxai1", gs, ignore.case = T, value = T)

#### Color palettes (continuous) ####
rocket_palette <- c("linen", viridis::viridis(direction = -1, option = "F", n=20))
solarExtraWhite_palette <- c("#3361A5","#248AF3","#14B3FF","#88CEEF","white","#EAD397","#FDB31A","#E42A2A","#A31D1D")
solarHalf_palette <- c("#F2F2F2","#C1D5DC","#DEE5AE","#EAD397","#FDB31A","#E42A2A","#A31D1D")
heatmap_palette <- colorRampPalette(c("royalblue4","royalblue2","deepskyblue",
                                                 "grey90","orangered2","#ce1300","darkred"), space="rgb")(1024)
color_blind_heatmap_palette <- colorRampPalette(c("dodgerblue4","dodgerblue3","dodgerblue2","dodgerblue1",
                                                  "grey98","darkorange1","darkorange2","darkorange3","darkorange4"), space="rgb")(120)
blue_gray_purple_diverging_palette <- c(rev(pals::brewer.pubu(10)[-c(1:2)]),'grey99',pals::brewer.purd(10)[-c(1:2)]) # pals::coolwarm(25),
godsnot_102 = c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#6A3A4C","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900",
                "#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75",
                "#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1",
                "#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72")

#### Utilities ####
# Map human gene symbols to mouse orthologs (Orthology.eg.db).
humanGeneList2MouseOrthologs <- function(humangenes){
  require(Orthology.eg.db)
  require(org.Mm.eg.db)
  require(org.Hs.eg.db)
  gns <- mapIds(org.Hs.eg.db, humangenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Mus_musculus","Homo_sapiens")
  naind <- is.na(mapped$Mus_musculus)
  msymb <- mapIds(org.Mm.eg.db, as.character(mapped$Mus_musculus[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Human_symbol = humangenes, mapped, Mouse_symbol = NA)
  out$Mouse_symbol[!naind] <- msymb
  print(head(out))
  out
}

## thresholded peaks
other.peak.cutoff <- "FDR <= 0.05 & abs(Log2FC) >= 0.5"

marker.cutOff <- "FDR <= 0.05 & Log2FC >= 0.5"
marker.sample.cutOff <- "FDR <= 0.05 & Log2FC >= 0.5"
abs.marker.sample.cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5"

marker.cluster.cutOff <- marker.sample.cutOff
abs.marker.cluster.cutOff = abs.marker.sample.cutOff

gene.marker.sample.cutOff = "FDR <= 0.01 & Log2FC >= 1.25"

gist_heat_r <- c('#ffffff', '#fffdfb', '#fffbf7', '#fff9f3', '#fff7ef', '#fff5eb', '#fff3e7', '#fff1e3', '#ffefdf', '#ffeddb', '#ffebd7', '#ffe9d3', '#ffe7cf', '#ffe5cb', '#ffe3c7', '#ffe1c3', '#ffdfbf', '#ffddbb', '#ffdbb7', '#ffd9b3', '#ffd7af', '#ffd5ab', '#ffd3a7', '#ffd1a3',
                 '#ffcf9f', '#ffcd9b', '#ffcb97', '#ffc993', '#ffc78f', '#ffc58b', '#ffc387', '#ffc183', '#ffbf7f', '#ffbd7b', '#ffbb77', '#ffb973', '#ffb76f', '#ffb56b', '#ffb367', '#ffb163', '#ffaf5f', '#ffad5b', '#ffab57', '#ffa953', '#ffa74f', '#ffa54b', '#ffa347', '#ffa143',
                 '#ff9f3f', '#ff9d3b', '#ff9b37', '#ff9933', '#ff972f', '#ff952b', '#ff9327', '#ff9123', '#ff8f1f', '#ff8d1b', '#ff8b17', '#ff8913', '#ff870f', '#ff850b', '#ff8307', '#ff8103', '#ff7f00', '#ff7d00', '#ff7b00', '#ff7900', '#ff7700', '#ff7500', '#ff7300', '#ff7100',
                 '#ff6f00', '#ff6d00', '#ff6b00', '#ff6900', '#ff6700', '#ff6500', '#ff6300', '#ff6100', '#ff5f00', '#ff5d00', '#ff5b00', '#ff5900', '#ff5700', '#ff5500', '#fe5300', '#fc5100', '#fb4f00', '#f94d00', '#f84b00', '#f64900', '#f44700', '#f34500', '#f24300', '#f04100',
                 '#ee3f00', '#ed3d00', '#ec3b00', '#ea3900', '#e93700', '#e73500', '#e53300', '#e43100', '#e22f00', '#e12d00', '#e02b00', '#de2900', '#dc2700', '#db2500', '#da2300', '#d82100', '#d61f00', '#d51d00', '#d41b00', '#d21900', '#d11700', '#cf1500', '#ce1300', '#cc1100',
                 '#cb0f00', '#c90d00', '#c80b00', '#c60900', '#c40700', '#c30500', '#c20300', '#c00100', '#be0000', '#bd0000', '#bc0000', '#ba0000', '#b90000', '#b70000', '#b50000', '#b40000', '#b20000', '#b10000', '#b00000', '#ae0000', '#ac0000', '#ab0000', '#aa0000', '#a80000',
                 '#a60000', '#a50000', '#a40000', '#a20000', '#a10000', '#9f0000', '#9e0000', '#9c0000', '#9b0000', '#990000', '#980000', '#960000', '#940000', '#930000', '#920000', '#900000', '#8e0000', '#8d0000', '#8c0000', '#8a0000', '#890000', '#870000', '#850000', '#840000',
                 '#820000', '#810000', '#800000', '#7e0000', '#7c0000', '#7b0000', '#7a0000', '#780000', '#760000', '#750000', '#740000', '#720000', '#710000', '#6f0000', '#6e0000', '#6c0000', '#6a0000', '#690000', '#680000', '#660000', '#640000', '#630000', '#620000', '#600000',
                 '#5e0000', '#5d0000', '#5c0000', '#5a0000', '#590000', '#570000', '#560000', '#540000', '#520000', '#510000', '#500000', '#4e0000', '#4c0000', '#4b0000', '#4a0000', '#480000', '#460000', '#450000', '#440000', '#420000', '#410000', '#3f0000', '#3d0000', '#3c0000',
                 '#3b0000', '#390000', '#380000', '#360000', '#340000', '#330000', '#310000', '#300000', '#2f0000', '#2d0000', '#2c0000', '#2a0000', '#290000', '#270000', '#250000', '#240000', '#230000', '#210000', '#200000', '#1e0000', '#1c0000', '#1b0000', '#190000', '#180000',
                 '#170000', '#150000', '#140000', '#120000', '#110000', '#0f0000', '#0d0000', '#0c0000', '#0b0000', '#090000', '#080000', '#060000', '#040000', '#030000', '#010000', '#000000')

solarHalf <- c("#F2F2F2","#C1D5DC","#DEE5AE","#EAD397","#FDB31A","#E42A2A","#A31D1D")

INPUT_LABELS = list(
  "MMO1_Control_multiome" = "Arid1a",
  "MMO1_E2_Estrogen_multiome" = "Arid1a E2",
  "MMO1_CRE_Control_multiome" = "Arid1a ko",
  "MMO1_CRE_E2_Estrogen_multiome" = "Arid1a ko E2"
  #"MMO2_CRE_E2_Estrogen_multiome" = "Arid1a(-/-) E2 mm02"
)

BIONAMES_LABELS = list(
  "MMO1_Control_multiome" = "Arid1a(+/+)",
  "MMO1_E2_Estrogen_multiome" = "Arid1a(+/+) E2",
  "MMO1_CRE_Control_multiome" = "Arid1a(-/-)",
  "MMO1_CRE_E2_Estrogen_multiome" = "Arid1a(-/-) E2"
  #"MMO2_CRE_E2_Estrogen_multiome" = "Arid1a(-/-) E2 mm02"
)

CONDITION_LABELS <- list(
  "MMO1_Control_multiome" = "WT",
  "MMO1_E2_Estrogen_multiome" = "WT",
  "MMO1_CRE_Control_multiome" = "KO",
  "MMO1_CRE_E2_Estrogen_multiome" = "KO"
  #"MMO2_CRE_E2_Estrogen_multiome" = "KO"
)
condition_palette <- c('KO'='#CBA00D','WT'='#272E6A')
condition_palette_list <- as.list(condition_palette)

TREATMENT_LABELS <- list(
  "MMO1_Control_multiome" = "none",
  "MMO1_E2_Estrogen_multiome" = "E2",
  "MMO1_CRE_Control_multiome" = "none",
  "MMO1_CRE_E2_Estrogen_multiome" = "E2"
  #"MMO2_CRE_E2_Estrogen_multiome" = "E2"
)
treatment_palette <- c("none"="orange", "E2"="blue")

bionames_palette_list <- sample_palette <- c('#03A45D', '#F9A619', '#EE1C24', '#203F9A')
names(bionames_palette_list) <- c('Arid1a(+/+)','Arid1a(+/+) E2','Arid1a(-/-)','Arid1a(-/-) E2')
levels(bionames_palette_list) <- names(bionames_palette_list)

luminal_basal_palette <- c("#222222","#F3C300","#875692","gray80") # is pals::kelly()-1
names(luminal_basal_palette) <- c("Mammary basal","Mammary Luminal mature","Mammary Luminal progenitor", "Unknown")
levels(luminal_basal_palette) <- names(luminal_basal_palette)

mammary_cluster_class_palette <- c("#F0A0FF", "#0075DC", "#006633", "#F3C300", "#875692", "#CCCCCC")
names(mammary_cluster_class_palette) <- c("Basal/MaSC","L-Alv","LA-Pro","L-Hor","LH-Pro", "Unknown")
mammary_cluster_class_palette_list <- as.list(mammary_cluster_class_palette)

mammary_refine_clusters_palette <- c("#F0A0FF","#0075DC","#F3C300","#006633","#875692","#bd0000") #, "#CCCCCC")
names(mammary_refine_clusters_palette) <- c("Basal/MaSC","L-Alv","L-Hor","LA-Pro","LH-Pro","L-Hor-Arid1a-") #, NA)
mammary_refine_clusters_palette_list <- as.list(mammary_refine_clusters_palette)

cluster_palette <- c(ArchRPalettes$stallion, ArchRPalettes$stallion2)
names(cluster_palette) <- paste0("C",seq_along(cluster_palette))
levels(cluster_palette) <- names(cluster_palette)
cluster_palette_list <- as.list(cluster_palette)

rna_cluster_palette <- list(
  '1'= "#F0A0FF",
  '2'= "#0075DC",
  '3'= "#993F00",
  '4'= "#4C005C",
  '5'= "#191919",
  '6'= "#005C31",
  '7'= "#2BCE48",
  '8'= "#FFCC99",
  '9'= "#808080",
  '10'= "#94FFB5",
  '11'="#8F7C00",
  '12'="#9DCC00",
  '13'="#C20088",
  '14'="#003380",
  '15'="#FFA405",
  '16'="#FFA8BB",
  '17'="#426600",
  '18'="#FF0010",
  '19'="#5EF1F2",
  '20'="#00998F",
  '21'="#E0FF66",
  '22'="#740AFF",
  '23'="#990000",
  '24'="#FFFF80",
  '25'="#FFE100",
  '26'="#FF5005",
  '-1' = '#CCCCCC'
)
rna_cluster_palette <- unlist(rna_cluster_palette)

#' Distinct colours for the TF families of expressed motifs
#'
#' @param proj ArchR project with a `VierstraMotifToGene` table (see Motif2GeneLookupDT()).
#' @return Named colour vector, one colour per motif family.
makeMotifFamilyPalette <- function(proj){
  library(Polychrome)
  dtM2Glookup <- proj@VierstraMotifToGene
  n <- length( unique(dtM2Glookup[expressed==T]$family) )
  # create your own color palette based on `seedcolors`
  family_palette = createPalette(n,  c("#ff0000", "#00ff00", "#0000ff"))
  #swatch(family_palette)
  names(family_palette) <- unique(dtM2Glookup[expressed==T]$family)
  family_palette
}

#' Read a GTF file into a data.table
#'
#' @param gtf.file Path to a (gzipped) GTF file.
#' @param typeof Feature type to keep (e.g. "gene", "exon").
importGTF <- function(gtf.file, typeof='gene'){
  library(rtracklayer) ## Import GTF
  #gtf <- import(gtf.file, colnames = c('type','gene_id','gene_type','gene_name','exon_number','transcript_id'))
  message("reading gtf file:",gtf.file)
  gtf <- import(gtf.file)
  gtf <- unique(as.data.table(gtf))
  gtf[type==typeof]
}

#' Downstream ATAC analysis of the multiome ArchR project
#'
#' Exploratory analysis notebook, organised in "####" sections: QC, clustering
#' and UMAPs, peak calling (MACS2), marker peaks and genes, differential
#' accessibility (KO vs WT, E2), motif enrichment and chromVAR deviations,
#' integration with the Seurat RNA object, and export of results. Requires the
#' project created by BuildProj().
#'
#' Run section by section in an interactive session. Figures are written to
#' `plot.dir`, tables to `spreadsheets.dir`.
main <- function(){

  # BuildProj()

  #### Load object ####
  proj <- loadArchRProject(path = archr.dir)

  family_palette <- makeMotifFamilyPalette(proj = proj)

  #### Save object ####
  saveArchRProject(ArchRProj = proj, load = TRUE)

  getAvailableMatrices(proj)
  ArchRBrowser(proj)

  # Look at data
  proj
  proj@cellColData
  table(proj@cellColData$bioNames)

  #### Plot final UMAPS ####
  plabel <- plotEmbedding(ArchRProj = proj,
                          colorBy = "cellColData",
                          name = "Clusters",
                          embedding = "UMAP") %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    addLabels = T,
                    title.size = 18,
                    panel.border.size = 0,
                    pal = NA,
                    label.size = 8)

  psample <- plotEmbedding(ArchRProj = proj,
                           colorBy = "cellColData",
                           name = "gex_bioNames",
                           embedding = "UMAP",
                           baseSize = 12) %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    axis.label = 'color',
                    pt.size = 1.1,
                    panel.border.size = 0,
                    #nrow=2,
                    title = "organoid samples",
                    pal=bionames_palette_list,
                    title.size = 18)

  pdf(file = file.path(plot.dir,"UMAP-ATAC-bioNames.pdf"), width = 6, height = 4,
      bg = 'transparent')
  psample + plabel$layers
  dev.off()

  pbioNamesfacet <- plotEmbedding(ArchRProj = proj,
                                   colorBy = "cellColData",
                                   name = "gex_bioNames",
                                   embedding = "UMAP",
                                   size=0.2,
                                   baseSize = 11) %>%
    UmapPlotEnhance(ArchRProj = proj, p1 = .,
                    addLabels = "Clusters",
                    axis.label = "gex_bioNames",
                    pt.size = .9,
                    nrow = 2,
                    order.random = F,
                    title.size = 18,
                    panel.border.size = 0,
                    legend.title = "snATAC",
                    legend.title.rotation = 0,
                    pal = bionames_palette_list)

  pdf(file = file.path(plot.dir,"UMAP-ATAC-gex_bioNames_facet.pdf"),
      width = 7, height = 6,
      bg = 'transparent')
  pbioNamesfacet
  dev.off()

  #Cell state
  pmammarystate <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                          name = "gex_mammary_gland_cluster_classification",
                          embedding = "UMAP", size=0.2, baseSize = 11) %>%
    UmapPlotEnhance(ArchRProj = proj, p1 = .,
                    addLabels = F,
                    pt.size = .9,
                    title.size = 18,
                    panel.border.size = 0,
                    legend.title = "snATAC",
                    legend.title.rotation = 0,
                    pal = mammary_cluster_class_palette)
  pdf(file = file.path(plot.dir,"UMAP-ATAC-mammary_cell_classification.pdf"), width = 6, height = 4,
      bg = 'transparent')
  pmammarystate + plabel$layers
  dev.off()

  #Refined cell state
  pcellstate <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                                   name = "gex_mammary_gland_refine_clusters",
                                   embedding = "UMAP", size=0.2, baseSize = 11) %>%
    UmapPlotEnhance(ArchRProj = proj, p1 = .,
                    addLabels = F,
                    pt.size = .9,
                    title.size = 18,
                    panel.border.size = 0,
                    #title = 'mammary cluster classification',
                    legend.title = "snATAC",
                    legend.title.rotation = 0,
                    pal = mammary_refine_clusters_palette)

  pdf(file = file.path(plot.dir,"UMAP-ATAC-mammary_cell_refine.pdf"), width = 6, height = 4,
      bg = 'transparent')
  pcellstate + plabel$layers
  dev.off()

  pcellstatefacet <- plotEmbedding(ArchRProj = proj,
                                   colorBy = "cellColData",
                                   name = "gex_mammary_gland_refine_clusters",
                                   embedding = "UMAP",
                                   size=0.2,
                                   baseSize = 11) %>%
    UmapPlotEnhance(ArchRProj = proj, p1 = .,
                    addLabels = "Clusters",
                    axis.label = "gex_mammary_gland_refine_clusters",
                    pt.size = .9,
                    nrow = 3,
                    order.random = F,
                    title.size = 18,
                    panel.border.size = 0,
                    legend.title = "snATAC",
                    legend.title.rotation = 0,
                    pal = mammary_refine_clusters_palette)

  pdf(file = file.path(plot.dir,"UMAP-ATAC-mammary_cell_refine_facet.pdf"),
      width = 7, height = 10,
      bg = 'transparent')
  pcellstatefacet
  dev.off()

  #Cluster
  pcluster <- plotEmbedding(ArchRProj = proj,
                colorBy = "cellColData",
                name = "Clusters",
                embedding = "UMAP") %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    pt.size = 1.1,
                    addLabels = F,
                    show.legend = T,
                    title.size = 18,
                    pal = cluster_palette,
                    title = 'organoid clusters',
                    label.size = 8)
  png(filename = file.path(plot.dir,"UMAP-ATAC-clusters.png"), width = 6, height = 4,
      bg = 'transparent', res = 150, units = 'in')
  pdf(file = file.path(plot.dir,"UMAP-ATAC-clusters.pdf"), width = 5.5, height = 4, bg = 'transparent')
  pcluster #+ plabel$layers
  dev.off()

  #scRNA overlap
  pcluster <- plotEmbedding(ArchRProj = proj,
                            colorBy = "cellColData",
                            name = "gex_seurat_clusters",
                            embedding = "UMAP") %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    pt.size = 1,
                    addLabels = F,
                    show.legend = T,
                    title.size = 18,
                    pal = seurat_cluster_palette,
                    title = 'RNA clusters on ATAC embedding',
                    label.size = 8)
  pdf(file = file.path(plot.dir,"UMAP-ATAC-RNAclusters.pdf"), width = 6, height = 4, bg = 'transparent')
  pcluster + plabel$layers
  dev.off()

  pgeno <- plotEmbedding(ArchRProj = proj,
                         colorBy = "cellColData",
                         name = "Condition",
                         embedding = "UMAP") %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    pt.size = 1.1,
                    addLabels = "Clusters",
                    axis.label = "Condition",
                    label.size = 7,
                    title.size = 18,
                    nrow=2,
                    title='organoid genotypes',
                    pal = condition_palette)
  pdf(file = file.path(plot.dir,"UMAP-ATAC-genotype_facet.pdf"), width = 4, height = 6, bg = 'transparent')
  pgeno
  dev.off()

  #condition
  pcondition <- plotEmbedding(ArchRProj = proj,
                           colorBy = "cellColData",
                           name = "Condition",
                           embedding = "UMAP",
                           baseSize = 12) %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    addLabels = F,
                    pt.size = 1.1,
                    panel.border.size = 0,
                    title = "organoid conditions",
                    pal=condition_palette,
                    title.size = 18)

  png(filename = file.path(plot.dir,"UMAP-ATAC-Condition.png"), width = 6, height = 4,
      bg = 'transparent', res = 150, units = 'in')
  pcondition + plabel$layers
  dev.off()

  pcondition <- plotEmbedding(ArchRProj = proj,
                              colorBy = "cellColData",
                              name = "Condition",
                              embedding = "FR",
                              baseSize = 12) %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    addLabels = F,
                    pt.size = 1.1,
                    panel.border.size = 0,
                    title = "organoid conditions",
                    pal=condition_palette,
                    title.size = 18)

  pcondition

  p <- plotEmbedding(ArchRProj = proj,
                     colorBy = "cellColData",
                     name = "gex_mammary_gland_refine_clusters",
                     embedding = "FR") %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    pt.size = 1.1,
                    outline.size = 0.02,
                    panel.border.size = 0,
                    addLabels = F,
                    label.size = 7,
                    title.size = 16,
                    order.random = F,
                    title = "Cell types",
                    pal =  mammary_refine_clusters_palette,
                    #pal =  bionames_palette_list,
                    legend.title = paste("cell type"),
                    legend.title.rotation = 0)
  pdf(file = file.path(plot.dir,'FA_gex_mammary_refine.pdf'), width = 5.5, height = 5, bg = 'transparent')
  p
  dev.off()

  ## Summary
  png(filename = file.path(plot.dir,'UMAP-summary.png'), width = 22, height = 5,
      bg = 'transparent', res = 150, units = 'in')
  grid.arrange(psample + plabel$layers, pcluster, pgeno, pmammarystate + plabel$layers, nrow = 1)
  dev.off()

  #nFrag
  dtl <- data.table(proj@embeddings$UMAP$df, log10(proj@cellColData$nFrags)) %>% setnames(.,c("x","y","z","color"))
  dtl[,z:=NULL]
  UmapPlotEnhance(ArchRProj = proj,
                  p1 = dtl,
                  addLabels = F,
                  pt.size = 1.1,
                  draw.border = T,
                  scale_color_bar = F,
                  title = "",
                  legend.title = "log10(nFrags)",
                  filename.prefix = "UMAP-ATAC-nFRAG")

  #### Density UMAP ####

  # lighter to darker
  g <- UmapPlotDensity(ArchrProj = proj, feature = 'Condition',
                       legend.height.factor = .5,
                       contour.color = NA,
                       filename.prefix = 'ConditionDensityPlots',
                       as.pdf = T,
                       fixed_coord = T,
                       nrow = 1, pt.size = 0.5)
  pdf(file = file.path(plot.dir,"ConditionDensityPlots2.pdf"), width=4, height=6)
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  # lighter to darker
  table(proj$bioNames)
  dense.pal <- list(pals::ocean.amp(10),
                    rev(pals::ocean.dense(10)),
                    rev(pals::ocean.algae(10)),
                    rev(pals::ocean.solar(10)))
  names(dense.pal) <- c('Arid1a(-/-)','Arid1a(-/-) E2','Arid1a(+/+)','Arid1a(+/+) E2')

  g <- UmapPlotDensity(ArchrProj = proj, feature = 'bioNames',
                       pal = dense.pal,
                       legend.height.factor = .5,
                       contour.color = NA,
                       filename.prefix = 'BioNamesDensityPlots',
                       as.pdf = T,
                       fixed_coord = T, nrow = 2, pt.size = 0.5)
  pdf(file = file.path(plot.dir,"BioNamesDensityPlots2.pdf"), width=9, height=7)
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  pd <- UmapPlotDensity(ArchrProj = proj,
                        feature = "gex_mammary_gland_refine_clusters",
                        embedding = "FR",
                        contour.color = NA,
                        legend.height.factor = 1.1)
  pdf(file = file.path(plot.dir,"FA_density_maps.pdf"), width = 16, height = 4, bg = "transparent")
  do.call(grid.arrange, c(pd, nrow=1))
  dev.off()

  library(MASS)
  library(plotly)
  df <- proj@embeddings$UMAP$df
  df$mammary_gland_cluster_classification <- proj1$mammary_gland_cluster_classification
  colnames(df) <- c('x','y','z','mammary_gland_cluster_classification')
  df$mammary_gland_cluster_classification <- factor(df$mammary_gland_cluster_classification, levels = names(mammary_cluster_class_palette))
  p <- plot_ly(df, x=df$x, y=df$y, z=df$z,
               color=~mammary_gland_cluster_classification,
               colors = mammary_cluster_class_palette ) %>% add_markers(size=.3)
  p

  table(p1$data$color)

  pcells <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                name = "mammary_gland_cluster_classification",
                embedding = "UMAP", size=0.2, baseSize = 11) %>%
    UmapPlotEnhance(ArchRProj = proj, p1 = ., addLabels = F,
                    pal = mammary_cluster_class_palette, pt.size = .8)
  png(filename = file.path(plot.dir,"UMAP-ATAC-mammary_gland_cluster_classification.png"), width = 7, height = 5,
      bg = 'transparent', res = 150, units = 'in')
  pcells + pcluster$layers
  dev.off()

  PlotBarVariableInClusters(ArchRProj = proj,
                            variable = 'gex_bioNames',
                            palette = bionames_palette_list,
                            orientation = "vertical",
                            height_scale_perc = .75,
                            width_scale_perc = .9,
                            addcounts = T,
                            filename = 'ATAC_bioNames_percent_in_clusters')

  PlotBarVariableInClusters(ArchRProj = proj,
                            variable = 'Condition',
                            palette = condition_palette_list,
                            orientation = "vertical",
                            height_scale_perc = .7,
                            width_scale_perc = .9,
                            addcounts = T,
                            filename = 'ATAC_condition_percent_in_clusters')

  PlotBarVariableInClusters(ArchRProj = proj,
                            variable = 'gex_mammary_gland_refine_clusters',
                            orientation = "vertical",
                            height_scale_perc = .8,
                            width_scale_perc = .9,
                            palette = mammary_refine_clusters_palette,
                            addcounts = T,
                            filename = 'ATAC_mammary_gland_refine_classification_percent_in_clusters')
  #------------------------------------------- Pseudobulk -------------------------------------------
  # Pseudo-Bulk for calling peaks per cluster #
  proj <- addGroupCoverages(ArchRProj = proj,
                            minReplicates = 2,
                            groupBy = "Clusters",
                            force = T)

  # norm by ReadsInTSS is default
  # tilesk$reads <- tilesk$reads * 10^4/sum(normBy[cellGroupi, 1])

  # Get BigWigs
  getGroupBW(ArchRProj = proj, groupBy = "bioNames", tileSize = 50, maxCells = 20000, normMethod = 'ReadsInTSS')
  getGroupBW(ArchRProj = proj, groupBy = "Clusters", tileSize = 50, maxCells = 20000, normMethod = 'ReadsInTSS')
  getGroupBW(ArchRProj = proj, groupBy = "gex_mammary_gland_refine_clusters", tileSize = 50, maxCells = 20000, normMethod = 'ReadsInTSS')
  # callMacs2Manually()

  #------------------------------------------- Call Peaks -------------------------------------------
  DT <- data.table(as.data.frame(proj@cellColData), keep.rownames = 'cell')
  table(DT[,list(bioNames, Clusters)])
  DT[,.N, by=Clusters][order(N)]

  # this is their implementation of Iterative Overlap
  # https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html
  proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Clusters",
    peakMethod = 'Macs2',
    reproducibility = "2", # Requires 2 samples to have peak call at this locus (chose 2 based on how these cluster)
    peaksPerCell = 500,
    maxPeaks = 150000,
    minCells = 25, # also used 250 here and same result
    shift = -75,
    extsize = 150,
    threads = 4,
    force = T,
    pathToMacs2 = pathToMacs2
  )
  proj@cellColData$Clusters2 <- proj@cellColData$Clusters

  #------------------------------------------- Peak metrics -------------------------------------------
  ## Peaks were called using MACS2 above in each cluster
  proj <- addPeakMatrix(proj)
  PlotPeakAnnotations(ArchRProj = proj)

  #---- FRIPS ----#
  dt <- as.data.table(proj@cellColData[,c("bioNames","Clusters","FRIP","nFrags")])
  dt[,mean(FRIP), by=bioNames]
  cairo_pdf(filename = file.path(plot.dir, 'Multiome_FRIPS_sample.pdf'), width = 6, height = 6, bg = 'transparent')
  ggplot(dt, aes(x=bioNames, y=FRIP, fill=bioNames)) +
    geom_violin(draw_quantiles = T, trim = F, scale = 'width') +
    geom_boxplot(fill='white', alpha=.6) +
    scale_fill_manual("", values = bionames_palette_list) +
    theme_ArchR(baseSize = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), plot.background = element_blank(), legend.position = 'none')
  dev.off()

  cairo_pdf(filename = file.path(plot.dir, 'Multiome_FRIPS_cluster.pdf'), width = 7, height = 6, bg = 'transparent')
  ggplot(dt, aes(x=Clusters, y=FRIP, fill=Clusters)) +
    geom_violin(draw_quantiles = T, trim = F, scale = 'width') +
    geom_boxplot(fill='white', alpha=.6) +
    scale_fill_manual("", values = cluster_palette) +
    theme_ArchR(baseSize = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), plot.background = element_blank(), legend.position = 'none')
  dev.off()

  saveArchRProject(ArchRProj = proj, load = TRUE)

  getAvailableMatrices(proj)

  #------------------------------------------- Accessible Regions -------------------------------------------
  #### Sample Accessible regions ####
  # Identify marker peaks Per BioName
  # Performs a wilcoxon test of each group vs all others
  #proj@cellColData$bioNames <- as.character(proj@cellColData$bioNames) # otherwise this procedure ends prematurely
  sampleMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "bioNames",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  # significant peaks
  marker.sample.cutOff
  sampleMarkerList <- getMarkers(sampleMarkersPeaks, cutOff = marker.sample.cutOff) # used in plotbrowser track

  #### Write bed file of sample peaks ####
  library(rtracklayer)
  table(names(proj@peakSet))
  export.bed(proj@peakSet, con = file.path(archr.dir, 'Annotations', 'All_Peaks.bed'))

  # Per bionames
  (outDir <- file.path(getOutputDirectory(proj),'Annotations'))
  for(i in names(sampleMarkerList)){
    if(nrow(sampleMarkerList[[i]])==0)
      next
    gr <- GRanges(sampleMarkerList[[i]])
    short.name <- gsub("\\(|\\/|\\)","",i)
    gr@ranges@NAMES <- short.name
    export.bed(gr, con = file.path(outDir, paste0(short.name,'.bed')))
  }

  #--------------------------------------- Plots Peaks in Browser ---------------------------------------
  #### Browse peaks ####
  ArchRBrowser(proj)

  p <- plotBrowserTrack(
    ArchRProj = proj,
    #region = GRanges(markerList$`C1-Luminal Mature`[3,]),
    groupBy = "Clusters",
    #geneSymbol = c("Foxa1","Trp63","Cd55","Ctcf","Krt19","Krt14","Krt18","Pgr","Gata3","Esr1"),
    geneSymbol = sort(unique(unname(unlist(markers_list)))),
    title = '(+- 20Kb)',
    #plotSummary = c("scTrack","bulkTrack","featureTrack","geneTrack"),
    #sizes = c(10, 10, 2, 2.5),
    plotSummary = c("bulkTrack","featureTrack","geneTrack"),
    sizes = c(10, 2, 2.5),
    # features =  getMarkers(clusterMarkersPeaks, cutOff = cutoff, returnGR = TRUE), # We can restrict to significant peaks or use all peaks
    upstream = 20000,
    downstream = 20000
  )
  grid.draw(p$Foxa1)
  grid.draw(p$Ar)
  grid.draw(p$Krt5)
  grid.draw(p$Krt8)

  #plotPDF(p, name = "Plot-Tracks-With-Features_Groups", width = 6, height = 12, ArchRProj = proj, addDOC = FALSE)
  plotPDF(p, name = "Plot-Tracks-With-Features_Groups", width = 6, height = 8, ArchRProj = proj, addDOC = FALSE)

  gmex <- ArchR:::.getFeatureDF(ArrowFiles = getArrowFiles(proj), subGroup = "GeneExpressionMatrix")
  p <- plotBrowserTrack(
    ArchRProj = proj,
    groupBy = "Clusters",
    geneSymbol = intersect(sort(unique(unname(unlist(markers_list)))), gmex$name),
    title = '(+- 2Kb)',
    useMatrix = 'GeneExpressionMatrix',
    log2Norm = F,
    plotSummary = c("bulkTrack","featureTrack","geneTrack"),
    sizes = c(10, 2, 2.5),
    upstream = 2000,
    downstream = 2000
  )
  grid.draw(p$Foxa1)
  grid.draw(p$Krt5)
  grid.draw(p$Krt8)
  plotPDF(p, name = "Plot-Expression_Tracks-With-Features_Groups", width = 6, height = 8, ArchRProj = proj, addDOC = FALSE)

  #### Disjoint bed and narrowPeak files of marker peaks ####
  mammaryRefineMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "gex_mammary_gland_refine_clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  marker.cutOff
  mammaryRefineMarkersList <- getMarkers(mammaryRefineMarkersPeaks, cutOff = marker.cutOff) # used in plotbrowser track
  # Export as narrowPeak
  # https://genome.ucsc.edu/FAQ/FAQformat.html#format12
  # chr, start, end, name, score, strand, signalvalue, pvalue, qvalue, peak(an offset from start)
  (annotation.dir <- file.path(getOutputDirectory(proj),'Annotations'))
  mammaryRefineMarkersList <- mammaryRefineMarkersList[names(sort(lengths(mammaryRefineMarkersList), decreasing = T))]
  # add cluster name to row
  mammaryRefineMarkersListDisjoint <- lapply( names(mammaryRefineMarkersList), function(x){
    mammaryRefineMarkersList[[x]]$Cluster <- x
    sort(mammaryRefineMarkersList[[x]])
  })
  mammaryRefineMarkersListDisjoint <- unlist(as(mammaryRefineMarkersListDisjoint, "GRangesList"), use.names = F)
  table(mammaryRefineMarkersListDisjoint$Cluster)
  ## sort by location and FDR, then get unique
  mammaryRefineMarkersListDisjoint <- sort(mammaryRefineMarkersListDisjoint, by=~seqnames+start+end+FDR) %>% unique(.)
  table(mammaryRefineMarkersListDisjoint$Cluster)

  (uqclusters <- unique(mammaryRefineMarkersListDisjoint$Cluster))
  lapply(uqclusters, function(i){
    short.name <- gsub(" ","_",i)
    sigLabel <- gsub("[<=&>=]","",marker.cutOff) %>% gsub("[ ]+","_",.)
    dt <- as.data.table(mammaryRefineMarkersListDisjoint[mammaryRefineMarkersListDisjoint$Cluster==i])
    ## Write Bed file
    fout <- file.path(annotation.dir, paste0(short.name,'_',sigLabel,'_disjoint.bed'))
    message("Writing file:",fout)
    write.table(dt[, list(seqnames, start, end)][order(seqnames, start, end)],
                file = fout, sep = "\t", row.names = F, col.names = F, quote = F)

    fout <- file.path(annotation.dir,paste0(short.name,'_',sigLabel,'_disjoint.narrowPeak'))
    message("Writing to ",fout)
    narrowpeak <- dt[, list(chr=seqnames, start, end, name=".",
                            score=Log2FC, strand='*', signalValue='.', pvalue='.',
                            qvalue=FDR, peak=round((end-start)/2)) ][order(chr, start, end)]
    write.table(narrowpeak, file = fout, sep = "\t", row.names = F, col.names = F, quote = F)
  })

  clustersMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  marker.cutOff
  clustersMarkersPeaksList <- getMarkers(clustersMarkersPeaks, cutOff = marker.cutOff) # used in plotbrowser track
  # Export as narrowPeak
  # https://genome.ucsc.edu/FAQ/FAQformat.html#format12
  # chr, start, end, name, score, strand, signalvalue, pvalue, qvalue, peak(an offset from start)
  (annotation.dir <- file.path(getOutputDirectory(proj),'Annotations'))
  clustersMarkersPeaksList <- clustersMarkersPeaksList[names(sort(lengths(clustersMarkersPeaksList), decreasing = T))]
  clustersMarkersPeaksList <- clustersMarkersPeaksList[which(lengths(clustersMarkersPeaksList)>0)]
  # add cluster name to row
  clustersMarkersPeaksListDisjoint <- lapply( names(clustersMarkersPeaksList), function(x){
    clustersMarkersPeaksList[[x]]$Cluster <- x
    sort(clustersMarkersPeaksList[[x]])
  })
  clustersMarkersPeaksListDisjoint <- unlist(as(clustersMarkersPeaksListDisjoint, "GRangesList"), use.names = F)
  table(clustersMarkersPeaksListDisjoint$Cluster)
  ## sort by location and FDR, then get unique
  clustersMarkersPeaksListDisjoint <- sort(clustersMarkersPeaksListDisjoint, by=~seqnames+start+end+FDR) %>% unique(.)
  table(clustersMarkersPeaksListDisjoint$Cluster)

  (uqclusters <- unique(clustersMarkersPeaksListDisjoint$Cluster))
  lapply(uqclusters, function(i){
    short.name <- gsub(" ","_",i)
    sigLabel <- gsub("[<=&>=]","",marker.cutOff) %>% gsub("[ ]+","_",.)
    dt <- as.data.table(clustersMarkersPeaksListDisjoint[clustersMarkersPeaksListDisjoint$Cluster==i])
    ## Write Bed file
    fout <- file.path(annotation.dir, paste0(short.name,'_',sigLabel,'_disjoint.bed'))
    message("Writing file:",fout)
    write.table(dt[, list(seqnames, start, end)][order(seqnames, start, end)],
                file = fout, sep = "\t", row.names = F, col.names = F, quote = F)

    fout <- file.path(annotation.dir,paste0(short.name,'_',sigLabel,'_disjoint.narrowPeak'))
    message("Writing to ",fout)
    narrowpeak <- dt[, list(chr=seqnames, start, end, name=".",
                            score=Log2FC, strand='*', signalValue='.', pvalue='.',
                            qvalue=FDR, peak=round((end-start)/2)) ][order(chr, start, end)]
    write.table(narrowpeak, file = fout, sep = "\t", row.names = F, col.names = F, quote = F)
  })

  #### Write fragment files per group ####
  FragmentFilesFromArrow(proj = proj, group.by = NA)

  FragmentFilesFromArrow(proj = proj)

  FragmentFilesFromArrow(proj = proj,
                         out.dir = file.path(data.dir,"fragments_by_gland_refine"),
                         group.by = "gex_mammary_gland_refine_clusters")

  #------------------------------------------- Gene Scores -------------------------------------------
  # 1. Accessibility within the entire gene body contributes to the gene score.
  # 2. An exponential weighting function that accounts for the activity of putative distal regulatory elements in a distance-dependent fashion.
  # 3. Imposed gene boundaries that minimizes the contribution of unrelated regulatory elements to the gene score.
  # https://www.archrproject.com/bookdown/identification-of-marker-features.html
  # identify features (here gene scores) specific to groupby=clusters.
  markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"), # parameters for background generation
    testMethod = "wilcoxon"
  )

  cat("using cutoff:",gene.marker.sample.cutOff)
  markerList <- getMarkers(markersGS, cutOff = gene.marker.sample.cutOff)
  heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = gene.marker.sample.cutOff,
    labelMarkers = unique(unlist(sapply(markerList, function(x){ x$name %>% head(10)} ))),
    transpose = T
  )
  heatmapGS
  plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-Top10-per-Cluster", width = 12, height = 6, ArchRProj = proj, addDOC = FALSE)

  heatmapGS_clusters <- lapply(names(markerList), function(x){
    genes2lab <- markerList[[x]]$name %>% head(50)
    heatmapGS <- plotMarkerHeatmap(
      seMarker = markersGS,
      cutOff = gene.marker.sample.cutOff,
      labelMarkers = genes2lab,
      transpose = T
    )
    if('column_title' %in% slotNames(heatmapGS))
      heatmapGS@column_title <- paste0("Top 50 Cluster ",x)
    heatmapGS
  })

  plotPDF(heatmapGS_clusters, name = "GeneScores-Marker-Heatmap-Top50-per-Cluster", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

  # UMAP marker genes no imputation
  markerGenes2Label  <- unique(unlist(sapply(markerList, function(x){ head(x$name,20) } )))
  markerGenes2Label <- sort(unique(c(markerGenes2Label, unlist(markers_list))))
  ptop <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes2Label,
    plotAs = 'points',
    embedding = "UMAP",
    size = 0.5,
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
  )

  #library(doParallel)
  dir.create(path = file.path(plot.dir,'marker_gene_scores_no_impute'), showWarnings = F)
  foreach(i = names(ptop)) %dopar%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = ptop[[i]],
                    addLabels = F,
                    pt.size = 1,
                    title = paste(i,"gene scores"),
                    legend.title = expression(Log[2]~('norm counts'+1)),
                    filename.prefix = paste0("/marker_gene_scores_no_impute/UMAP-ATAC-Marker-Genes-No-Imputation-",i))

  # UMAP Imputed marker genes
  proj <- addImputeWeights(proj)
  # UMAP marker genes no imputation
  addi <- c('Arid5a','Bach1','Foxa1','Foxb2','Foxi1','Foxl1','Hif3a','Nfkb1','Spdef','Trp63')
  markerGenes2Label  <- unique(unlist(sapply(markerList, function(x){ head(x$name,20) } )))
  markerGenes2Label <- c(markerGenes2Label, addi)
  markerGenes2Label <- sort(unique(c(markerGenes2Label, unlist(markers_list))))
  ptop <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes2Label,
    plotAs = 'points',
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(proj)
  )

  dir.create(path = file.path(plot.dir,'marker_gene_scores_magic_impute'), showWarnings = F)
  foreach(i = names(ptop)) %dopar%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = ptop[[i]],
                    addLabels = F,
                    pt.size = 1,
                    title = paste("UMAP of InterativeLSI colored by\nGeneScoreMatrix:",i),
                    legend.title = expression(Log[2]~('norm counts'+1)),
                    filename.prefix = paste0("/marker_gene_scores_magic_impute/UMAP-ATAC-Marker-Genes-Magic-Imputation-",i))

  #### Genome tracks ####
  p <- plotBrowserTrack(
    ArchRProj = proj,
    groupBy = "Clusters",
    geneSymbol = c('Arid5a','Bach1','Foxa1','Foxb2','Foxi1','Foxl1','Hif3a','Nfkb1','Spdef','Trp63'),
    upstream = 50000,
    downstream = 50000
  )

  plotPDF(plotList = p,
          name = "Plot-Tracks-Marker-Genes.pdf",
          ArchRProj = proj,
          addDOC = FALSE, width = 5, height = 5)

  #------------------------------------------- scRNA integrate multiomic -------------------------------------------
  proj <- loadArchRProject(path = archr.dir)
  library(SingleCellExperiment)
  library(Seurat)

  # pbmc[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA.
  # To assist with visualization and interpretation. we also convert Pearson residuals back to ‘corrected’ UMI counts.
  # You can interpret these as the UMI counts we would expect to observe if all cells were sequenced to the same depth. If you want to see exactly how we do this, please look at the correct function here.
  # The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@counts.
  # We store log-normalized versions of these corrected counts in pbmc[["SCT"]]@data, which are very helpful for visualization.

  arid1a.multiome <- readRDS(file = multiRNA.file)

  png(filename = file.path(plot.dir,"Dist_test.png"), width = 5, height = 10, units = 'in', res = 150)
  par(mfrow=c(2,1))
  plot(density(arid1a.multiome@assays$RNA@data["Pou2f1",]), main="RNA data")
  plot(density(arid1a.multiome@assays$SCT@data["Pou2f1",]), main="SCT data")
  dev.off()

  adata <- SingleCellExperiment(assays = list(#counts = GetAssayData(object = arid1a.multiome, slot = "counts"),
    #data = GetAssayData(object = arid1a.multiome, slot = "data"), # log1p transformed normalized counts
    #normcounts = GetAssayData(object = arid1a.multiome, slot = "scale.data"),
    counts = GetAssayData(object = arid1a.multiome, assay='RNA', slot = "counts"),
    rnadata = GetAssayData(object = arid1a.multiome, assay='RNA', slot = "data"),
    sctdata = GetAssayData(object = arid1a.multiome, assay='SCT', slot = "data")),
    colData=DataFrame(arid1a.multiome@meta.data),
    #rowData=DataFrame(arid1a.multiome@assays$RNA@meta.features),
    reducedDims = arid1a.multiome@reductions)
  head(colnames(adata))
  table(adata@colData$seurat_clusters)

  table(gsub("_multiome.*","",rownames(proj@cellColData)))
  table(gsub("_multiome.*","",colnames(adata)))
  table(gsub("_multiome.*","", setdiff(colnames(adata), rownames(proj@cellColData)) ))

  message("# cells in common:", length(intersect(colnames(adata), rownames(proj@cellColData))))
  message("# cells in RNA not in ATAC:", length(setdiff(colnames(adata), rownames(proj@cellColData))))
  message("# cells in ATAC not in RNA:", length(setdiff(rownames(proj@cellColData), colnames(adata))))

  unique(gsub("#.*","",rownames(proj@cellColData)))
  unique(gsub("#.*","",colnames(adata)))

  #expand the rna matrix to contain 0 for missing gene entries
  mat <- Matrix::Matrix(assay(adata, 'rnadata'), sparse = TRUE) # log1p transformed normalized counts
  #mat <- Matrix::Matrix(assay(adata, 'sctdata'), sparse = TRUE) # log1p transformed normalized counts
  dim(mat)
  (missingcells <- setdiff(rownames(proj@cellColData), colnames(mat)))
  if(length(missingcells) > 0){
    message("There are missing cells!")
    missingmat <- Matrix::Matrix(data = 0,
                                 nrow = nrow(mat),
                                 ncol = length(missingcells),
                                 sparse = T, dimnames = list(rownames(adata), missingcells))
    mat <- cbind(mat, missingmat)
  }
  dim(mat)

  #expand the colData
  dtfound <- data.table(arid1a.multiome@meta.data, keep.rownames = 'cellid')
  dtrna <- dtfound
  dtmiss <- data.table(orig.ident=gsub("#.*","",missingcells), cellid = missingcells)
  if(nrow(dtmiss))
    dtrna <- rbindlist(list(dtmiss, dtfound), fill = T)

  # Integrate ATAC and sce
  addGeneExpressionMatrixUpdate(input = proj,
                                seRNA = SingleCellExperiment(list(counts=mat), colData=dtrna),
                                scaleTo = NULL,
                                normalize = F,
                                gtfFile = file.path(annot.dir, 'genes.gtf'),
                                matrixName = "GeneExpressionMatrix",
                                force = T)

  #subset for existing atacdata
  dtupdate <- dtrna[cellid %in% rownames(proj@cellColData)]
  dtupdate
  # sort it
  dtupdate <- dtupdate[order(match(cellid, rownames(proj@cellColData)))]
  setnames(dtupdate, paste0("gex_",colnames(dtupdate)))
  input <- proj
  if(identical(dtupdate$gex_cellid, rownames(input@cellColData))){
    for (i in colnames(dtupdate)){
      #message("Checking column ",i)
      if(!i %in% colnames(input@cellColData)){
        message("Adding column ",i)
        input <- addCellColData(ArchRProj = input, data = as.vector(unlist(dtupdate[,i,with=F], use.names = F)), name = i, cells = paste0(dtupdate$gex_cellid))
      }
    }
  }
  proj <- input

  # check a few entries
  # check a few entries
  getAvailableMatrices(proj)

  gp <- lapply(c("Gata3","Pou2f1"), function(g){
    UmapPlotEnhance(ArchRProj = proj,
                    colorBy = "GeneExpressionMatrix",
                    outline.color = "gray80",
                    feature = g,
                    title = g,
                    addLabels = F,
                    scale_color_bar = F,
                    break.unit = .5,
                    vmin = 0,
                    legend.title = "Normalized expression",
                    panel.border.size = 0,
                    legend.height.factor = 0.6,
                    pal = brewer.bupu(15))
  })
  png(filename = file.path(plot.dir,"Gata3_Pou2f1_ATAC_exp_umap.png"), width = 10, height = 5, units = 'in', res = 150)
  do.call(grid.arrange, c(gp, nrow=1))
  dev.off()

  ## Check import vs Seurat expression
  expassay <- getMatrixFromProject(ArchRProj = proj, useMatrix = 'GeneExpressionMatrix', binarize = F, threads = 2)
  range(expassay@assays@data$GeneExpressionMatrix)
  range(mat)

  for(g in c("Foxa1","Krt5","Krt8","Hif3a")){
    a <- range(arid1a.multiome@assays$RNA@data[g,])
    b <- range(expassay[which(rowData(expassay)$name==g),]@assays@data$GeneExpressionMatrix)
    message("Range in seurat ",g," ",paste(a,collapse = ":"))
    message("Range in arch-R ",g," ",paste(a,collapse = ":"))
  }

  if( length(setdiff(rownames(proj), rownames(expassay@colData))) > 0 )
    warnings("It appears several cells are missing and this might be a problem later to plot dim reductions!")

  rm(list = c('missingmat', 'mat', 'arid1a.multiome', 'expassay'))
  gc()

  saveArchRProject(ArchRProj = proj)

  proj <- loadArchRProject(path = archr.dir)

  #### Refine Saeki classes ####
  proj$gex_mammary_gland_refine_clusters <- proj$gex_mammary_gland_cluster_classification
  proj@cellColData[proj@cellColData$gex_bioNames %in% c("Arid1a(-/-)","Arid1a(-/-) E2") &
                     proj@cellColData$gex_mammary_gland_cluster_classification=="L-Hor",]$gex_mammary_gland_refine_clusters <- "L-Hor-Arid1a-"
  #proj@cellColData[proj@cellColData$gex_mammary_gland_cluster_classification=="Basal/MaSC",]$gex_mammary_gland_refine_clusters <- "Basal_MaSC"
  proj@cellColData[proj@cellColData$gex_mammary_gland_refine_clusters=="Basal_MaSC",]$gex_mammary_gland_refine_clusters <- "Basal/MaSC"
  table(proj$gex_mammary_gland_refine_clusters)

  #------------------------------------------- Plot scRNA integration -------------------------------------------

  PlotIntegrationHeatmap_ATAC_RNA(proj, show.percentages = T)
  ## % bar plots (Integration heatmap is better)

  p1 <- plotEmbedding(
    proj,
    baseSize = 12,
    size = 1.1,
    shape = 21,
    stroke = .05,
    color = 'black',
    plotAs = "points",
    colorBy = "cellColData",
    name = c("gex_seurat_clusters","gex_bioNames"),
    imputeWeights = NULL,
    pal = rna_cluster_palette,
    quantCut = c(0,1),
    threads = 4
  )

  # Plot the integration UMAP
  p1$gex_seurat_clusters$data$color <- as.integer(gsub("-.*","",p1$gex_seurat_clusters$data$color))
  p1$gex_seurat_clusters$data <- p1$gex_seurat_clusters$data[p1$gex_seurat_clusters$data$color > 0,]
  p1$gex_seurat_clusters$data <- p1$gex_seurat_clusters$data[order(p1$gex_seurat_clusters$data$color),]
  p1$gex_seurat_clusters$data$color <- factor(x = p1$gex_seurat_clusters$data$color,
                                              levels = unique(p1$gex_seurat_clusters$data$color), ordered = T)
  UmapPlotEnhance(ArchRProj = proj,
                  p1 = p1$gex_seurat_clusters,
                  addLabels = F,
                  title.size = 18,
                  title = 'scATAC colored by RNA clusters',
                  show.legend = T,
                  panel_width_in = 4,
                  scale_color_bar = F,
                  pal = rna_cluster_palette,
                  filename.prefix = 'UMAP_RNA_Integration_cell_clusters')

  p1$gex_bioNames$data <- p1$gex_bioNames$data[order(p1$gex_bioNames$data$color),]
  table(p1$gex_bioNames$data$color)
  p1$gex_bioNames$data$color <- gsub("^[0-9]+-","",p1$gex_bioNames$data$color)
  table(p1$gex_bioNames$data$color)
  bionames_palette_list
  p1$gex_bioNames$data$color <- factor(x = p1$gex_bioNames$data$color, levels = unique(p1$gex_bioNames$data$color), ordered = T)
  UmapPlotEnhance(ArchRProj = proj,
                  p1 = p1$gex_bioNames,
                  addLabels = F,
                  title = 'scATAC colored by RNA samples',
                  pal=bionames_palette_list,
                  show.legend = T,
                  filename.prefix = 'Plot-UMAP-RNA-Integration-cells')
  ## UMAP per cluster
  table(p1$gex_seurat_clusters$data$color)
  cplots <- lapply(sort(unique(p1$gex_seurat_clusters$data$color)), function(i){
    dt <- data.table(p1$gex_seurat_clusters$data)
    dt[,clust := as.integer(gsub("^C","",gsub("-.*","",color)))]
    dt <- dt[order(clust),]
    dtI <- dt[color==i]
    dtnI <- dt[color!=i]
    dt <- rbindlist(list(dtnI, dtI))
    dt[,color := factor(x = color, levels = unique(color), ordered = T)]

    median.labels <- dt[,.(x=median(x),y=median(y)), by=color]
    median.labels <- median.labels[order(color)]
    median.labels[color=='Luminal Mature', x:=x-4]

    pal.highlight <- rna_cluster_palette
    pal.highlight[names(pal.highlight) %ni% i] <- NA

    ##
    p <- ggplot(dt, aes(x=x, y=y, fill=color), color='gray80') +
      geom_point(size=1.5, shape=21, stroke=.05) + # deintensified stroke
      theme_ArchR() +
      scale_fill_manual("", values=pal.highlight, na.value = 'gray80') +
      xlab("") + ylab("") +
      coord_fixed(clip = 'off') +
      theme(#legend.direction = "vertical",
        legend.position = 'none',
        text = element_text(size = 15),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 12),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(x = c(0,0,0,.5), units = 'cm'))
    p <- p + geom_shadowtext(data = median.labels[color==i], aes(x=x, y=y, label=color), color='white', size=8)
  })

  ggsave(filename = file.path(plot.dir, 'Plot-UMAP-RNA-Integration_Clusters.png'),
         plot = do.call(grid.arrange, c(cplots, nrow=4)),
         width = 11, height = 10, units = 'in',
         device = 'png',
         bg = 'transparent',
         dpi = 200)

  #### Plot UMAP of expression ####
  markerRNA <- c("Smarca4")
  p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneExpressionMatrix",
    name = sort(markerRNA),
    embedding = "UMAP",
    pal = solarHalf_palette,
    plotAs = 'points',
    log2Norm = F,
    size = 0.9,
    imputeWeights = NULL #getImputeWeights(proj)
  )
  UmapPlotEnhance(ArchRProj = proj, p1 = p,
                  title = gsub(".*: ","",p$labels$title),
                  legend.title = 'scRNA normalized',
                  pal = solarHalf_palette,
                  title.size = 14,
                  scale_color_bar = F,
                  as.pdf = T,
                  filename.prefix = paste0("UMAP_",gsub(".*: ","",p$labels$title)))
  p$Smarca4
  pNolog <- copy(p)
  grid.arrange(p$Gata3, pNolog$Gata3)
  do.call(grid.arrange, c(p, nrow=3))

  pENL <- lapply(pNolog, function(x){
    UmapPlotEnhance(ArchRProj = proj, p1 = x, title = x$labels$title,
                    legend.title = 'scRNA normalized',
                    pal = solarHalf_palette, title.size = 14, scale_color_bar = F)
  })

  pE <- lapply(p, function(x){
    UmapPlotEnhance(ArchRProj = proj, p1 = x, title = x$labels$title,
                    legend.title = 'scRNA log2(normalized)',
                    pal = solarHalf_palette, title.size = 14, scale_color_bar = F)
  })

  grid.arrange(pENL$Gata3, pE$Gata3)
  message("Prefer to plot normalized over log2(normalized)")

  zscans <- grep("Zscan", getFeatures(ArchRProj = proj, useMatrix = "GeneExpressionMatrix"), value = T)
  p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneExpressionMatrix",
    name = zscans,
    embedding = "UMAP",
    plotAs = 'points',
    size = 0.5,
    imputeWeights = NULL #getImputeWeights(proj)
  )
  p1 <- lapply(p, function(x){
    UmapPlotEnhance(ArchRProj = proj, p1 = x, title = x$labels$title,
                    legend.title = 'scRNA expression', title.size = 14,
                    scale_color_bar = F, pal = rocket_palette)
  })
  png(filename = file.path(plot.dir,'UMAP-scRNA_zscan.png'), width = 24, height = 10,
      units = 'in', bg = 'transparent', res = 150)
  do.call(grid.arrange, c(p1, nrow=2, top="scRNA expression per cell (no imputation)"))
  dev.off()

  #### Plot single expression ####
  genes <- getFeatures(proj, useMatrix = "GeneExpressionMatrix")
  sort(grep("Elk4", genes, ignore.case = T, value = T))
  markerRNA <- c("Elk4")
  UmapPlotEnhance(ArchRProj = proj,
                  colorBy = "GeneExpressionMatrix",
                  feature = markerRNA,
                  title = markerRNA,
                  addLabels = F,
                  scale_color_bar = F,
                  break.unit = 1,
                  legend.title = "Normalized expression",
                  panel.border.size = .2,
                  pal = solarHalf_palette)

  markerRNA <- c("Trp63")
  p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneExpressionMatrix",
    name = sort(markerRNA),
    embedding = "UMAP",
    plotAs = 'points',
    size = 0.5,
    imputeWeights = NULL #getImputeWeights(proj)
  )
  UmapPlotEnhance(ArchRProj = proj,
                  p1 = p,
                  title = paste0("Expression of ",markerRNA),
                  legend.title = 'scRNA expression',
                  title.size = 14,
                  panel.border.size = 0,
                  #pt.size = 2,
                  pal = solarHalf_palette, # brewer.bupu(20),
                  scale_color_bar = F,
                  #outline.size = .2,
                  filename.prefix = paste0(markerRNA,'_expression'))

  #------------------------------------------- Plot mammary cell type umaps -------------------------------------------
  p <- plotEmbedding(ArchRProj = proj,
                     colorBy = "cellColData",
                     name = "gex_mammary_gland_cluster_classification",
                     embedding = "UMAP",
                     plotAs = 'points') %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    pt.size = 1,
                    panel.border.size = 0,
                    addLabels = F,
                    #title = 'Luminal / Basal',
                    label.size = 14,
                    title.size = 20,
                    break.unit = .1,
                    scale_color_bar = F,
                    pal = mammary_cluster_class_palette, #c('#1f77b4','gray','#FFD966'),
                    legend.title = '',
                    filename = 'ATAC_gex_mammary_gland_cluster_classification_umap')

  #------------------------------------------- mammary cell type umaps with label -------------------------------------------
  plabel <- plotEmbedding(ArchRProj = proj,
                            colorBy = "cellColData",
                            name = "Clusters",
                            embedding = "UMAP",
                            plotAs = 'points')
  # smaller point size to show cells
  plabel <- UmapPlotEnhance(ArchRProj = proj,
                            p1 = plabel,
                            pt.size = .6,
                            panel.border.size = 0,
                            addLabels = T,
                            #title = 'Luminal / Basal',
                            label.size = 8,
                            title.size = 20,
                            break.unit = .1,
                            scale_color_bar = F,
                            pal = NULL,
                            legend.title = '')

  p <- plotEmbedding(ArchRProj = proj,
                     colorBy = "cellColData",
                     name = "gex_mammary_gland_cluster_classification",
                     embedding = "UMAP",
                     plotAs = 'points') %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    pt.size = 1,
                    panel.border.size = 0,
                    addLabels = F,
                    label.size = 14,
                    title.size = 20,
                    break.unit = .1,
                    scale_color_bar = F,
                    pal = mammary_cluster_class_palette,
                    legend.title = '')
  png(filename = file.path(plot.dir,'ATAC_gex_mammary_gland_cluster_classification_umap.png'), width = 5.5, height = 4, units = 'in', res = 150, bg = "transparent")
  #pdf(file = file.path(plot.dir,'ATAC_luminal_basal_umap_label.pdf'), width = 5, height = 5, bg = "transparent")
  p + plabel$layers
  dev.off()

  #### Export to Scanpy ####
  ## Export for Scanpy
  ArchRProjectToSeurat(build = T,
                       ArchRProj = proj,
                       matrix = 'VierstraMatrix',
                       file.path(robj.dir,'arid1a_archR_atac_peaks_multiome_sce.rds'))

  #------------------------------------------- Ar targets Umaps -------------------------------------------
  p <- myplotEmbedding(ArchRProj = proj,
                       colorBy = "GeneExpressionMatrix",
                       name = "Ar",
                       embedding = "UMAP",
                       plotAs = 'points',
                       size=0.2, baseSize = 14)

  UmapPlotEnhance(ArchRProj = proj,
                  p1 = p,
                  pt.size = 1,
                  panel.border.size = 0,
                  addLabels = F,
                  title = 'Ar',
                  label.size = 14,
                  title.size = 20,
                  break.unit = 1,
                  scale_color_bar = F,
                  pal = solarHalf_palette, # as.vector(pals::ocean.balance(20)),
                  legend.title = 'Ar expression',
                  filename = 'ATAC_Ar_umap')

  p <- myplotEmbedding(ArchRProj = proj,
                       colorBy = "cellColData",
                       name = "gex_ar_targets",
                       embedding = "UMAP",
                       plotAs = 'points',
                       size=0.2, baseSize = 14)

  UmapPlotEnhance(ArchRProj = proj,
                  p1 = p,
                  pt.size = 1,
                  panel.border.size = 0,
                  addLabels = F,
                  title = 'Ar targets',
                  label.size = 14,
                  title.size = 20,
                  break.unit = .1,
                  scale_color_bar = F,
                  legend.size = 20,
                  pal = ArchRPalettes$solarExtra, # as.vector(pals::ocean.balance(20)),
                  legend.title = 'Ar expression target scores',
                  filename = 'ATAC_Artargets_umap')

  #--------------------------------------------Pseudotime--------------------------------------------------

  # Import pseudotime from Palantir
  pt <- fread(file.path(data.dir, 'scanpy', 'scATAC_pseuodtime.txt'))
  setnames(pt, "V1","cells")
  setequal(proj$cellNames, pt$cells)
  pt <- pt[order(match(cells, proj$cellNames))]
  if(identical(proj$cellNames, pt$cells)){
    proj$palantir_pseudotime <- pt$palantir_pseudotime
    proj$palantir_entropy <- pt$palantir_entropy
    proj$`palantir_fate_probabilities_L_Hor_Arid1a` <- pt$`palantir_fate_probabilities_L-Hor-Arid1a-`
    proj$`palantir_fate_probabilities_L_Hor` <- pt$`palantir_fate_probabilities_L-Hor`
    proj$`palantir_fate_probabilities_L_Alv` <- pt$`palantir_fate_probabilities_L-Alv`
  }

  dt <- as.data.table(proj@cellColData[,c("bioNames","Clusters","gex_mammary_gland_refine_clusters","palantir_pseudotime")])
  dt[,mean(palantir_pseudotime), by=gex_mammary_gland_refine_clusters]
  dt[,gex_mammary_gland_refine_clusters := factor(x = gex_mammary_gland_refine_clusters,
                                                  levels = c('L-Alv','LA-Pro','Basal/MaSC','L-Hor','L-Hor-Arid1a-')) ]
  cairo_pdf(filename = file.path(plot.dir, 'Multiome_Pseudotime_Violin.pdf'), width = 4.5, height = 5, bg = 'transparent')
  ggplot(dt, aes(x=gex_mammary_gland_refine_clusters, y=palantir_pseudotime, fill=gex_mammary_gland_refine_clusters)) +
    geom_violin(draw_quantiles = T, trim = F, scale = 'width') +
    geom_boxplot(fill='white', alpha=.6, size = .9, linewidth = .5, outliers = T) +
    #geom_beeswarm_rast(dodge.width = 1, size = .05, cex = .1, scale=.2) +
    scale_fill_manual("", values = mammary_refine_clusters_palette) +
    ylab("Pseudotime") + xlab("") +
    theme_ArchR(baseSize = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), plot.background = element_blank(), legend.position = 'none')
  dev.off()

  ppalantir_pseudotime <- plotEmbedding(ArchRProj = proj,
                                        colorBy = "cellColData",
                                        name = "palantir_pseudotime",
                                        plotAs = "points", # Need this to return dataframe
                                        embedding = "FR") %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    addLabels = F,
                    pt.size = 2,
                    outline.size = .01,
                    panel.border.size = 0,
                    legend.height.factor = .35,
                    legend.size = 18,
                    title = "Pseudotime",
                    pal = pals::cividis(100), #rev(pals::brewer.rdylbu(30)),
                    scale_color_bar = F,
                    vmax = 1,
                    title.size = 18)
  ppalantir_pseudotime

  pseudotimes <- c("palantir_pseudotime","palantir_entropy")
  gp <- lapply(pseudotimes, function(x){
    p <- plotEmbedding(ArchRProj = proj,
                       colorBy = "cellColData",
                       name = x,
                       plotAs = "points", # Need this to return dataframe
                       embedding = "FR") %>%
      UmapPlotEnhance(ArchRProj = proj,
                      p1 = .,
                      legend.height.factor = .35,
                      legend.size = 18,
                      addLabels = F,
                      pt.size = 1.5,
                      outline.size = .01,
                      panel.border.size = 0,
                      title = x,
                      pal = rev(pals::brewer.rdylbu(30)),
                      scale_color_bar = F,
                      title.size = 18)
  }) %>% set_names(pseudotimes)
  pdf(file = file.path(plot.dir,"pseudotimes.pdf"), width = 10, height = 5)
  grid.arrange(gp$palantir_pseudotime,
               gp$palantir_entropy, nrow=1)
  dev.off()

  fates <- c("palantir_fate_probabilities_L_Hor_Arid1a",
                   "palantir_fate_probabilities_L_Hor",
                   "palantir_fate_probabilities_L_Alv")
  gp <- lapply(fates, function(x){
    p <- plotEmbedding(ArchRProj = proj,
                       colorBy = "cellColData",
                       name = x,
                       plotAs = "points", # Need this to return dataframe
                       embedding = "FR") %>%
      UmapPlotEnhance(ArchRProj = proj,
                      p1 = .,
                      legend.height.factor = .35,
                      legend.size = 18,
                      addLabels = F,
                      pt.size = 1.5,
                      outline.size = .01,
                      panel.border.size = 0,
                      title = x,
                      pal = pals::plasma(30), #rev(pals::brewer.rdylbu(30)),
                      scale_color_bar = F,
                      title.size = 18)
  }) %>% set_names(fates)
  pdf(file = file.path(plot.dir,"fate_probabilities.pdf"), width = 14, height = 5)
  grid.arrange(gp$palantir_fate_probabilities_L_Hor_Arid1a,
               gp$palantir_fate_probabilities_L_Hor,
               gp$palantir_fate_probabilities_L_Alv, nrow=1)
  dev.off()

  dt <- data.table(data.frame(proj@cellColData), keep.rownames = "cells")
  dt <- dt[order(palantir_pseudotime, gex_mammary_gland_refine_clusters)]

  ggplot(dt, aes(x=palantir_pseudotime, y=gex_mammary_gland_refine_clusters, color=gex_mammary_gland_refine_clusters)) +
    geom_point() + scale_color_manual(values = mammary_refine_clusters_palette)

  #------------------------------------------- 1. Motif calling -------------------------------------------
  getAvailableMatrices(proj)

  # custom motifset (less redundant)
  pwmListNR <- NonRedundantMotifs(ArchRProj = proj,
                                  motifs = chromVARmotifs::mouse_pwms_v2,
                                  similarityMethod = 'Pearson',
                                  combine.dist.thresh = .999)
  ## My NR Matrix

  # ViestraMatrix motifset (less archetypes)
  # https://www.vierstra.org/resources/motif_clustering
  # https://www.nature.com/articles/s41586-020-2528-x

  # default motifset (redundant)
  library(TFBSTools)

  system(paste("wget https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra-Human-Motifs.rds -O", file.path(robj.dir,"Vierstra-Human-Motifs.rds")))
  mPWMs <- readRDS(file.path(robj.dir, "Vierstra-Human-Motifs.rds"))
  names(mPWMs)
  # update names
  names(mPWMs) <- gsub(":","-",names(mPWMs))
  motifPWMs <- lapply(mPWMs, function(i){
    i@name <- gsub(":","-",i@name)
    i
  })
  motifPWMs <- do.call(PWMatrixList, motifPWMs)

  vierstramotifs <- gsub("_.*","",names(motifPWMs)) %>% unique %>% toupper()
  cisbpmotifs <- gsub(".*_LINE[0-9]+_","",names(chromVARmotifs::mouse_pwms_v2)) %>% gsub("_.*","",.) %>% unique %>% toupper()
  missing <- setdiff(cisbpmotifs, vierstramotifs) %>% sort
  ## Add cisbp entries to Vierestra
  motifPWMs[['Hif3a_869#HIF-bHLH-61-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000004328_LINE324_Hif3a_I_N1']]

  motifPWMs[['Trp63_853#P53-like/1-P53L-257-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000022510_LINE6881_Trp63_I_N2']]
  motifPWMs[['Trp73_691#P53-like/1-P53L-257-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000029026_LINE1639_Trp73_D']]

  motifPWMs[['Smarcc1_843#AP1/1-BZIP-843-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000032481_LINE6526_Smarcc1_I_N2']]
  motifPWMs[['Smarcc2_646#AP1/1-BZIP-646-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000025369_LINE1534_Smarcc2_I']]

  motifPWMs[['Gata1_387#GATA-C2H2-242-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000031162_LINE1118_Gata1_D']]
  motifPWMs[['Gata2_383#GATA-C2H2-242-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000015053_LINE1112_Gata2_I']]
  motifPWMs[['Gata3_384#GATA-C2H2-242-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000015619_LINE1113_Gata3_D']]
  motifPWMs[['Gata4_386#GATA-C2H2-242-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000021944_LINE1116_Gata4_D_N1']]
  motifPWMs[['Gata5_385#GATA-C2H2-242-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000015627_LINE1114_Gata5_D']]
  motifPWMs[['Gata6_382#GATA-C2H2-242-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000005836_LINE1110_Gata6_D']]

  motifPWMs[['Foxa1_1065#FOX/4-forkhead-79-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000035451_LINE1065_Foxa1_I']]
  motifPWMs[['Foxa2_1066#FOX/4-forkhead-79-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000037025_LINE1066_Foxa2_D_N1']]
  motifPWMs[['Foxa3_3490#FOX/4-forkhead-79-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000040891_LINE3490_Foxa3_I_N47']]

  motifPWMs[['Foxj1_305#FOX/5-forkhead-80-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000034227_LINE1061_Foxj1_D']]
  motifPWMs[['Foxj2_298#FOX/5-forkhead-80-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000003154_LINE1043_Foxj2_D_N2']]
  motifPWMs[['Foxj3_303#FOX/5-forkhead-80-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000032998_LINE1051_Foxj3_D_N3']]

  motifPWMs[['Arid3a_7#HD/18-homeodomain-24-1']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000019564_LINE13_Arid3a_D_N3']]
  motifPWMs[['Arid3b_6#HD/2-homeodomain-2-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000004661_LINE10_Arid3b_D']]
  motifPWMs[['Arid3c_10#HD/2-homeodomain-2-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000066224_LINE24_Arid3c_I']]
  motifPWMs[['Arid5a_9#HD/2-homeodomain-2-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000037447_LINE18_Arid5a_D']]
  motifPWMs[['Arid5b_8#HD/2-homeodomain-2-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000019947_LINE15_Arid5b_D']]
  motifPWMs[['Arid2_11#HD/2-homeodomain-2-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000033237_LINE27_Arid2_I']]
  motifPWMs[['Arid2_11#HD/2-homeodomain-2-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000033237_LINE27_Arid2_I']]

  motifPWMs[['Tcf7_726#TCF/LEF-HMG-29-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000000782_LINE1702_Tcf7_D_N1']]
  motifPWMs[['Tcf7l1_744#TCF/LEF-HMG-29-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000055799_LINE1755_Tcf7l1_D']]
  motifPWMs[['Tcf7l2_732#TCF/LEF-HMG-29-0']] <- chromVARmotifs::mouse_pwms_v2[['ENSMUSG00000024985_LINE1715_Tcf7l2_D']]

  motifPWMs[['ZSCAN4_274#ZSCAN4-C2H2-286-2']] <- chromVARmotifs::mouse_pwms_v2[["EDL381201_LINE274_Zscan4_D_N3"]]
  motifPWMs[['ZSCAN4_353#ZSCAN4-C2H2-286-3']] <- chromVARmotifs::mouse_pwms_v2[["ENSMUSG00000023902_LINE353_Zscan10_D_N2"]]
  motifPWMs[['ZSCAN4_629#ZSCAN4-C2H2-286-4']] <- chromVARmotifs::mouse_pwms_v2[["ENSMUSG00000054272_LINE629_Zscan4c_D"]]
  motifPWMs[['ZSCAN4_719#ZSCAN4-C2H2-286-5']] <- chromVARmotifs::mouse_pwms_v2[["ENSMUSG00000061894_LINE719_Zscan20_D"]]
  motifPWMs[['ZSCAN4_795#ZSCAN4-C2H2-286-6']] <- chromVARmotifs::mouse_pwms_v2[["ENSMUSG00000070826_LINE795_Zscan4ps2_I"]]
  motifPWMs[['ZSCAN4_796#ZSCAN4-C2H2-286-7']] <- chromVARmotifs::mouse_pwms_v2[["ENSMUSG00000070828_LINE796_Zscan4f_I"]]

    saveRDS(object = motifPWMs, file = file.path(robj.dir,"Vierstra-Human-Motifs-Update.rds"))

  proj <- addMotifAnnotations(ArchRProj = proj,
                              motifSet = "cisbp",
                              cutOff = 5e-05, # p-value cutoff against bg set
                              width = 7, # motif width
                              name = "cisbp",
                              force = T)

  # Scan for all motifs (human+mouse)
  motifPWMs <- readRDS(file = file.path(robj.dir, 'Vierstra-Human-Motifs-Update.rds'))
  names(motifPWMs)
  proj <- addMotifAnnotations(proj,
                              motifPWMs = motifPWMs,
                              annoName = "Vierstra",
                              cutOff = 5e-05, # p-value cutoff against bg set
                              width = 7, # motif width
                              force = T)

  saveArchRProject(ArchRProj = proj, overwrite = T)

  getAvailableMatrices(proj)

  #### Motif 2 Gene Lookup table ####
  # Filter for motifs that are not well defined or expression is missing
  if(!'VierstraMotifToGene' %in% names(attributes(proj))){
    dtM2Glookup <- Motif2GeneLookupDT(proj = proj, motifpwm = proj@peakAnnotation@listData$Vierstra$motifs)
    dtM2Glookup$family <- gsub(".*#","", dtM2Glookup$motifID) %>% gsub("-[0-9]+-[01].*","",.)
    attr(proj, 'VierstraMotifToGene') <- dtM2Glookup
  }
  (dtM2Glookup <- proj@VierstraMotifToGene)
  fwrite(x = dtM2Glookup, file = file.path(data.dir,"Motif2GeneLookup.txt"), quote = F, sep = "\t", row.names = F, col.names = T )

  #------------------------------------------- chromVar TF per peak -------------------------------------------
  # chromVar is designed for predicting enrichment of TF activity on a per-cell basis from sparse chromatin accessibility data. The two primary outputs of chromVAR are:
  # 1. “deviations” - A deviation is a bias-corrected measurement of how far the per-cell accessibility of a given feature (i.e motif)
  # deviates from the expected accessibility based on the average of all cells or samples.
  # 2. “z-score” - The z-score, is the z-score for each bias-corrected “deviation score”  across all cells.
  #.    The absolute value of the deviation score is correlated with the per-cell read depth.
  #     This is because, with more reads, you have higher confidence that the difference in per-cell accessibility of the given feature (i.e. motif) from the expectation is greater than would occur by chance.
  ## stores both deviations and z-scores
  proj <- addDeviationsMatrix(proj, peakAnnotation = "Vierstra", force = TRUE, threads = 8)

  proj <- addDeviationsMatrix(proj, peakAnnotation = "cisbp", force = T, threads = 8)

  saveArchRProject(ArchRProj = proj, overwrite = T)

  cisbpMotifNames <- data.table(rownames(proj@peakAnnotation$cisbp$motifSummary)[
    match(names(chromVARmotifs::mouse_pwms_v2), proj@peakAnnotation$cisbp$motifSummary$ID)
  ], proj@peakAnnotation$cisbp$motifSummary$ID)

  #------------------------------------------- ChromVAR top deviations -------------------------------------------
  # chromVar is designed for predicting enrichment of TF activity on a per-cell basis from sparse chromatin accessibility data.
  # The two primary outputs of chromVAR are:
  # 1. “deviations” - A deviation is a bias-corrected measurement of how far the per-cell accessibility of a given feature (i.e motif) deviates from the expected accessibility based on the average of all cells or samples.
  # 2. “z-score” - The z-score, also known as a “deviation score” is the z-score for each bias-corrected deviation across all cells.
  #.    The absolute value of the deviation score is correlated with the per-cell read depth. This is because, with more reads, you have higher confidence that the difference in per-cell accessibility of the given feature (i.e. motif) from the expectation is greater than would occur by chance.
  ## stores both deviations and z-scores (deviation score)

  getAvailableMatrices(proj)

  #### Motif 2 Gene Lookup ####
  # Filter for motifs that are not well defined or expression is missing
  if(!'VierstraMotifToGene' %in% names(attributes(proj))){
    attr(proj, 'VierstraMotifToGene') <- Motif2GeneLookupDT(proj = proj, motifpwm = proj@peakAnnotation@listData$Vierstra$motifs)
  }
  dtM2Glookup <- proj@VierstraMotifToGene

  saveArchRProject(ArchRProj = proj, overwrite = T)

  ## Get deviations ##
  # Vierstra
  plotVierstraVarDev <- getVarDeviations(proj, name = "VierstraMatrix", plot = TRUE, n = 100)
  dtVierstraVarDev <- data.table(plotVierstraVarDev$data) %>% setnames(., 'name','motifID')
  (PlotDeviations(dt = dtVierstraVarDev, n=50, text.adjust = T)) # Plots Top deviations

  # Conclusion too much redundancy let's try to filter for expressed genes

  #### Remove some redundancy in motifs ####
  dtVierstraVarDevFilt <- merge(dtVierstraVarDev, dtM2Glookup, by='motifID')
  dtVierstraVarDevFilt[expressed==F]
  dtVierstraVarDevFilt <- dtVierstraVarDevFilt[expressed==T]
  dtVierstraVarDevFilt$family <- gsub(".*#","", dtVierstraVarDevFilt$motifID) %>% gsub("\\|[0-9]+\\|[01]$","",.)
  dtVierstraVarDevFilt[,family := gsub("-[0-9]+$","",family)]
  dtVierstraVarDevFilt <- dtVierstraVarDevFilt[order(rank, decreasing = F)]
  dtVierstraVarDevFilt[,idx := seq_along(.I)] # update idx
  dtVierstraVarDevFilt
  (PlotDeviations(dt = dtVierstraVarDevFilt, n=50, text.adjust = T)) # Plots Top deviations with expression

  writexl::write_xlsx(x = dtVierstraVarDevFilt,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra.xlsx'), col_names = T, format_headers = T)
  dtVierstraVarDevFilt <- data.table(readxl::read_xlsx(file.path(spreadsheets.dir,
                                                                 'archr_multiome_atac_motifs_vierstra.xlsx')))

  #------------------------------------------- Peak Enrichment -------------------------------------------
  # Per Sample - get enriched peaks and overlap with motifs
  sampleMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "bioNames",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )

  # Per Sample - get motifs within enriched peaks
  sampleEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = sampleMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  # Per Cluster - get enriched peaks
  proj@cellColData$Clusters <- as.character(proj@cellColData$Clusters) # otherwise this procedure ends prematurely
  clusterMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )

  # Per Cluster - get motifs within enriched peaks
  clusterEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = clusterMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.cluster.cutOff
  )

  ## By condition
  conditionMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Condition",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )

  # get motifs within enriched peaks
  conditionEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = conditionMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## By mammary cluster class
  mammary_classificationMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "gex_mammary_gland_cluster_classification",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )

  # get motifs within enriched peaks
  mammary_classificationEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = mammary_classificationMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## By mammary refine class
  mammary_refineMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "gex_mammary_gland_refine_clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )

  # get motifs within enriched peaks
  mammary_refineEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = mammary_refineMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## By cell class
  cell_classificationMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "gex_cell_classification",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )

  # get motifs within enriched peaks
  cell_classificationEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = cell_classificationMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## Pairwise C1C3
  C1C3 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C1",
    bgdGroups = "C3"
  )
  C3C1 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C3",
    bgdGroups = "C1"
  )
  # get motifs within enriched peaks
  C1_C3_MotifVierstra <- peakAnnoEnrichment(
    seMarker = cbind(C1C3, C3C1),
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## Pairwise L-Hor and L-Hor-Arid1a
  LHorLArid1a <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "gex_mammary_gland_refine_clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "L-Hor",
    bgdGroups = "L-Hor-Arid1a-"
  )
  LArid1aLHor <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "gex_mammary_gland_refine_clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "L-Hor-Arid1a-",
    bgdGroups = "L-Hor"
  )
  # get motifs within enriched peaks
  LHor_LArid1a_MotifVierstra <- peakAnnoEnrichment(
    seMarker = cbind(LHorLArid1a, LArid1aLHor),
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## Pairwise C1C2
  C1C2 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C1",
    bgdGroups = "C2"
  )
  C2C1 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C2",
    bgdGroups = "C1"
  )
  # get motifs within enriched peaks
  C1_C2_MotifVierstra <- peakAnnoEnrichment(
    seMarker = cbind(C1C2, C2C1),
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## Pairwise C2C3
  C2C3 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C2",
    bgdGroups = "C3"
  )
  C3C2 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C3",
    bgdGroups = "C2"
  )
  # get motifs within enriched peaks
  C2_C3_MotifVierstra <- peakAnnoEnrichment(
    seMarker = cbind(C2C3, C3C2),
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## Pairwise C4C5
  C4C5 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C4",
    bgdGroups = "C5"
  )
  C5C4 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C5",
    bgdGroups = "C4"
  )
  # get motifs within enriched peaks
  C4_C5_MotifVierstra <- peakAnnoEnrichment(
    seMarker = cbind(C4C5, C5C4),
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  ## GroupWise
  C1C2C3 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "C1",
    bgdGroups = c("C2","C3")
  )
  C2C3C1 <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = c("C2","C3"),
    bgdGroups = "C1"
  )
  # get motifs within enriched peaks
  C1C2C3_MotifVierstra <- peakAnnoEnrichment(
    seMarker = cbind(C1C2C3, C2C3C1),
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.sample.cutOff
  )

  # Write peaks and FDR to spreadsheets #
  dtSampleEnrichMotifVierstra <- data.table("motifID"=names(sampleEnrichMotifVierstra),
                                            "mlog10Padj"=do.call(cbind, assays(sampleEnrichMotifVierstra)$mlog10Padj),
                                            "enrichment"=do.call(cbind, assays(sampleEnrichMotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtSampleEnrichMotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_samplepeaks.xlsx'), col_names = T, format_headers = T)

  dtClusterEnrichMotifVierstra <- data.table("motifID"=names(clusterEnrichMotifVierstra),
                                             "mlog10Padj"=do.call(cbind, assays(clusterEnrichMotifVierstra)$mlog10Padj),
                                             "enrichment"=do.call(cbind, assays(clusterEnrichMotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtClusterEnrichMotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_clusterpeaks.xlsx'), col_names = T, format_headers = T)

  dtConditionEnrichMotifVierstra <- data.table("motifID"=names(conditionEnrichMotifVierstra),
                                               "mlog10Padj"=do.call(cbind, assays(conditionEnrichMotifVierstra)$mlog10Padj),
                                               "enrichment"=do.call(cbind, assays(conditionEnrichMotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtConditionEnrichMotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_condition.xlsx'), col_names = T, format_headers = T)

  dtMammaryClassificationEnrichMotifVierstra <- data.table("motifID"=names(mammary_classificationEnrichMotifVierstra),
                                                        "mlog10Padj"=do.call(cbind, assays(mammary_classificationEnrichMotifVierstra)$mlog10Padj),
                                                        "enrichment"=do.call(cbind, assays(mammary_classificationEnrichMotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtMammaryClassificationEnrichMotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_mammary_classification.xlsx'), col_names = T, format_headers = T)

  dtMammaryRefineEnrichMotifVierstra <- data.table("motifID"=names(mammary_refineEnrichMotifVierstra),
                                                   "mlog10Padj"=do.call(cbind, assays(mammary_refineEnrichMotifVierstra)$mlog10Padj),
                                                   "enrichment"=do.call(cbind, assays(mammary_refineEnrichMotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtMammaryRefineEnrichMotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_mammary_refine.xlsx'), col_names = T, format_headers = T)

  dtCellClassificationEnrichMotifVierstra <- data.table("motifID"=names(cell_classificationEnrichMotifVierstra),
                                                        "mlog10Padj"=do.call(cbind, assays(cell_classificationEnrichMotifVierstra)$mlog10Padj),
                                                        "enrichment"=do.call(cbind, assays(cell_classificationEnrichMotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtCellClassificationEnrichMotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_cell_classification.xlsx'), col_names = T, format_headers = T)

  dtLHor_LArid1a_MotifVierstra <- data.table("motifID"=names(LHor_LArid1a_MotifVierstra),
                                             "mlog10Padj"=do.call(cbind, assays(LHor_LArid1a_MotifVierstra)$mlog10Padj),
                                             "enrichment"=do.call(cbind, assays(LHor_LArid1a_MotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtLHor_LArid1a_MotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_LHor_vs_LHorArid1a.xlsx'), col_names = T, format_headers = T)

  # KO vs WT in L-Hor cells C1, C3
  dtC1_C3_MotifVierstra <- data.table("motifID"=names(C1_C3_MotifVierstra),
                                                 "mlog10Padj"=do.call(cbind, assays(C1_C3_MotifVierstra)$mlog10Padj),
                                                 "enrichment"=do.call(cbind, assays(C1_C3_MotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtC1_C3_MotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_C1_vs_C3.xlsx'), col_names = T, format_headers = T)

  dtC1_C2_MotifVierstra <- data.table("motifID"=names(C1_C2_MotifVierstra),
                                                   "mlog10Padj"=do.call(cbind, assays(C1_C2_MotifVierstra)$mlog10Padj),
                                                   "enrichment"=do.call(cbind, assays(C1_C2_MotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtC1_C2_MotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_C1_vs_C2.xlsx'), col_names = T, format_headers = T)

  dtC2_C3_MotifVierstra <- data.table("motifID"=names(C2_C3_MotifVierstra),
                                                    "mlog10Padj"=do.call(cbind, assays(C2_C3_MotifVierstra)$mlog10Padj),
                                                    "enrichment"=do.call(cbind, assays(C2_C3_MotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtC2_C3_MotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_C2_vs_C3.xlsx'), col_names = T, format_headers = T)

  dtC4_C5_MotifVierstra <- data.table("motifID"=names(C4_C5_MotifVierstra),
                                       "mlog10Padj"=do.call(cbind, assays(C4_C5_MotifVierstra)$mlog10Padj),
                                       "enrichment"=do.call(cbind, assays(C4_C5_MotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtC4_C5_MotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_C4_vs_C5.xlsx'), col_names = T, format_headers = T)

  dtC1C2C3_MotifVierstra <- data.table("motifID"=names(C1C2C3_MotifVierstra),
                                         "mlog10Padj"=do.call(cbind, assays(C1C2C3_MotifVierstra)$mlog10Padj),
                                         "enrichment"=do.call(cbind, assays(C1C2C3_MotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtC1C2C3_MotifVierstra,
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_motifs_vierstra_enriched_in_C1_vs_C2C3.xlsx'), col_names = T, format_headers = T)

  #### Accessible regions Heatmaps ####
  # samples
  htB <- DoHeatmapEnhance(marker.peaks = sampleMarkersPeaks,
                          plotLog2FC = T,
                          orientation.vertical = F,
                          marker.cutoff = marker.sample.cutOff, #"FDR <= 0.05 & Log2FC >= 0.5",
                          sample.feature.ratio = .8,
                          show.group.dendrogram = F,
                          feature.fontsize = 8,
                          archRSortDescending = T,
                          archRSort = T, # False so we keep order largest to smallest
                          clustering.within = F,
                          group.colors = bionames_palette_list,
                          heatmap.colors = ArchR::ArchRPalettes$solarExtra,
                          spreadsheets.path = file.path(spreadsheets.dir,'BioNames-Peak-Marker-Heatmap.xlsx'),
                          title = paste0("Peaks per sample"))
  pdf(file = file.path(getOutputDirectory(proj),'Plots','BioNames-One_vs_all-Peak-Marker-Heatmap.pdf'), width = 5.5, height = 2.5, bg = 'transparent')
  draw(htB, heatmap_legend_side = "bot")
  dev.off()

  ## clusters
  htC <- DoHeatmapEnhance(marker.peaks = clusterMarkersPeaks,
                          plotLog2FC = T,
                          orientation.vertical = F,
                          marker.cutoff = marker.sample.cutOff, #"FDR <= 0.05 & Log2FC >= 0.5",
                          sample.feature.ratio = 1,
                          show.group.dendrogram = F,
                          feature.fontsize = 8,
                          clustering.within = F,
                          archRSortDescending = T,
                          archRSort = T,
                          group.colors = cluster_palette,
                          heatmap.colors = ArchR::ArchRPalettes$solarExtra,
                          spreadsheets.path = file.path(spreadsheets.dir,'Cluster-Peak-Marker-Heatmap.xlsx'),
                          title = paste0("Peaks per cluster"))
  pdf(file = file.path(getOutputDirectory(proj),'Plots','Cluster-One_vs_all-Peak-Marker-Heatmap.pdf'), width = 5.5, height = 2.5, bg = 'transparent')
  draw(htC, heatmap_legend_side = "bot")
  dev.off()

  ## conditions
  htX <- DoHeatmapEnhance(marker.peaks = conditionMarkersPeaks,
                          plotLog2FC = T,
                          orientation.vertical = F,
                          marker.cutoff = marker.sample.cutOff, #"FDR <= 0.05 & Log2FC >= 0.5",
                          sample.feature.ratio = 1,
                          show.group.dendrogram = F,
                          feature.fontsize = 8,
                          clustering.within = F,
                          archRSort = T,
                          archRSortDescending = T,
                          group.colors = condition_palette,
                          heatmap.colors = ArchR::ArchRPalettes$solarExtra,
                          spreadsheets.path = file.path(spreadsheets.dir,'Condition-Peak-Marker-Heatmap.xlsx'),
                          title = paste0("Peaks per condition"))
  pdf(file = file.path(getOutputDirectory(proj),'Plots','Condition-One_vs_all-Peak-Marker-Heatmap.pdf'), width = 5.5, height = 1.8, bg = 'transparent')
  draw(htX, heatmap_legend_side = "bot")
  dev.off()

  ## mammary celltype
  htT <- DoHeatmapEnhance(marker.peaks = mammary_refineMarkersPeaks,
                          plotLog2FC = T,
                          orientation.vertical = F,
                          marker.cutoff = marker.sample.cutOff, #"FDR <= 0.05 & Log2FC >= 0.5",
                          sample.feature.ratio = 1,
                          show.group.dendrogram = F,
                          feature.fontsize = 8,
                          clustering.within = F,
                          archRSort = T,
                          archRSortDescending = F,
                          group.colors = mammary_refine_clusters_palette,
                          heatmap.colors = ArchR::ArchRPalettes$solarExtra,
                          spreadsheets.path = file.path(spreadsheets.dir,'Mammary-Peak-Marker-Heatmap.xlsx'),
                          title = paste0("Peaks per mammary cell state"))
  pdf(file = file.path(getOutputDirectory(proj),'Plots','Mammary-One_vs_all-Peak-Marker-Heatmap.pdf'), width = 5.5, height = 2.2, bg = 'transparent')
  draw(htT, heatmap_legend_side = "bot")
  dev.off()

  ## celltype
  htT <- DoHeatmapEnhance(marker.peaks = cell_classificationMarkersPeaks,
                          orientation.vertical = F,
                          marker.cutoff = marker.sample.cutOff, #"FDR <= 0.05 & Log2FC >= 0.5",
                          sample.feature.ratio = 1,
                          show.group.dendrogram = F,
                          group.fontsize = 6,
                          clustering.within = F,
                          archRSort = T,
                          archRSortDescending = F,
                          group.colors = luminal_basal_palette,
                          heatmap.colors = ArchR::ArchRPalettes$solarExtra,
                          spreadsheets.path = file.path(spreadsheets.dir,'CellState-Peak-Marker-Heatmap.xlsx'),
                          title = paste0("Peaks per cell state"))
  pdf(file = file.path(getOutputDirectory(proj),'Plots','CellState-One_vs_all-Peak-Marker-Heatmap.pdf'), width = 7, height = 2.5, bg = 'transparent')
  draw(htT, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 20), "mm")) #tlbr
  dev.off()

  ## C3C1 - not more informative than just looking at cluster vs others
  htX <- DoHeatmapEnhance(marker.peaks = cbind(C1C3, C3C1),
                          orientation.vertical = F,
                          marker.cutoff = marker.sample.cutOff, #"FDR <= 0.05 & Log2FC >= 0.5",
                          sample.feature.ratio = 1,
                          show.group.dendrogram = F,
                          feature.fontsize = 12,
                          clustering.within = F,
                          archRSort = T,
                          archRSortDescending = T,
                          group.colors = cluster_palette,
                          heatmap.colors = ArchR::ArchRPalettes$solarExtra,
                          spreadsheets.path = file.path(spreadsheets.dir,'C1_C3-Peak-Marker-Heatmap.xlsx'),
                          title = paste0("Peaks per cluster"))
  pdf(file = file.path(getOutputDirectory(proj),'Plots','C1_C3-One_vs_all-Peak-Marker-Heatmap.pdf'), width = 6, height = 2.5, bg = 'transparent')
  draw(htX, heatmap_legend_side = "bot")
  dev.off()

  #### Write bed file of cluster peaks ####
  marker.cluster.cutOff
  clusterMarkerList <- getMarkers(clusterMarkersPeaks, cutOff = marker.cluster.cutOff) # used in plotbrowser track

  lengths(clusterMarkerList)
  library(rtracklayer)
  (outDir <- file.path(getOutputDirectory(proj),'Annotations'))
  for(i in names(clusterMarkerList)){
    if(nrow(clusterMarkerList[[i]])==0)
      next
    gr <- GRanges(clusterMarkerList[[i]])
    short.name <- gsub("-.*","",i)
    gr@ranges@NAMES <- short.name
    export.bed(gr, con = file.path(outDir, paste0(short.name,'.bed')))
  }

  #### Write bed file of condition peaks ####
  marker.cluster.cutOff
  conditionMarkerList <- getMarkers(conditionMarkersPeaks, cutOff = marker.cluster.cutOff) # used in plotbrowser track

  lengths(conditionMarkerList)
  library(rtracklayer)
  (outDir <- file.path(getOutputDirectory(proj),'Annotations'))
  for(i in names(conditionMarkerList)){
    if(nrow(conditionMarkerList[[i]])==0)
      next
    gr <- GRanges(conditionMarkerList[[i]])
    short.name <- gsub("-.*","",i)
    gr@ranges@NAMES <- short.name
    export.bed(gr, con = file.path(outDir, paste0(short.name,'.bed')))
  }

  #### Write bed file of cell classification peaks ####
  marker.cluster.cutOff
  cell_classificationMarkerList <- getMarkers(cell_classificationMarkersPeaks, cutOff = marker.cluster.cutOff) # used in plotbrowser track

  lengths(cell_classificationMarkerList)
  library(rtracklayer)
  (outDir <- file.path(getOutputDirectory(proj),'Annotations'))
  for(i in names(cell_classificationMarkerList)){
    if(nrow(cell_classificationMarkerList[[i]])==0)
      next
    gr <- GRanges(cell_classificationMarkerList[[i]])
    short.name <- gsub("-.*","",i)
    gr@ranges@NAMES <- short.name
    export.bed(gr, con = file.path(outDir, paste0(short.name,'.bed')))
  }

  #### Ma plots of peak markers ####
  cat("Using cutoff=",abs.marker.cluster.cutOff)
  pma <- mclapply( unique(proj$Clusters)[order(as.integer(gsub("^C|-.*","",unique(proj$Clusters))))], function(x){
    plotMarkers(seMarker = clusterMarkersPeaks, name = x, cutOff = abs.marker.cluster.cutOff, plotAs = "MA") +
      labs(subtitle = paste("Cluster",x)) + theme(plot.margin=unit(c(0,0.1,0,0), "cm"), legend.position = 'none')
  })
  for( i in seq_along(pma) ){
    pma[[i]]$labels$title <- gsub(" Percent","%", pma[[i]]$labels$title)
  }

  png(file = file.path(getOutputDirectory(proj),'Plots','MA_Peaks_in_clusters_one_vs_all.png'),
      width = 12, height = 6, units = 'in', bg='transparent', res = 200)
  do.call(grid.arrange, c(pma, nrow=2))
  dev.off()

  #------------------------------------------- Sequence logo collections -------------------------------------------
  # motifPWMs <- readRDS(file = file.path(robj.dir,"Vierstra-Human-Motifs-Update.rds"))

  cismotifPWMs <- proj@peakAnnotation@listData$cisbp$motifs

  motifPWMs <- proj@peakAnnotation@listData$Vierstra$motifs

  g <- 'Ari'
  pgL <- PlotMotifSeqLogo(motifList = cismotifPWMs[grep(g, names(cismotifPWMs), ignore.case = T, value = T)])
  length(pgL)
  cairo_pdf(filename = file.path(plot.dir,paste0('MotifLogo_',g,'Like.pdf')), width = 14, height = 5)
  do.call(grid.arrange, c(pgL, nrow=2))
  dev.off()

  g <- 'Ari'
  pgL <- PlotMotifSeqLogo(motifList = motifPWMs[grep(g, names(motifPWMs), ignore.case = T, value = T)])
  length(pgL)
  cairo_pdf(filename = file.path(plot.dir,paste0('MotifLogo_',g,'Like.pdf')), width = 14, height = 5)
  do.call(grid.arrange, c(pgL, nrow=2))
  dev.off()

  pBHlh <- PlotMotifSeqLogo(motifList = motifPWMs[grep("Bhlh", ignore.case = T, names(motifPWMs), value = T)])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_BhlhLike1.pdf'), width = 35, height = 12, onefile = T)
  do.call(grid.arrange, c(pBHlh[1:50], nrow=5))
  do.call(grid.arrange, c(pBHlh[51:100], nrow=5))
  do.call(grid.arrange, c(pBHlh[101:150], nrow=5))
  do.call(grid.arrange, c(pBHlh[151:length(pBHlh)], nrow=5))
  dev.off()

  pETSh <- PlotMotifSeqLogo(motifList = motifPWMs[grep("ETS", ignore.case = T, names(motifPWMs), value = T)])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_ETSLike1.pdf'), width = 35, height = 10, onefile = T)
  do.call(grid.arrange, c(pETSh[1:50], nrow=5))
  do.call(grid.arrange, c(pETSh[51:length(pETSh)], nrow=5))
  dev.off()

  pPouh <- PlotMotifSeqLogo(motifList = motifPWMs[grep("POU", ignore.case = T, names(motifPWMs), value = T)])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_PouLike.pdf'), width = 35, height = 15)
  do.call(grid.arrange, c(pPouh, nrow=6))
  dev.off()

  pArid <- PlotMotifSeqLogo(motifList = motifPWMs[grep("Arid", ignore.case = T, names(motifPWMs), value = T)])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_AridLike.pdf'), width = 8, height = 10)
  do.call(grid.arrange, c(pArid, nrow=5))
  dev.off()

  pZscan4L <- PlotMotifSeqLogo(motifList = motifPWMs[grep("Zscan4", ignore.case = T, names(motifPWMs), value = T)])
  pZscan4 <- PlotMotifSeqLogo(chromVARmotifs::mouse_pwms_v2[grep('Zscan4', ignore.case = T, names(chromVARmotifs::mouse_pwms_v2))])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_Zscan4Like.pdf'), width = 8, height = 10)
  do.call(grid.arrange, c(pZscan4L, pZscan4, nrow=5))
  dev.off()

  pfoxcis <- PlotMotifSeqLogo(chromVARmotifs::mouse_pwms_v2[grep('Foxa', ignore.case = T, names(chromVARmotifs::mouse_pwms_v2))])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_Foxa1CisLike.pdf'), width = 8, height = 10)
  do.call(grid.arrange, c(pfoxcis, nrow=3))
  dev.off()

  pSmarccL <- PlotMotifSeqLogo(motifList = motifPWMs[grep("Smarcc", ignore.case = T, names(motifPWMs), value = T)])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_smarcc.pdf'), width = 10, height = 3)
  do.call(grid.arrange, c(pSmarccL, nrow=1))
  dev.off()

  pFox4 <- PlotMotifSeqLogo(motifList = motifPWMs[grep("FOX/4", ignore.case = T, names(motifPWMs), value = T)])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_Fox4.pdf'), width = 16, height = 10)
  do.call(grid.arrange, c(pFox4, nrow=6))
  dev.off()

  pFox5 <- PlotMotifSeqLogo(motifList = motifPWMs[grep("FOX/5", ignore.case = T, names(motifPWMs), value = T)])
  cairo_pdf(filename = file.path(plot.dir,'MotifLogo_Fox5.pdf'), width = 24, height = 10)
  do.call(grid.arrange, c(pFox5, nrow=6))
  dev.off()

  #------------------------------------------- Motif Significant Enrichment Heatmaps -------------------------------------------
  nrow(sampleEnrichMotifVierstra)
  gv <- plotEnrichHeatmap(sampleEnrichMotifVierstra, n = 100, transpose = TRUE, returnMatrix = F)
  gvm <- plotEnrichHeatmap(sampleEnrichMotifVierstra, n = 50, transpose = TRUE, returnMatrix = T)
  gc <- plotEnrichHeatmap(sampleEnrichMotifCisbp, n = 50, transpose = TRUE, returnMatrix = F)
  gcond <- plotEnrichHeatmap(conditionEnrichMotifVierstra, n = 50, transpose = TRUE, returnMatrix = F)

  ComplexHeatmap::draw(gv, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  ComplexHeatmap::draw(gc, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  ComplexHeatmap::draw(gcond, heatmap_legend_side = "bot", annotation_legend_side = "bot")

  #-----sample enrichment summary ----#
  htH <- plotEnrichMotifHeatmap(sampleEnrichMotifVierstra[,c('Arid1a E2','Arid1a','Arid1a ko E2','Arid1a ko')],
                                n = Inf,
                                mlog10Padj.cutOff = 4,
                                title = 'Enriched T.F. motifs per sample',
                                remove.redundant = T,
                                cluster.columns = F,
                                labelMotifs = F,
                                col = rev(pals::cubehelix()))
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Sample-Heatmap-Horizontal-Vierstra.pdf"),
            width = 7, height = 3, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 20), "mm")) #tlbr
  dev.off()

  ## Motif hits Horizontal - grouped by family
  htH <- plotEnrichMotifHeatmap(sampleEnrichMotifVierstra,
                                n = 30,
                                title = paste("Top",30,"enriched transcription factor motifs per family"),
                                remove.redundant = T,
                                labelMotifClusters= '',
                                cutOff = 20,
                                labelMotifs = T,
                                row_name_fontsize = 8,
                                column_name_fontsize = 9) #rev(pals::cubehelix())
  dims <- dim(htH@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-Sample-Horizontal-Vierstra.pdf"),
            width = ceiling(dims[2]/4), height = 5, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 10), "mm")) #tlbr
  dev.off()

  ## Motif hits vertical - grouped by family
  htv <- plotEnrichMotifHeatmap(sampleEnrichMotifVierstra,
                                n = 50, #Inf,
                                remove.redundant = T,
                                cutOff = 10,
                                orientation_horizontal = F,
                                labelMotifs = T,
                                labelMotifN = T,
                                title = 'Enriched TF motifs per sample',
                                #col = as.vector(rev(pals::cubehelix())),
                                labelMotifClusters= '',
                                row_name_fontsize = 7,
                                outline = 0.1
  )
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Sample-Heatmap-Cluster-Vierstra.pdf"),
            width = 6.5, height = 22, bg = 'transparent', onefile = T, family = 'Helvetica', pointsize = 14)
  draw(htv, padding = unit(c(1, 8, 2, 4), "mm")) #tlbr #heatmap_legend_side = "right",
  dev.off()

  htH <- plotEnrichMotifHeatmap(sampleEnrichMotifVierstra,
                                n = Inf,
                                remove.redundant = T,
                                labelMotifClusters= '',
                                cutOff = 10,
                                orientation_horizontal = T,
                                labelMotifs = F,
                                labelMotifN = T,
                                title = 'Enriched TF motifs per sample',
                                #col = as.vector(rev(pals::cubehelix())),
                                outline = 0.1
  )
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-Sample-Vierstra-nolabel.pdf"),
            width = 12, height = 6, bg = 'transparent', onefile = T, family = 'Helvetica', pointsize = 14)
  htH
  dev.off()

  #-----cluster enrichment summary ----#
  clusterHeatmapEM <- plotEnrichHeatmap(clusterEnrichMotifVierstra, n = 30, transpose = TRUE, returnMatrix = F)
  #ComplexHeatmap::draw(clusterHeatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  htH <- plotEnrichMotifHeatmap(clusterEnrichMotifVierstra[,c("C2","C3","C1","C4","C5")],
                                n = Inf,
                                column_title = 'Enriched T.F. motifs per cluster',
                                remove.redundant = T,
                                cluster.columns = F,
                                cutOff = 10,
                                labelMotifs = F,
                                col = rev(pals::cubehelix()))
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-cluster-Heatmap-Horizontal-Vierstra.pdf"),
            width = 8, height = 4, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 20), "mm")) #tlbr
  dev.off()

  ## Motif hits Horizontal - grouped by family
  htH <- plotEnrichMotifHeatmap(clusterEnrichMotifVierstra,
                                n = 20,
                                column_title = paste("Top",20,"enriched transcription factor motifs per family"),
                                remove.redundant = T,
                                labelMotifClusters= '',
                                cutOff = 50,
                                labelMotifs = T,
                                row_name_fontsize = 14,
                                column_name_fontsize = 9)
  #col = rev(pals::cubehelix()))
  dims <- dim(htH@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-cluster-Horizontal-Vierstra.pdf"),
            width = dims[2]/6, height = 6, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 10), "mm")) #tlbr
  dev.off()

  ## Motif hits Vertical - grouped by family
  htv <- plotEnrichMotifHeatmap(clusterEnrichMotifVierstra,
                                n = Inf,
                                remove.redundant = T,
                                labelMotifClusters= '',
                                cutOff = 30,
                                orientation_horizontal = F,
                                labelMotifs = T,
                                labelMotifN = T,
                                title = 'Enriched TF motifs per sample',
                                #col = as.vector(rev(pals::cubehelix())),
                                outline = 0.1
  )
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-cluster-Vierstra-nolabel.pdf"),
            width = 6, height = 16, bg = 'transparent', onefile = T, family = 'Helvetica', pointsize = 14)
  htv
  dev.off()

  #-----condition enrichment summary ----#
  ## Motif hits Horizontal - grouped by family
  htH <- plotEnrichMotifHeatmap(conditionEnrichMotifVierstra,
                                n = 50,
                                column_title = paste("Top",50,"enriched transcription factor motifs per family"),
                                remove.redundant = T,
                                labelMotifClusters= '',
                                cutOff = 20,
                                labelMotifs = T,
                                row_name_fontsize = 8,
                                scaling = 'sample',
                                column_name_fontsize = 9) #rev(pals::cubehelix())
  dims <- dim(htH@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-Condition-Horizontal-Vierstra.pdf"),
            width = ceiling(dims[2]/4), height = 4, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 10), "mm")) #tlbr
  dev.off()

  ## Motif hits Vertical - grouped by family
  htv <- plotEnrichMotifHeatmap(conditionEnrichMotifVierstra,
                                n = 50, #Inf,
                                remove.redundant = T,
                                cutOff = 10,
                                orientation_horizontal = F,
                                labelMotifs = T,
                                labelMotifN = T,
                                column_title = 'Enriched TF motifs per sample',
                                #col = as.vector(rev(pals::cubehelix())),
                                labelMotifClusters= '',
                                row_name_fontsize = 7,
                                outline = 0.1
  )
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Condition-Heatmap-Cluster-Vierstra.pdf"),
            width = 6.5, height = 22, bg = 'transparent', onefile = T, family = 'Helvetica', pointsize = 14)
  draw(htv, padding = unit(c(1, 8, 2, 4), "mm")) #tlbr #heatmap_legend_side = "right",
  dev.off()

  htv <- plotEnrichMotifHeatmap(sampleEnrichMotifVierstra,
                                n = Inf,
                                remove.redundant = T,
                                labelMotifClusters= '',
                                cutOff = 100,
                                orientation_horizontal = F,
                                labelMotifs = F,
                                labelMotifN = T,
                                column_title = 'Enriched TF motifs per sample',
                                #col = as.vector(rev(pals::cubehelix())),
                                outline = 0.1
  )
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-Sample-Vierstra-nolabel.pdf"),
            width = 6, height = 16, bg = 'transparent', onefile = T, family = 'Helvetica', pointsize = 14)
  htv
  dev.off()

  #-----cell classification enrichment summary ----#
  htH <- plotEnrichMotifHeatmap(cell_classificationEnrichMotifVierstra[,c('Arid1a E2','Arid1a','Arid1a ko E2','Arid1a ko')],
                                n = Inf,
                                column_title = 'Enriched T.F. motifs per sample',
                                remove.redundant = T, cluster.columns = F,
                                cutOff = 10,
                                labelMotifs = F,
                                col = rev(pals::cubehelix()))
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-CellClassification-Heatmap-Horizontal-Vierstra.pdf"),
            width = 7, height = 3, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 20), "mm")) #tlbr
  dev.off()

  htH <- plotEnrichMotifHeatmap(cell_classificationEnrichMotifVierstra,
                                n = 50,
                                column_title = paste("Top",50,"enriched transcription factor motifs per family"),
                                remove.redundant = T,
                                labelMotifClusters= '',
                                cutOff = 20,
                                labelMotifs = T,
                                row_name_fontsize = 8,
                                scaling = 'motif', #sample
                                column_name_fontsize = 9) #rev(pals::cubehelix())
  dims <- dim(htH@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-CellClassification-Horizontal-Vierstra.pdf"),
            width = ceiling(dims[2]/4), height = 4, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 10), "mm")) #tlbr
  dev.off()

  htH <- plotEnrichMotifHeatmap(cell_classificationEnrichMotifVierstra,
                                n = Inf,
                                remove.redundant = T,
                                labelMotifClusters= '',
                                cutOff = 20,
                                orientation_horizontal = T,
                                labelMotifs = F,
                                labelMotifN = T,
                                legend_ncolumns=2,
                                column_title = 'Enriched TF motifs per cell type'
  )
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-CellClassification-Vierstra-nolabel.pdf"),
            width = 16, height = 3, bg = 'transparent', onefile = T, family = 'Helvetica', pointsize = 14)
  htH
  dev.off()

  #------#
  saveArchRProject(ArchRProj = proj, load = T)

  #------#
  #### Bed files near motif ####
  sampleMarkersPeaksList <- getMarkers(
      getMarkerFeatures(
        ArchRProj = proj,
        useMatrix = "PeakMatrix",
        groupBy = "bioNames",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
      ),
      returnGR = F,
      cutOff = marker.sample.cutOff # significantly different peaks
  )
  motifmatch <- readRDS(proj@peakAnnotation@listData$Vierstra$Matches)
  a <- assays(motifmatch)$matches #

  #------------------------------------------- Plot chromVAR -------------------------------------------
  # chromVar is designed for predicting enrichment of TF activity on a per-cell basis from sparse chromatin accessibility data. The two primary outputs of chromVAR are:
  # 1. “deviations” - A deviation is a bias-corrected measurement of how far the per-cell accessibility of a given feature (i.e motif) deviates from the expected accessibility based on the average of all cells or samples.
  # 2. “z-score” - The z-score, also known as a “deviation score” is the z-score for each bias-corrected deviation across all cells.
  #.    The absolute value of the deviation score is correlated with the per-cell read depth. This is because, with more reads, you have higher confidence that the difference in per-cell accessibility of the given feature (i.e. motif) from the expectation is greater than would occur by chance.
  ## stores both deviations and z-scores (deviation score)

  getAvailableMatrices(proj)

  #### Motif <=> Gene ####
  # Filter for motifs that are not well defined or expression is missing
  dtM2Glookup <- Motif2GeneLookupDT(proj = proj, motifpwm = proj@peakAnnotation@listData$Vierstra$motifs)
  # Add as attribute to object
  attr(proj, 'VierstraMotifToGene') <- dtM2Glookup
  # get lookup
  dtM2Glookup <- proj@VierstraMotifToGene

  #------#
  saveArchRProject(ArchRProj = proj, load = T)
  #------#

  #### Top deviations ####
  # CisBP
  plotCisBPVarDev <- getVarDeviations(proj, name = "cisbpMatrix", plot = TRUE, n = 100)
  dtCisBPVarDevFilt <- data.table(plotCisBPVarDev$data)
  setnames(dtCisBPVarDevFilt, 'name','motifID')
  (PlotDeviations(dt = dtCisBPVarDevFilt, n=50, text.adjust = T)) # Plots Top deviations

  # Conclusion too much redundancy let's try to filter for expressed genes
  #### Remove some redundancy in motifs ####
  plotVierstraVarDev <- getVarDeviations(proj, name = "VierstraMatrix", plot = TRUE, n = 100)
  dtVierstraVarDevFilt <- data.table(plotVierstraVarDev$data) %>% setnames(., 'name','motifID')
  dtVierstraVarDevFilt <- merge(dtVierstraVarDevFilt, proj@VierstraMotifToGene, by='motifID')
  dtVierstraVarDevFilt[expressed==T]

  dtVierstraVarDevFilt <- dtVierstraVarDevFilt[expressed==T]
  dtVierstraVarDevFilt$family <- gsub(".*#","", dtVierstraVarDevFilt$motifID) %>% gsub("\\|[0-9]+\\|[01]$","",.)
  dtVierstraVarDevFilt[,family := gsub("-[0-9]+$","",family)]
  dtVierstraVarDevFilt <- dtVierstraVarDevFilt[order(rank, decreasing = F)]
  dtVierstraVarDevFilt[,idx := seq_along(.I)] # update idx
  (PlotDeviations(dt = dtVierstraVarDevFilt, n=50, text.adjust = T)) # Plots Top deviations

  #------------------------------------------- Correlate TF and Gene -------------------------------------------
  # vassay contains chromVar deviations and z-scores
  vassay <- getMatrixFromProject(ArchRProj = proj, useMatrix = 'VierstraMatrix', binarize = F, threads = 6)
  # expassay contains gene expression from Seurat
  expassay <- getMatrixFromProject(ArchRProj = proj, useMatrix = 'GeneExpressionMatrix', binarize = F, threads = 6)

  # We do not assign any significance to motif enrichment because we are not comparing groups
  dtMotifGeneCor <- correlateAccessibilityExpression(build=T,
                                                     ArchRproj = proj,
                                                     dtInputChromVar = dtVierstraVarDevFilt,
                                                     vassay = vassay,
                                                     expassay = expassay,
                                                     min.expression = 0,
                                                     cor.method = 'spearman',
                                                     xlsx='dtMotifGeneCorNoImpute.xlsx',
                                                     ncores = 24)

  png(filename = file.path(plot.dir,'Access_expression_corr.png'), width = 10, height = 8, bg = 'transparent', res = 150, units = 'in')
  ggplot(dtMotifGeneCor, aes(x=cor_mlog10Padj, y=coefficient)) +
    geom_label_repel(data = . %>% dplyr::mutate(label = ifelse(abs(coefficient) < 0.2, "", gsub("#.*","",motifID))),
                     aes(label = label),
                     box.padding   = 0.01,
                     point.padding = 0.01,
                     label.padding = 0.01,
                     label.size = 0,
                     color='blue',
                     segment.color = 'grey50',
                     show.legend = FALSE, #this removes the 'a' from the legend)
                     size=2) +
    geom_vline(xintercept = -log(0.01+.Machine$double.xmin)/log(10), col='red', linetype='dotted') +
    geom_point(aes(color=avg_accessibility), shape=20) +
    scale_color_viridis_c("Avg. accessibility", option = 'C') +
    scale_x_continuous(breaks = seq(from = 0, to = max(round(dtMotifGeneCor$cor_mlog10Padj)), by=50)) +
    xlab("mlog10 pearson adjusted pvalue") +
    ylab("pearson's r") +
    theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=1))
  dev.off()

  #------------------------------------------- Variance Rank Filter Grid by condition-------------------------------------------
  # Instead of chromvar let's use the summary statistic
  dtVierstraVarDevFilt <- setDT(read_xlsx(file.path(spreadsheets.dir, 'archr_multiome_atac_motifs_vierstra.xlsx')))
  # Enrichment per sample or cluster (created in 2. Motif Significant Enrichment)
  dtConditionEnrichMotifVierstra <- setDT(read_xlsx(file.path(spreadsheets.dir, 'archr_multiome_atac_motifs_vierstra_enriched_in_condition.xlsx')))
  dtMotifGeneCorConditionNoImputeSpearman <- correlateAccessibilityExpression(build=T,
                                                                              ArchRproj = proj,
                                                                              dtInputChromVar = dtVierstraVarDevFilt,
                                                                              vassay = vassay,
                                                                              expassay = expassay,
                                                                              #min.expression = 0.00,
                                                                              #min.dev.zscore = 1,
                                                                              groupby = 'Condition',
                                                                              ncores = 16,
                                                                              cor.method = 'spearman',
                                                                              xlsx='dtMotifGeneCorConditionNoImputeSpearman.xlsx')
  dtCor <- dtMotifGeneCorConditionNoImputeSpearman
  # This simply combines correlation and chromvar statistics with sample enrichment
  dt <- dtConditionEnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtConditionEnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Condition', value.name = 'mlog10Padj')

  ms <- merge(dtCor, dt, by=c('motifID','Condition'))
  ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]
  ms[,score:= mlog10Padj * sign(avg_accessibility) * frac_exp]
  ms <- ms[order(abs(score), m10log_fisher_p, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("/.*","",family)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  pval.cutoff <- 1e-3
  (min.cutoff <- round(-log10(pval.cutoff)))
  # Gives top short one short motifID and short family ID
  # Then filters by criteria and returns top N short family and luminal/basal
  ms[,plot := F]
  mark <- ms[cor_padj < 0.01 &
               score >= min.cutoff &
               frac_both > 0.1 &
               sign(avg_accessibility) == sign(avg_expression) &
               coefficient > 0]
  #mark <- mark[mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID, bioNames)][,ridx]
  mark <- mark[mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID)][,ridx]
  mark <- ms[ridx %in% mark]
  setkey(ms, motifID, family, motifGene, geneCorrelated)
  mark <- lapply(mark$ridx, function(i){
    r <- mark[ridx==i, list(motifID, family, motifGene, geneCorrelated)]
    ms[r, unlist(ridx)]
  }) %>% unlist

  ms[ridx %in% mark, plot := T]
  nrow(ms[plot==T])
  unique(ms[plot==T, motifID])
  ms[plot==T]
  ms[motifID == 'ANDR_18#FOX/4-forkhead-79-1' & plot==T]
  ms[motifID == 'ANDR_19#NR/20-nuclearreceptor-248-1' & plot==T]
  ms[grep("gata3", motifID, ignore.case = T)][plot==T]
  ms[plot==T & short_motifID %in% c("P53","P63","TP63","TP73") & score > 3]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'Condition_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)

  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'Condition_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))
  # Plots the 'best' motif in grid
  unique(ms[plot==T, motifID])

  msmat <- as.matrix(dcast.data.table(ms[plot==T], motifID ~ Condition, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- NA
  msmat[is.infinite(msmat)] <- NA

  # set up annotation
  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, Condition)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Condition)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Condition, value.var = 'same'), rownames="motifID")
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))
  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                cell.outline = T,
                                labelMotifClusters= '',
                                labelMotifN = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Enriched TF motifs per condition\npadj <",pval.cutoff),
                                title.legend = "Enrichment score",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                na.col = 'grey80',
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  cairo_pdf(filename = file.path(plot.dir, paste0("Motif-Heatmap-condition-vertical_Vierstra_family_top_",top_familyN,".pdf")),
            width = 5.5, height = 9, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="Condition",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = T,
                                    show.group.bar.legend=T,
                                    title = paste("Enriched TF motifs per condition\npadj <",pval.cutoff),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 12,
                                    na.col = 'gray90',
                                    legend.ncol = 4,
                                    col = solarExtraWhite_palette,
                                    group.colors = condition_palette,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())
  fname <- file.path(plot.dir, paste0("Motif-SCHeatmap-condition-vertical_Vierstra_family.png"))
  message("Outputting file:",fname)
  res <- 200
  w <- res*11 # res * inches
  h <- res*10
  fig <- image_graph(width = w, height = h, res = res)
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 10, 2), "mm")) #bltr
  dev.off()
  image_write(image = fig, path = fname, format = 'png')

  pdf(file.path(plot.dir, paste0("Motif-SCHeatmap-condition-vertical_Vierstra_family.pdf")), width=11, height=10)
  htV
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[plot==T]

  #registerDoMC(8)
  mclapply(unique(dt$bioNames), function(i){
    toplot <- dt[bioNames==i & sign(avg_accessibility) > 0]
    dir.create(path = file.path(plot.dir,i), recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal1 = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              imputeValues = F,
                                              nrows = nrow)
      png(filename = file.path(plot.dir,i,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png')),
          width = (lengths(gpGene)*3.4)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 8)

  ### correlation scatter plot ###
  # Correlation pvalue ~ correlation
  g <- lapply(names(bionames_palette_list), function(m){
    toplot <- ms[bioNames==m]
    toplot[,label:=ifelse(abs(coefficient) < 0.2, "", paste(gsub("#.*","",motifID),geneCorrelated, sep='/'))]
    ggplot(toplot, aes(x=coefficient, y=cor_mlog10Padj, label=label)) +
      geom_vline(xintercept = 0, linetype='dashed') +
      geom_text_repel(
        data          = subset(toplot, coefficient > 0.2 & cor_mlog10Padj > 2),
        nudge_x       = .3 + subset(toplot, coefficient > 0.2 & cor_mlog10Padj > 2)$coefficient,
        seed = 3,
        segment.size  = 0.2,
        color='blue',
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        size         = 3,
        show.legend = F,
        max.overlaps = 10,
      ) +
      geom_text_repel(
        data          = subset(toplot, coefficient < -0.2 & cor_mlog10Padj > 2),
        nudge_x       = -.1 - subset(toplot, coefficient < -0.2 & cor_mlog10Padj > 2)$coefficient,
        seed = 3,
        segment.size  = 0.2,
        color='blue',
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        size          = 3,
        max.overlaps = 10,
        show.legend = F
      ) +
      geom_hline(yintercept = -log(0.01)/log(10), col='red', linetype='dotted') +
      geom_point(aes(color=bioNames), shape=20, size=.8, show.legend = F) +
      scale_color_manual("", values = bionames_palette_list) +
      facet_wrap(~bioNames, scales = 'free_y', nrow = 1) +
      scale_y_continuous(breaks = seq(from = 0, to = max(ceiling(toplot$cor_mlog10Padj)), by=20)) +
      ylab("-log10 adjusted pvalues (Pearson correlation)") +
      xlab("pearson's r") +
      theme_bw() + theme(panel.grid = element_blank())
  })
  png(filename = file.path(plot.dir,'Access_bionames_expression_cor.png'), width = 25, height = 16,
      bg = 'transparent', res = 150, units = 'in')
  do.call(grid.arrange, c(g, nrow=3))
  dev.off()

  #------------------------------------------- ChromVar nominate by condition -------------------------------------------
  (clu <- unique(proj$Condition))
  # for each cell type
  ztests <- lapply(clu, function(x){
    xcells <- rownames(proj@cellColData[proj@cellColData$Condition %in% x,])
    othercells <- rownames(proj@cellColData[!proj@cellColData$Condition %in% x,])
    zx <- assays(vassay)$z[,xcells]
    zothers <- assays(vassay)$z[,othercells]
    ztest <- mclapply(rownames(zx), function(m){
      test <- wilcox.test(zx[m,], zothers[m,]) # we want sig. up and down and filter away later
      data.table(cluster=x, motifID=m, zdev_p_value=test$p.value)
    }, mc.cores = 8) %>% rbindlist()
    # every motif is tests so we adjust pvalue pvalue
    ztest[,zdev_adjust := p.adjust(p = zdev_p_value, method = "fdr")]
  }) %>% rbindlist()
  setnames(x = ztests, old = "cluster", new = "Condition")
  # each cluster is tested against all others so we adjust for that test
  ztests[,zdev_adjust := p.adjust(p = zdev_adjust, method = "fdr"), by=Condition]

  write.table(x = ztests, file = file.path(data.dir,"ztests_Condition.txt"), sep = "\t", row.names = F, col.names = T, quote = T)

  dtVierstraVarDevFilt <- setDT(read_xlsx(file.path(spreadsheets.dir, 'archr_multiome_atac_motifs_vierstra.xlsx')))
  dtConditionEnrichMotifVierstra <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'archr_multiome_atac_motifs_vierstra_enriched_in_condition.xlsx')))
  dtConditionNoImputeSpearman <- correlateAccessibilityExpression(build=T,
                                                                  ArchRproj = proj,
                                                                  dtInputChromVar = dtVierstraVarDevFilt,
                                                                  vassay = vassay,
                                                                  expassay = expassay,
                                                                  min.expression = -1,
                                                                  groupby = 'Condition',
                                                                  ncores = 16,
                                                                  cor.method = 'spearman',
                                                                  xlsx='dtMotifGeneCorConditionNoImputeSpearman.xlsx')

  dtCor <- dtConditionNoImputeSpearman
  dt <- dtConditionEnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtConditionEnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Condition', value.name = 'mlog10Padj')

  ms <- merge(dtCor, dt, by=c('motifID','Condition'))
  ms <- merge(ms, ztests, by=c("motifID","Condition"))
  ms[,mlog10_zdev_adjust := -log10(max(c(zdev_adjust, .Machine$double.xmin))),by=1:nrow(ms)]
  ms[,score:= mlog10_zdev_adjust * sign(avg_accessibility) * frac_exp]
  ms <- ms[order(abs(score), mlog10_zdev_adjust, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("-.*","",family)]
  ms[,short_familyID := gsub("/.*","",short_familyID)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  pval.cutoff <- 10^-10
  (min.cutoff <- round(-log10(pval.cutoff)))
  ms[,plot := F]
  msfilt <- ms[cor_padj < 0.05 &
                 score >= min.cutoff &
                 frac_both > 0.1 &
                 sign(avg_accessibility) == sign(avg_expression) &
                 coefficient > 0 &
                 avg_expression > 0.1]
  msfilt <- msfilt[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  mark <- msfilt[msfilt[, .I[head(top_familyN)], by=list(short_familyID, Condition)]$V1]
  mark <- msfilt[ridx %in% mark$ridx]
  mark <- mark[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  setkey(mark, motifID, family, motifGene, geneCorrelated)

  ms[ridx %in% mark$ridx, plot := T]
  nrow(ms[plot==T])
  ms[family == 'FOX/4-forkhead-79' & plot==T,]
  ms[motifID=="MEIS1_1031#HD/14-homeodomain-19-0" & plot==T]
  ms[short_familyID=="HD/14" & plot==T]
  ms[short_familyID=="MIES" & plot==T]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'Condition_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'Condition_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Top n motif per group
  topn <- 30
  msplot <- ms[plot==T][order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  # These are top
  msplot <- msplot[msplot [, .I[seq_len(topn)], by = .(Condition)]$V1]

  # Plots the 'best' motif in grid
  msmat <- as.matrix(dcast.data.table(ms[motifID %in% msplot$motifID], motifID ~ Condition, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(score)], by=list(motifID, Condition)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Condition)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Condition, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(rownames(msmat)),]
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = F,
                                cluster.rows = T,
                                add.score.to.name = F,
                                labelMotifClusters= '',
                                labelMotifN = F,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Top",topn,"enriched TF motifs per cell type"),
                                title.legend = paste0("padj <",pval.cutoff),
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap_Condition-chromvar-spearman.pdf"),
            width = 5, height = 14, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="Condition",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    show.group.bar.legend=T,
                                    title = paste("Top",topn,"enriched TF motifs per cell type"),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = condition_palette,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-condition-chromvar.png")
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 10*res, height = 19*res, res = res, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  f <- file.path(plot.dir,"Motif-SCHeatmap-condition-chromvar.pdf")
  message("Plotting to:",f)
  pdf(file = f, width = 10, height = 19, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$Condition), function(i){
    toplot <- dt[Condition==i]
    plotdir <- file.path(plot.dir,gsub("/","_",i))
    dir.create(path = plotdir, recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plotdir,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access_chromvar.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- Variance Rank Filter Grid by Sample-------------------------------------------
  # Instead of chromvar let's use the summary statistic
  dtVierstraVarDevFilt <- setDT(read_xlsx(file.path(spreadsheets.dir, 'archr_multiome_atac_motifs_vierstra.xlsx')))
  # Enrichment per sample or cluster (created in 2. Motif Significant Enrichment)
  dtSampleEnrichMotifVierstra <- setDT(read_xlsx(file.path(spreadsheets.dir, 'archr_multiome_atac_motifs_vierstra_enriched_in_samplepeaks.xlsx')))
  dtMotifGeneCorSampleNoImputeSpearman <- correlateAccessibilityExpression(build=T,
                                                                           ArchRproj = proj,
                                                                           dtInputChromVar = dtVierstraVarDevFilt,
                                                                           vassay = vassay,
                                                                           expassay = expassay,
                                                                           #min.expression = 0.00,
                                                                           #min.dev.zscore = 1,
                                                                           groupby = 'bioNames',
                                                                           ncores = 16,
                                                                           cor.method = 'spearman',
                                                                           xlsx='dtMotifGeneCorSampleNoImputeSpearman.xlsx')
  dtCor <- dtMotifGeneCorSampleNoImputeSpearman
  # This simply combines correlation and chromvar statistics with sample enrichment
  dt <- dtSampleEnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtSampleEnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'bioNames', value.name = 'mlog10Padj')

  ms <- merge(dtCor, dt, by=c('motifID','bioNames'))
  # combine p-values
  #ms[,padj := exp(-mlog10Padj*log(10))]
  ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]

  ms[,score:= mlog10Padj * sign(avg_accessibility) * frac_exp]
  ms <- ms[order(abs(score), m10log_fisher_p, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("/.*","",family)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  pval.cutoff <- 1e-3
  (min.cutoff <- round(-log10(pval.cutoff)))
  # Gives top short one short motifID and short family ID
  # Then filters by criteria and returns top N short family and luminal/basal
  ms[,plot := F]
  mark <- ms[cor_padj < 0.01 &
               score >= min.cutoff &
               frac_both > 0.1 &
               sign(avg_accessibility) == sign(avg_expression) &
               coefficient > 0]
  #mark <- mark[mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID, bioNames)][,ridx]
  mark <- mark[mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID)][,ridx]
  mark <- ms[ridx %in% mark]
  setkey(ms, motifID, family, motifGene, geneCorrelated)
  mark <- lapply(mark$ridx, function(i){
    r <- mark[ridx==i, list(motifID, family, motifGene, geneCorrelated)]
    ms[r, unlist(ridx)]
  }) %>% unlist

  ms[ridx %in% mark, plot := T]
  nrow(ms[plot==T])
  unique(ms[plot==T, motifID])
  ms[plot==T]
  ms[motifID == 'ANDR_18#FOX/4-forkhead-79-1' & plot==T]
  ms[motifID == 'ANDR_19#NR/20-nuclearreceptor-248-1' & plot==T]
  ms[grep("gata3", motifID, ignore.case = T)][plot==T]
  ms[plot==T & short_motifID %in% c("P53","P63","TP63","TP73") & score > 3]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'BioNames_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)

  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'BioNames_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))
  # Plots the 'best' motif in grid
  unique(ms[plot==T, motifID])

  msmat <- as.matrix(dcast.data.table(ms[plot==T], motifID ~ bioNames, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- NA
  msmat[is.infinite(msmat)] <- NA

  # set up annotation
  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, bioNames)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, bioNames)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ bioNames, value.var = 'same'), rownames="motifID")
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))
  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                cell.outline = T,
                                labelMotifClusters= '',
                                labelMotifN = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Enriched TF motifs per sample\npadj <",pval.cutoff),
                                title.legend = "Enrichment score",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                na.col = 'grey80',
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  cairo_pdf(filename = file.path(plot.dir, paste0("Motif-Heatmap-sample-vertical_Vierstra_family_top_",top_familyN,".pdf")),
            width = 5.5, height = 9, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="bioNames",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = T,
                                    show.group.bar.legend=T,
                                    title = paste("Enriched TF motifs per sample\npadj <",pval.cutoff),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 12,
                                    na.col = 'gray90',
                                    legend.ncol = 4,
                                    col = solarExtraWhite_palette,
                                    group.colors = bionames_palette_list,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())

  fname <- file.path(plot.dir, paste0("Motif-SCHeatmap-bioNames-vertical_Vierstra_family.png"))
  message("Outputting file:",fname)
  res <- 200
  w <- res*11 # res * inches
  h <- res*10
  fig <- image_graph(width = w, height = h, res = res)
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 10, 2), "mm")) #bltr
  dev.off()
  image_write(image = fig, path = fname, format = 'png')

  htH <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = T,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                labelMotifClusters= '',
                                labelMotifN = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Enriched TF motifs per sample\npadj <",pval.cutoff),
                                title.legend = "Enrichment score",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                col = blue_gray_purple_diverging_palette,
                                na.col = "gray99",
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  dims <- dim(htH@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir, paste0("Motif-Heatmap-sample-horizontal_Vierstra_family_top_",top_familyN,".pdf")),
            width = ceiling(dims[2]/2.5), height = dims[1]*.8, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "right", padding = unit(c(1, 1, 1, 1), "mm"))
  decorate_heatmap_body(names(htH)[1], {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  htH <- plotSingleCellMotifHeatmap(dev,
                                    group.by="bioNames",
                                    orientation_horizontal = T,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = T,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    show.group.bar.legend=T,
                                    title = paste("Enriched TF motifs per sample\npadj <",pval.cutoff),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 12,
                                    na.col = 'gray90',
                                    col = color_blind_heatmap_palette,
                                    group.colors = bionames_palette_list,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())
  png(filename = file.path(plot.dir, paste0("Motif-SCHeatmap-bioNames-horizontal_Vierstra_family.png")),
      width = 15, height = 12, bg = 'transparent', family = 'Helvetica', units = 'in', res = 200)
  draw(htH, padding = unit(c(5, 2, 10, 8), "mm")) #bltr
  decorate_heatmap_body(names(htH)[1], {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  gp <-   ggplot(dtCor, aes(x=cor_mlog10Padj, y=coefficient)) +
    geom_label_repel(data = . %>% dplyr::mutate(label = ifelse(abs(coefficient) < 0.2 & cor_padj > 0.05 & cor_mlog10Padj > 50, "", gsub("#.*","",motifID))),
                     aes(label = label),
                     box.padding   = 0.01,
                     point.padding = 0.01,
                     label.padding = 0.01,
                     label.size = 0,
                     color='blue',
                     segment.color = 'grey50',
                     show.legend = FALSE, #this removes the 'a' from the legend)
                     size=2) +
    geom_vline(xintercept = -log(0.01+.Machine$double.xmin)/log(10), col='red', linetype='dotted') +
    geom_point(aes(color=avg_accessibility), shape=20) +
    scale_color_viridis_c("Avg. accessibility", option = 'C') +
    scale_x_continuous(breaks = seq(from = 0, to = max(round(dtMotifGeneCor$cor_mlog10Padj)), by=50)) +
    xlab("mlog10 pearson adjusted pvalue") +
    ylab("pearson's r") +
    facet_wrap( ~ bioNames) +
    theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=1))

  png(filename = file.path(plot.dir,'Access_expression_corr_by_sample.png'), width = 15, height = 10, bg = 'transparent', res = 150, units = 'in')
  gp
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[plot==T]

  #registerDoMC(8)
  mclapply(unique(dt$bioNames), function(i){
    toplot <- dt[bioNames==i & sign(avg_accessibility) > 0]
    dir.create(path = file.path(plot.dir,i), recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal1 = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              imputeValues = F,
                                              nrows = nrow)
      png(filename = file.path(plot.dir,i,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png')),
          width = (lengths(gpGene)*3.4)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 8)

  ### correlation scatter plot ###
  # Correlation pvalue ~ correlation
  g <- lapply(names(bionames_palette_list), function(m){
    toplot <- ms[bioNames==m]
    toplot[,label:=ifelse(abs(coefficient) < 0.2, "", paste(gsub("#.*","",motifID),geneCorrelated, sep='/'))]
    ggplot(toplot, aes(x=coefficient, y=cor_mlog10Padj, label=label)) +
      geom_vline(xintercept = 0, linetype='dashed') +
      geom_text_repel(
        data          = subset(toplot, coefficient > 0.2 & cor_mlog10Padj > 2),
        nudge_x       = .3 + subset(toplot, coefficient > 0.2 & cor_mlog10Padj > 2)$coefficient,
        seed = 3,
        segment.size  = 0.2,
        color='blue',
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        size         = 3,
        show.legend = F,
        max.overlaps = 10,
      ) +
      geom_text_repel(
        data          = subset(toplot, coefficient < -0.2 & cor_mlog10Padj > 2),
        nudge_x       = -.1 - subset(toplot, coefficient < -0.2 & cor_mlog10Padj > 2)$coefficient,
        seed = 3,
        segment.size  = 0.2,
        color='blue',
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        size          = 3,
        max.overlaps = 10,
        show.legend = F
      ) +
      geom_hline(yintercept = -log(0.01)/log(10), col='red', linetype='dotted') +
      geom_point(aes(color=bioNames), shape=20, size=.8, show.legend = F) +
      scale_color_manual("", values = bionames_palette_list) +
      facet_wrap(~bioNames, scales = 'free_y', nrow = 1) +
      scale_y_continuous(breaks = seq(from = 0, to = max(ceiling(toplot$cor_mlog10Padj)), by=20)) +
      ylab("-log10 adjusted pvalues (Pearson correlation)") +
      xlab("pearson's r") +
      theme_bw() + theme(panel.grid = element_blank())
  })
  png(filename = file.path(plot.dir,'Access_bionames_expression_cor.png'), width = 25, height = 16,
      bg = 'transparent', res = 150, units = 'in')
  do.call(grid.arrange, c(g, nrow=3))
  dev.off()

  ### correlation scatter plot ###
  # Fisher pval ~ correlation
  g <- lapply(unique(ms$bioNames), function(m){
    #toplot <- ms[bioNames==m & cor_padj < 0.05]
    toplot <- ms[bioNames==m]
    toplot[,label:=""]
    toplot <- toplot[order(m10log_fisher_p, decreasing = T)]
    toplot[,idx:=1:nrow(toplot)]
    topn <- 20

    tolabel <- toplot[rho > .2 & cor_padj < 0.05, head(.SD, 1), by=list(family)]$idx

    toplot[idx %in% tolabel, label:=paste(gsub("#.*","",motifID),geneCorrelated, sep='/')]

    tolabel <- toplot[rho < -.2 & cor_padj < 0.05, head(.SD, 1), by=list(family)]$idx

    toplot[idx %in% tolabel, label:=paste(gsub("#.*","",motifID),geneCorrelated, sep='/')]

    ggplot(toplot, aes(x=rho, y=m10log_fisher_p, label=label)) +
      geom_vline(xintercept = 0, linetype='dashed') +
      geom_vline(xintercept = .2, linetype='dotted', color='grey70') +
      geom_vline(xintercept = -.2, linetype='dotted', color='grey70') +
      geom_hline(yintercept = -log(0.01)/log(10), col='red', linetype='dotted') +
      geom_point(data = toplot[cor_padj>0.05], color='gray80', shape=20, size=.8, show.legend = F) +
      geom_point(data = toplot[cor_padj<0.05], aes(color=bioNames), shape=20, size=.8, show.legend = F) +
      geom_text_repel(
        data          = subset(toplot, rho > 0.2 & m10log_fisher_p > 2),
        nudge_x       = .5 - subset(toplot, rho > 0.2 & m10log_fisher_p > 2)$rho,
        force        = 0.5,
        seed = 3,
        segment.size  = 0.2,
        color='black',
        segment.color = "grey50",
        segment.curvature = -1e-20,
        direction     = "y",
        hjust         = 0,
        size         = 3,
        show.legend = F,
        max.overlaps = 15,
      ) +
      geom_text_repel(
        data          = subset(toplot, rho < -0.2 & m10log_fisher_p > 2),
        nudge_x       = -.5 - subset(toplot, rho < -0.2 & m10log_fisher_p > 2)$rho,
        force        = 0.5,
        seed = 3,
        segment.size  = 0.2,
        color='black',
        segment.color = "grey50",
        segment.curvature = -1e-20,
        direction     = "y",
        hjust         = 1,
        size          = 3,
        max.overlaps = 15,
        show.legend = F
      ) +
      scale_color_manual("", values = bionames_palette_list) +
      facet_wrap(~bioNames, scales = 'free_y', nrow = 1) +
      scale_x_continuous(limits = c(-1,1)) +
      scale_y_continuous(breaks = seq(from = 0, to = max(ceiling(toplot$m10log_fisher_p)), by=50)) +
      ylab("-log10 adjusted pvalues (fisher combined)") +
      xlab("pearson's r") +
      theme_bw() + theme(panel.grid = element_blank())
  })
  png(filename = file.path(plot.dir,'Access_bionames_expression_corr_fisher.png'), width = 16, height = 8,
      #png(filename = file.path(plot.dir,'Access_bionames_expression_corr_s.png'), width = 25, height = 8,
      bg = 'transparent', res = 150, units = 'in')
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  ### correlation accessibility and scatter plot ###
  g <- lapply(unique(ms$bioNames), function(m){
    #toplot <- ms[bioNames==m & cor_padj < 0.05]
    toplot <- ms[bioNames==m]
    toplot[,label:=""]
    toplot <- toplot[order(mlog10Padj, decreasing = T)]
    toplot[,idx:=1:nrow(toplot)]
    topn <- 20

    tolabel <- toplot[rho > .2 & cor_padj < 0.05, head(.SD, 1), by=list(family)]$idx

    toplot[idx %in% tolabel, label:=paste(gsub("#.*","",motifID),geneCorrelated, sep='/')]

    tolabel <- toplot[rho < -.2 & cor_padj < 0.05, head(.SD, 1), by=list(family)]$idx

    toplot[idx %in% tolabel, label:=paste(gsub("#.*","",motifID),geneCorrelated, sep='/')]

    ggplot(toplot, aes(x=rho, y=mlog10Padj, label=label)) +
      geom_vline(xintercept = 0, linetype='dashed') +
      geom_vline(xintercept = .2, linetype='dotted', color='grey70') +
      geom_vline(xintercept = -.2, linetype='dotted', color='grey70') +
      geom_hline(yintercept = -log(0.01)/log(10), col='red', linetype='dotted') +
      geom_point(data = toplot[cor_padj>0.05], color='gray80', shape=20, size=.8, show.legend = F) +
      geom_point(data = toplot[cor_padj<0.05], aes(color=bioNames), shape=20, size=.8, show.legend = F) +
      geom_text_repel(
        data          = subset(toplot, rho > 0.2 & mlog10Padj > 2),
        nudge_x       = .5 - subset(toplot, rho > 0.2 & mlog10Padj > 2)$rho,
        force        = 0.5,
        seed = 3,
        segment.size  = 0.2,
        color='black',
        segment.color = "grey50",
        segment.curvature = -1e-20,
        direction     = "y",
        hjust         = 0,
        size         = 3,
        show.legend = F,
        max.overlaps = 15,
      ) +
      geom_text_repel(
        data          = subset(toplot, rho < -0.2 & mlog10Padj > 2),
        nudge_x       = -.5 - subset(toplot, rho < -0.2 & mlog10Padj > 2)$rho,
        force        = 0.5,
        seed = 3,
        segment.size  = 0.2,
        color='black',
        segment.color = "grey50",
        segment.curvature = -1e-20,
        direction     = "y",
        hjust         = 1,
        size          = 3,
        max.overlaps = 15,
        show.legend = F
      ) +
      scale_color_manual("", values = bionames_palette_list) +
      facet_wrap(~bioNames, scales = 'free_y', nrow = 1) +
      scale_x_continuous(limits = c(-1,1)) +
      scale_y_continuous(breaks = seq(from = 0, to = max(ceiling(toplot$mlog10Padj)), by=50)) +
      ylab("-log10 adjusted pvalues (accessibility)") +
      xlab("pearson's r") +
      theme_bw() + theme(panel.grid = element_blank())
  })
  png(filename = file.path(plot.dir,'Access_bionames_expression_corr_accessibility.png'), width = 16, height = 8,
      #png(filename = file.path(plot.dir,'Access_bionames_expression_corr_s.png'), width = 25, height = 8,
      bg = 'transparent', res = 150, units = 'in')
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  #------------------------------------------- ChromVar nominate by bioNames -------------------------------------------
  clu <- unique(proj$bioNames)
  # for each cell type
  ztests <- lapply(clu, function(x){
    xcells <- rownames(proj@cellColData[proj@cellColData$bioNames %in% x,])
    othercells <- rownames(proj@cellColData[!proj@cellColData$bioNames %in% x,])
    zx <- assays(vassay)$z[,xcells]
    zothers <- assays(vassay)$z[,othercells]
    ztest <- mclapply(rownames(zx), function(m){
      test <- wilcox.test(zx[m,], zothers[m,]) # we want sig. up and down and filter away later
      data.table(cluster=x, motifID=m, zdev_p_value=test$p.value)
    }, mc.cores = 8) %>% rbindlist()
    # every motif is tests so we adjust pvalue pvalue
    ztest[,zdev_adjust := p.adjust(p = zdev_p_value, method = "fdr")]
  }) %>% rbindlist()
  setnames(x = ztests, old = "cluster", new = "bioNames")
  # each cluster is tested against all others so we adjust for that test
  ztests[,zdev_adjust := p.adjust(p = zdev_adjust, method = "fdr"), by=bioNames]

  write.table(x = ztests, file = file.path(data.dir,"ztests_bioNames.txt"), sep = "\t", row.names = F, col.names = T, quote = T)

  ztests <- fread(file = file.path(data.dir,"ztests_bioNames.txt"), sep = "\t")

  dtbioNamesMotifVierstra <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'archr_multiome_atac_motifs_vierstra_enriched_in_samplepeaks.xlsx')))
  dtMotifGeneCorbioNamesNoImputeSpearman <- correlateAccessibilityExpression(build=F,
                                                                             ArchRproj = proj,
                                                                             dtInputChromVar = dtVierstraVarDevFilt,
                                                                             vassay = vassay,
                                                                             expassay = expassay,
                                                                             min.expression = -1,
                                                                             groupby = 'bioNames',
                                                                             ncores = 16,
                                                                             cor.method = 'spearman',
                                                                             xlsx='dtMotifGeneCorbioNamesNoImputeSpearman.xlsx')
  dtCor <- dtMotifGeneCorbioNamesNoImputeSpearman
  dt <- dtbioNamesMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtbioNamesMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'bioNames', value.name = 'mlog10Padj')

  ms <- merge(dtCor, dt, by=c('motifID','bioNames'))
  ms <- merge(ms, ztests, by=c("motifID","bioNames"))
  ms[,mlog10_zdev_adjust := -log10(max(c(zdev_adjust, .Machine$double.xmin))),by=1:nrow(ms)]
  #ms[,score:= mlog10_zdev_adjust * sign(avg_accessibility) * frac_exp]
  ms[,score:= mlog10_zdev_adjust * avg_accessibility]
  ms <- ms[order(abs(score), mlog10_zdev_adjust, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("-.*","",family)]
  ms[,short_familyID := gsub("/.*","",short_familyID)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  pval.cutoff <- 10^-2
  (min.cutoff <- round(-log10(pval.cutoff)))
  ms[,plot := F]

  # Require gene and TF same
  ms <- ms[motifGene == geneCorrelated]
  ms[avg_expression <=0, score := 0]

  msfilt <- ms[cor_padj < 0.05 &
                 mlog10_zdev_adjust >= min.cutoff &
                 frac_both > 0.1 &
                 sign(avg_accessibility) == sign(avg_expression) &
                 coefficient > 0 &
                 avg_expression > 0.1]
  msfilt <- msfilt[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  # by family and class
  mark <- msfilt[msfilt[, .I[head(top_familyN)], by=list(short_familyID, bioNames)]$V1]
  # just family
  #mark <- msfilt[msfilt[, .I[head(top_familyN)], by=list(short_familyID)]$V1]
  mark <- msfilt[ridx %in% mark$ridx]
  mark <- mark[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  setkey(mark, motifID, family, motifGene, geneCorrelated)
  # Unmark duplicate motifGene
  mark <- mark[!duplicated(motifGene)]

  ms[ridx %in% mark$ridx, plot := T]
  nrow(ms[plot==T])
  ms[family == 'FOX/4-forkhead-79' & plot==T,]
  ms[motifID=="MEIS1_1031#HD/14-homeodomain-19-0" & plot==T]
  ms[short_familyID=="MIES" & plot==T]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'mammary_gland_bioNames_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'mammary_gland_bioNames_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Top n motif per group
  topn <- 60
  msplot <- ms[plot==T][order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  # These are top
  msplot <- msplot[msplot [, .I[seq_len(topn)], by = .(bioNames)]$V1]
  msplot <- msplot[!is.na(motifID)]

  # Plots the 'best' motif in grid
  msmat <- as.matrix(dcast.data.table(ms[motifID %in% msplot$motifID], motifID ~ bioNames, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(score)], by=list(motifID, bioNames)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, bioNames)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ bioNames, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(rownames(msmat)),]
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                add.score.to.name = F,
                                labelMotifClusters= '',
                                labelMotifN = F,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Top",topn,"enriched TF motifs per cell type"),
                                title.legend = paste0("padj <",pval.cutoff),
                                scaling = 'motif',
                                row_name_fontsize = 13,
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap_bioNames-chromvar-spearman.pdf"),
            width = 5, height = 14, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="bioNames",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    show.group.bar.legend=T,
                                    title = paste("Top",topn,"enriched TF motifs per cell type"),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = bionames_palette_list,
                                    column_name_fontsize = 12,
                                    label.col = family_palette) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-bioNames-chromvar.pdf")
  message("Plotting to:",f)
  pdf(file = f, width = 15, height = 11, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "right", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$bioNames), function(i){
    toplot <- dt[bioNames==i]
    plotdir <- file.path(plot.dir,gsub("/","_",i))
    dir.create(path = plotdir, recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plotdir,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access_chromvar.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- ChromVar nominate by Cluster -------------------------------------------
  clu <- unique(proj$Clusters)
  # for each cell type
  ztests <- lapply(clu, function(x){
    # compare zscores in x vs all others
    xcells <- rownames(proj@cellColData[proj@cellColData$Clusters %in% x,])
    othercells <- rownames(proj@cellColData[!proj@cellColData$Clusters %in% x,])
    zx <- assays(vassay)$z[,xcells]
    zothers <- assays(vassay)$z[,othercells]
    ztest <- mclapply(rownames(zx), function(m){
      test <- wilcox.test(zx[m,], zothers[m,])
      data.table(cluster=x, motifID=m, zdev_p_value=test$p.value)
    }, mc.cores = 8) %>% rbindlist()
    ztest[,zdev_adjust := p.adjust(p = zdev_p_value, method = "fdr")]
  }) %>% rbindlist()
  setnames(x = ztests, old = "cluster", new = "Clusters")
  write.table(x = ztests, file = file.path(data.dir,"ztests_clusters.txt"), sep = "\t", row.names = F, col.names = T, quote = T)
  ztests <- fread(file.path(data.dir,"ztests_clusters.txt"))

  c3_cells <- rownames(proj@cellColData[proj@cellColData$Clusters %in% c('C3'),])
  edata <- expassay@assays@data$GeneExpressionMatrix[,c3_cells]
  vdata <- assays(vassay[,c3_cells])$deviations

  c1_cells <- rownames(proj@cellColData[proj@cellColData$Clusters %in% c('C1'),])
  edata <- expassay@assays@data$GeneExpressionMatrix[,c1_cells]
  vdata <- assays(vassay[,c1_cells])$deviations

  lhor_cells <- rownames(proj@cellColData[proj@cellColData$gex_mammary_gland_cluster_classification %in% c('L-Hor'),])
  edata <- expassay@assays@data$GeneExpressionMatrix[,lhor_cells]
  vdata <- assays(vassay[,lhor_cells])$deviations

  rownames(edata) <- rowData(expassay)$name
  identical(names(vdata["MEIS1_1031#HD/14-homeodomain-19-0",]), colnames(edata))
  cor.test(x = vdata["MEIS1_1031#HD/14-homeodomain-19-0",], y=edata["Meis1",], method = "spearman")
  plot(x = vdata["MEIS1_1031#HD/14-homeodomain-19-0",], y=edata["Meis1",])

  cor.test(x = vdata["Gata3_384#GATA-C2H2-242-0",], y=edata["Gata3",])
  plot(x = vdata["Gata3_384#GATA-C2H2-242-0",], y=edata["Gata3",])

  dtVierstraVarDevFilt <- data.table(read_xlsx(file.path(spreadsheets.dir,
                                                         'archr_multiome_atac_motifs_vierstra.xlsx')))
  dtClusterEnrichMotifVierstra <- readxl::read_xlsx(path = file.path(spreadsheets.dir,
                                                                     'archr_multiome_atac_motifs_vierstra_enriched_in_clusterpeaks.xlsx')) %>% setDT()
  dtMotifGeneCorClusterNoImputeSpearman <- correlateAccessibilityExpression(build=T,
                                                                            ArchRproj = proj,
                                                                            dtInputChromVar = dtVierstraVarDevFilt,
                                                                            vassay = vassay,
                                                                            expassay = expassay,
                                                                            min.expression = 0.0,
                                                                            groupby = 'Clusters',
                                                                            ncores = 16,
                                                                            cor.method = 'spearman',
                                                                            xlsx='dtMotifGeneCorClusterNoImputeSpearman.xlsx')
  dtCor <- dtMotifGeneCorClusterNoImputeSpearman
  dt <- dtClusterEnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtClusterEnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Clusters', value.name = 'mlog10Padj')
  ms <- merge(dtCor, dt, by=c('motifID','Clusters'))
  ms <- merge(ms, ztests, by=c("motifID","Clusters"))

  ms[,mlog10_zdev_adjust := -log10(max(c(zdev_adjust, .Machine$double.xmin))),by=1:nrow(ms)]
  ms[,score:= mlog10_zdev_adjust * sign(avg_accessibility) * frac_exp]

  #ms[,score := mlog10Padj * sign(avg_accessibility) * frac_exp]
  ms <- ms[order(abs(score), mlog10_zdev_adjust, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  #ms[,short_familyID := gsub("/.*","",family)]
  ms[,short_familyID := gsub("-.*","",family)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  pval.cutoff <- 10^-30
  (min.cutoff <- round(-log10(pval.cutoff)))
  ms[,plot := F]
  mark <- ms[cor_padj < 0.01 &
               score >= min.cutoff &
               frac_both > 0.1 &
               sign(avg_accessibility) == sign(avg_expression) &
               coefficient > 0]
  mark <- mark[order(score, mlog10_zdev_adjust, mlog10Padj, cor_mlog10Padj, decreasing = T)][mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID, Clusters)][,ridx]
  mark <- ms[ridx %in% mark]
  setkey(ms, motifID, family, motifGene, geneCorrelated)
  mark <- lapply(mark$ridx, function(i){
    r <- mark[ridx==i, list(motifID, family, motifGene, geneCorrelated)]
    ms[r, unlist(ridx)]
  }) %>% unlist

  ms[ridx %in% mark, plot := T]
  nrow(ms[plot==T])
  ms[plot==T]
  ms[family == 'FOX/4-forkhead-79' & plot==T,]
  ms[motifID=="MEIS1_1031#HD/14-homeodomain-19-0" & plot==T]
  ms[motifID == 'ANDR_18#FOX/4-forkhead-79-1' & plot==T,]
  ms[motifID == 'ANDR_19#NR/20-nuclearreceptor-248-1' & plot==T,]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'Clusters_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'Clusters_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Plots the 'best' motif in grid
  msmat <- as.matrix(dcast.data.table(ms[plot==T], motifID ~ Clusters, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, Clusters)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Clusters)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Clusters, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(ms[plot==T, motifID]),]
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = F,
                                cluster.rows = T,
                                add.score.to.name = T,
                                labelMotifClusters= '',
                                labelMotifN = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = "Enriched TF motifs in mammary gland cells",
                                title.legend = paste0("padj <",pval.cutoff),
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap_clusters-chromvar-spearman.pdf"),
            width = 5, height = 14, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="Clusters",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    show.group.bar.legend=T,
                                    title = paste("Enriched TF motifs per type\npadj <",pval.cutoff),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = cluster_palette,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-Clusters-chromvar.png")
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 10*res, height = 18*res, res = res, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  f <- file.path(plot.dir,"Motif-SCHeatmap-Clusters-chromvar.pdf")
  message("Plotting to:",f)
  pdf(file = f, width = 10, height = 18, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$Clusters), function(i){
    toplot <- dt[Clusters==i]
    dir.create(path = file.path(plot.dir,i), recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plot.dir,i,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access_chromvar.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- Variance Rank Filter Grid by Cluster -------------------------------------------
  dtVierstraVarDevFilt <- data.table(read_xlsx(file.path(spreadsheets.dir,
                                                         'archr_multiome_atac_motifs_vierstra.xlsx')))
  dtClusterEnrichMotifVierstra <- readxl::read_xlsx(path = file.path(spreadsheets.dir,
                                                                     'archr_multiome_atac_motifs_vierstra_enriched_in_clusterpeaks.xlsx')) %>% setDT()
  dtMotifGeneCorClusterNoImpute <- correlateAccessibilityExpression(build=T,
                                                                    ArchRproj = proj,
                                                                    dtInputChromVar = dtVierstraVarDevFilt,
                                                                    vassay = vassay,
                                                                    expassay = expassay,
                                                                    #matGE = matGE,
                                                                    min.expression = 0.0,
                                                                    groupby = 'Clusters',
                                                                    ncores = 12,
                                                                    cor.method = 'spearman',
                                                                    xlsx='dtMotifGeneCorClusterNoImputeSpearman.xlsx')
  dtCor <- dtMotifGeneCorClusterNoImpute
  # This simply combines correlation and chromvar statistics with sample enrichment
  dt <- dtClusterEnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtClusterEnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Clusters', value.name = 'mlog10Padj')
  ms <- merge(dtCor, dt, by=c('motifID','Clusters'))
  # combine p-values (not used anymore)
  #ms[,padj := exp(-mlog10Padj*log(10))]
  ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]
  #ms[,score:= (mlog10Padj + cor_mlog10Padj) * sign(avg_accessibility) * frac_exp]

  ms[,score:= mlog10Padj * sign(avg_accessibility) * frac_exp]
  ms <- ms[order(abs(score), m10log_fisher_p, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("/.*","",family)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  pval.cutoff <- 0.05
  min.cutoff <- round(-log10(pval.cutoff))
  ms[,plot := F]
  mark <- ms[cor_padj < 0.01 &
               score >= min.cutoff &
               frac_both > 0.1 &
               sign(avg_accessibility) == sign(avg_expression) &
               coefficient > 0]
  mark <- mark[order(score, decreasing = T)][mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID, Clusters)][,ridx]
  #mark <- mark[mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID)][,ridx]
  mark <- ms[ridx %in% mark]
  setkey(ms, motifID, family, motifGene, geneCorrelated)
  mark <- lapply(mark$ridx, function(i){
    r <- mark[ridx==i, list(motifID, family, motifGene, geneCorrelated)]
    ms[r, unlist(ridx)]
  }) %>% unlist

  ms[ridx %in% mark, plot := T]
  nrow(ms[plot==T])
  ms[plot==T]
  ms[motifID == 'ANDR_18#FOX/4-forkhead-79-1' & plot==T,]
  ms[motifID == 'ANDR_19#NR/20-nuclearreceptor-248-1' & plot==T,]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'Clusters_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'Clusters_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Plots the 'best' motif in grid
  #(msTop <- sort(unique(ms[plot==T, motifID])))

  msmat <- as.matrix(dcast.data.table(ms[plot==T], motifID ~ Clusters, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, Clusters)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Clusters)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Clusters, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(ms[plot==T, motifID]),]
  msgroup[is.na(msgroup)] <- ""
  #msgroup <- msgroup[rownames(msgroup) %in% msTop,]

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                labelMotifClusters= '',
                                labelMotifN = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Enriched TF motifs per sample\npadj <",pval.cutoff),
                                title.legend = "Enrichment score",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                na.col = 'white',
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102) #pals::alphabet(26)
  cairo_pdf(filename = file.path(plot.dir,paste0("Motif-Heatmap-cluster-vertical_Vierstra_top",top_familyN,"_family.pdf")),
            width = 4.5, height = 8, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  range(z)
  quantile(z, seq(0.75,1,by=.01))
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="Clusters",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = T,
                                    quantile.min = .10,
                                    quantile.max = .90,
                                    show.group.bar.legend=T,
                                    title = paste("Enriched TF motifs per cluster\npadj <",pval.cutoff),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = cluster_palette,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-Cluster-vertical_Vierstra_family.png")
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 10*res, height = 9*res, res = res, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  htH <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = T,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                labelMotifClusters= '',
                                labelMotifN = F,
                                labelMotifFrac = T,
                                mlog10Padj.cutOff = 0,
                                labelMotifs = T,
                                add.score.to.name=T,
                                label.outline = .5,
                                title = paste("Enriched TF motifs per sample\npadj <",pval.cutoff),
                                title.legend = "Average accessibility",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                na.col = 'white',
                                col = solarExtraWhite_palette, # blue_gray_purple_diverging_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102) #pals::alphabet(26)
  dims <- dim(htH@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir,paste0("Motif-Heatmap-cluster-horizontal_Vierstra_top",top_familyN,"_family_avg_access.pdf")),
            width = ceiling(dims[2]*.35), height = dims[1]*.4, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "right", padding = unit(c(1, 1, 4, 5), "mm")) # bltr
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$Clusters), function(i){
    toplot <- dt[Clusters==i]
    dir.create(path = file.path(plot.dir,i), recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plot.dir,i,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- Variance Rank Filter by mammary gland class -------------------------------------------
  dtMammaryClassificationEnrichMotifVierstra <- setDT(read_xlsx(path = file.path(spreadsheets.dir,
                                                                        'archr_multiome_atac_motifs_vierstra_enriched_in_mammary_classification.xlsx')))
  dtMotifGeneCorMammaryGlandClusterNoImputeSpearman <- correlateAccessibilityExpression(build=T,
                                                                                 ArchRproj = proj,
                                                                                 dtInputChromVar = dtVierstraVarDevFilt,
                                                                                 vassay = vassay,
                                                                                 expassay = expassay,
                                                                                 min.expression = 0.0,
                                                                                 groupby = 'gex_mammary_gland_cluster_classification',
                                                                                 ncores = 16,
                                                                                 cor.method = 'spearman',
                                                                                 xlsx='dtMotifGeneCorMammaryGlandClusterNoImputeSpearman.xlsx')

  dtCor <- dtMotifGeneCorMammaryGlandClusterNoImputeSpearman
  dt <- dtMammaryClassificationEnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtMammaryClassificationEnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'gex_mammary_gland_cluster_classification', value.name = 'mlog10Padj')
  ms <- merge(dtCor, dt, by=c('motifID','gex_mammary_gland_cluster_classification'))

  ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]

  ms[,score := mlog10Padj * sign(avg_accessibility) * frac_exp]
  ms <- ms[order(abs(score), m10log_fisher_p, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("/.*","",family)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  pval.cutoff <- 0.05
  min.cutoff <- round(-log10(pval.cutoff))
  ms[,plot := F]
  mark <- ms[cor_padj < 0.01 &
               score >= min.cutoff &
               frac_both > 0.1 &
               sign(avg_accessibility) == sign(avg_expression) &
               coefficient > 0]
  mark <- mark[order(score, decreasing = T)][mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID, gex_mammary_gland_cluster_classification)][,ridx]
  mark <- ms[ridx %in% mark]
  setkey(ms, motifID, family, motifGene, geneCorrelated)
  mark <- lapply(mark$ridx, function(i){
    r <- mark[ridx==i, list(motifID, family, motifGene, geneCorrelated)]
    ms[r, unlist(ridx)]
  }) %>% unlist

  ms[ridx %in% mark, plot := T]
  nrow(ms[plot==T])
  ms[plot==T]
  ms[family == 'FOX/4-forkhead' & plot==T,]
  ms[motifID == 'ANDR_18#FOX/4-forkhead-79-1' & plot==T,]
  ms[motifID == 'ANDR_19#NR/20-nuclearreceptor-248-1' & plot==T,]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'mammary_gland_cluster_classification_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'mammary_gland_cluster_classification_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Plots the 'best' motif in grid
  msmat <- as.matrix(dcast.data.table(ms[plot==T], motifID ~ gex_mammary_gland_cluster_classification, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, gex_mammary_gland_cluster_classification)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, gex_mammary_gland_cluster_classification)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ gex_mammary_gland_cluster_classification, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(ms[plot==T, motifID]),]
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = F,
                                cluster.rows = T,
                                add.score.to.name = T,
                                labelMotifClusters= '',
                                labelMotifN = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = "Enriched TF motifs in mammary gland cells",
                                title.legend = paste0("padj <",pval.cutoff),
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap_mammary_gland_classification-vertical-Vierstra-spearman.pdf"),
            width = 4.2, height = 7, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="gex_mammary_gland_cluster_classification",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    show.group.bar.legend=T,
                                    title = paste("Enriched TF motifs per type\npadj <",pval.cutoff),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = mammary_cluster_class_palette,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-mammary-cluster-class-vertical_Vierstra_family.png")
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 10*res, height = 11*res, res = res, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  htH <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = T,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                labelMotifClusters= '',
                                labelMotifN = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Enriched TF motifs in luminal basal cells\npadj <",pval.cutoff),
                                title.legend = "Enrichment score",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  dims <- dim(htH@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap-LuminalBasal-horizontal-Vierstra-spearman.pdf"),
            width = ceiling(dims[2]*.35), height = dims[1]*2, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "right", padding = unit(c(1, 1, 1, 1), "mm")) # bltr
  decorate_heatmap_body(names(htH)[1], {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  # Avg access
  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = F,
                                cluster.rows = T,
                                labelMotifClusters= '',
                                labelMotifN = F,
                                labelMotifFrac = T,
                                mlog10Padj.cutOff = 0,
                                labelMotifs = T,
                                add.score.to.name = F,
                                title = paste("Enriched TF motifs in luminal basal cells\npadj <",pval.cutoff),
                                title.legend = "Average accessibility",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                col = blue_gray_purple_diverging_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  dims <- dim(htV@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap-LuminalBasal-vertical-Vierstra-spearman.pdf"),
            width = dims[2]*2, height = 13, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  decorate_heatmap_body(names(htV)[1], {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  htH <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = T,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                labelMotifClusters= '',
                                labelMotifN = F,
                                labelMotifFrac = T,
                                mlog10Padj.cutOff = 0,
                                labelMotifs = T,
                                add.score.to.name = F,
                                title = paste("Enriched TF motifs in luminal basal cells\npadj <",pval.cutoff),
                                title.legend = "Average accessibility",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  dims <- dim(htH@ht_list[[1]])
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap-LuminalBasal-horizontal-Vierstra-spearman.pdf"),
            width = ceiling(dims[2]*.3), height = dims[1]*1.8, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "right", padding = unit(c(1, 1, 1, 1), "mm")) # bltr
  decorate_heatmap_body(names(htH)[1], {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$gex_mammary_gland_cluster_classification), function(i){
    toplot <- dt[gex_mammary_gland_cluster_classification==i]
    dir.create(path = file.path(plot.dir,i), recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plot.dir,i,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- ChromVar nominate by mammary gland class -------------------------------------------
  clu <- unique(proj$gex_mammary_gland_cluster_classification)
  # for each cell type
  ztests <- lapply(clu, function(x){
    xcells <- rownames(proj@cellColData[proj@cellColData$gex_mammary_gland_cluster_classification %in% x,])
    othercells <- rownames(proj@cellColData[!proj@cellColData$gex_mammary_gland_cluster_classification %in% x,])
    zx <- assays(vassay)$z[,xcells]
    zothers <- assays(vassay)$z[,othercells]
    ztest <- mclapply(rownames(zx), function(m){
      test <- wilcox.test(zx[m,], zothers[m,]) # we want sig. up and down and filter away later
      data.table(cluster=x, motifID=m, zdev_p_value=test$p.value)
    }, mc.cores = 8) %>% rbindlist()
    # every motif is tests so we adjust pvalue pvalue
    ztest[,zdev_adjust := p.adjust(p = zdev_p_value, method = "fdr")]
  }) %>% rbindlist()
  setnames(x = ztests, old = "cluster", new = "gex_mammary_gland_cluster_classification")
  # each cluster is tested against all others so we adjust for that test
  ztests[,zdev_adjust := p.adjust(p = zdev_adjust, method = "fdr"), by=gex_mammary_gland_cluster_classification]

  write.table(x = ztests, file = file.path(data.dir,"ztests_gex_mammary_gland_classification_clusters.txt"), sep = "\t", row.names = F, col.names = T, quote = T)

  dtMammaryClassificationEnrichMotifVierstra <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'archr_multiome_atac_motifs_vierstra_enriched_in_mammary_classification.xlsx')))
  dtMotifGeneCorMammaryGlandClusterNoImputeSpearman <- correlateAccessibilityExpression(build=T,
                                                                                       ArchRproj = proj,
                                                                                       dtInputChromVar = dtVierstraVarDevFilt,
                                                                                       vassay = vassay,
                                                                                       expassay = expassay,
                                                                                       min.expression = -1,
                                                                                       groupby = 'gex_mammary_gland_cluster_classification',
                                                                                       ncores = 16,
                                                                                       cor.method = 'spearman',
                                                                                       xlsx='dtMotifGeneCorMammaryGlandClassificationNoImputeSpearman.xlsx')

  dtCor <- dtMotifGeneCorMammaryGlandClusterNoImputeSpearman
  dt <- dtMammaryRefineEnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtMammaryRefineEnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'gex_mammary_gland_cluster_classification', value.name = 'mlog10Padj')

  ms <- merge(dtCor, dt, by=c('motifID','gex_mammary_gland_cluster_classification'))
  ms <- merge(ms, ztests, by=c("motifID","gex_mammary_gland_cluster_classification"))
  ms[,mlog10_zdev_adjust := -log10(max(c(zdev_adjust, .Machine$double.xmin))),by=1:nrow(ms)]
  ms[,score:= mlog10_zdev_adjust * sign(avg_accessibility) * frac_exp]
  ms <- ms[order(abs(score), mlog10_zdev_adjust, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("-.*","",family)]
  ms[,short_familyID := gsub("/.*","",short_familyID)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  pval.cutoff <- 10^-10
  (min.cutoff <- round(-log10(pval.cutoff)))
  ms[,plot := F]
  msfilt <- ms[cor_padj < 0.05 &
                 score >= min.cutoff &
                 frac_both > 0.1 &
                 sign(avg_accessibility) == sign(avg_expression) &
                 coefficient > 0 &
                 avg_expression > 0.1]
  msfilt <- msfilt[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  mark <- msfilt[msfilt[, .I[head(top_familyN)], by=list(short_familyID, gex_mammary_gland_cluster_classification)]$V1]
  mark <- msfilt[ridx %in% mark$ridx]
  mark <- mark[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  setkey(mark, motifID, family, motifGene, geneCorrelated)

  ms[ridx %in% mark$ridx, plot := T]
  nrow(ms[plot==T])
  ms[family == 'FOX/4-forkhead-79' & plot==T,]
  ms[motifID=="MEIS1_1031#HD/14-homeodomain-19-0" & plot==T]
  ms[short_familyID=="HD/14" & plot==T]
  ms[short_familyID=="MIES" & plot==T]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'mammary_gland_cluster_classification_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'mammary_gland_cluster_classification_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Top n motif per group
  topn <- 30
  msplot <- ms[plot==T][order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  # These are top
  msplot <- msplot[msplot [, .I[seq_len(topn)], by = .(gex_mammary_gland_cluster_classification)]$V1]

  # Plots the 'best' motif in grid
  msmat <- as.matrix(dcast.data.table(ms[motifID %in% msplot$motifID], motifID ~ gex_mammary_gland_cluster_classification, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(score)], by=list(motifID, gex_mammary_gland_cluster_classification)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, gex_mammary_gland_cluster_classification)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ gex_mammary_gland_cluster_classification, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(rownames(msmat)),]
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = F,
                                cluster.rows = T,
                                add.score.to.name = F,
                                labelMotifClusters= '',
                                labelMotifN = F,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Top",topn,"enriched TF motifs per cell type"),
                                title.legend = paste0("padj <",pval.cutoff),
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap_mammary_gland_cluster-chromvar-spearman.pdf"),
            width = 5, height = 14, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="gex_mammary_gland_cluster_classification",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    show.group.bar.legend=T,
                                    title = paste("Top",topn,"enriched TF motifs per cell type"),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = mammary_refine_clusters_palette,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-mammary-cluster-chromvar.png")
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 10*res, height = 19*res, res = res, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  f <- file.path(plot.dir,"Motif-SCHeatmap-mammary-cluster-chromvar.pdf")
  message("Plotting to:",f)
  pdf(file = f, width = 10, height = 19, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$gex_mammary_gland_cluster_classification), function(i){
    toplot <- dt[gex_mammary_gland_cluster_classification==i]
    plotdir <- file.path(plot.dir,gsub("/","_",i))
    dir.create(path = plotdir, recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plotdir,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access_chromvar.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- ChromVar nominate by gex_mammary_gland_refine_clusters -------------------------------------------
  clu <- unique(proj$gex_mammary_gland_refine_clusters)
  # for each cell type
  ztests <- lapply(clu, function(x){
    xcells <- rownames(proj@cellColData[proj@cellColData$gex_mammary_gland_refine_clusters %in% x,])
    othercells <- rownames(proj@cellColData[!proj@cellColData$gex_mammary_gland_refine_clusters %in% x,])
    zx <- assays(vassay)$z[,xcells]
    zothers <- assays(vassay)$z[,othercells]
    ztest <- mclapply(rownames(zx), function(m){
      test <- wilcox.test(zx[m,], zothers[m,]) # we want sig. up and down and filter away later
      data.table(cluster=x, motifID=m, zdev_p_value = test$p.value )
    }, mc.cores = 8) %>% rbindlist()
    # every motif is tests so we adjust pvalue pvalue
    ztest[,zdev_adjust := p.adjust(p = zdev_p_value, method = "fdr")]
  }) %>% rbindlist()
  setnames(x = ztests, old = "cluster", new = "gex_mammary_gland_refine_clusters")
  # each cluster is tested against all others so we adjust for that test
  ztests[,zdev_adjust := p.adjust(p = zdev_adjust, method = "fdr"), by=gex_mammary_gland_refine_clusters]

  ztests[motifID=="FOXA1_480#FOX/4-forkhead-79-0"]

  write.table(x = ztests, file = file.path(data.dir,"ztests_gex_mammary_gland_refine_clusters.txt"), sep = "\t", row.names = F, col.names = T, quote = T)

  ztests <- fread(file = file.path(data.dir,"ztests_gex_mammary_gland_refine_clusters.txt"), sep = "\t")

  dtMammaryRefineEnrichMotifVierstra <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'archr_multiome_atac_motifs_vierstra_enriched_in_mammary_refine.xlsx')))
  dtMotifGeneCorMammaryGlandRefineNoImputeSpearman <- correlateAccessibilityExpression(build=T,
                                                                                        ArchRproj = proj,
                                                                                        dtInputChromVar = dtVierstraVarDevFilt,
                                                                                        vassay = vassay,
                                                                                        expassay = expassay,
                                                                                        #min.expression = -1,
                                                                                        groupby = 'gex_mammary_gland_refine_clusters',
                                                                                        ncores = 16,
                                                                                        cor.method = 'spearman',
                                                                                        xlsx='dtMotifGeneCorMammaryGlandRefineNoImputeSpearman.xlsx')
  dtCor <- dtMotifGeneCorMammaryGlandRefineNoImputeSpearman
  dt <- dtMammaryRefineEnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtMammaryRefineEnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'gex_mammary_gland_refine_clusters', value.name = 'mlog10Padj')

  ms <- merge(dtCor, dt, by=c('motifID','gex_mammary_gland_refine_clusters'))
  ms <- merge(ms, ztests, by=c("motifID","gex_mammary_gland_refine_clusters"))
  ms[,mlog10_zdev_adjust := -log10(max(c(zdev_adjust, .Machine$double.xmin))),by=1:nrow(ms)]

  ms[,score:= mlog10_zdev_adjust * avg_accessibility]

  ms <- ms[order(abs(score), mlog10_zdev_adjust, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("-.*","",family)]
  ms[,short_familyID := gsub("/.*","",short_familyID)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  #pval.cutoff <- 10^-5
  pval.cutoff <- 10^-2
  (min.cutoff <- round(-log10(pval.cutoff)))
  ms[,plot := F]

  # Require gene and TF same
  ms <- ms[motifGene == geneCorrelated]
  ms[avg_expression <=0, score := 0]

  msfilt <- ms[cor_padj < 0.05 &
                 mlog10_zdev_adjust >= min.cutoff &
                 #score >= min.cutoff &
                 frac_both > 0.1 &
                 sign(avg_accessibility) == sign(avg_expression) &
                 coefficient > 0 &
                 avg_expression > 0.1]
  msfilt <- msfilt[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  # by family and class
  mark <- msfilt[msfilt[, .I[head(top_familyN)], by=list(short_familyID, gex_mammary_gland_refine_clusters)]$V1]
  # just family
  #mark <- msfilt[msfilt[, .I[head(top_familyN)], by=list(short_familyID)]$V1]
  mark <- msfilt[ridx %in% mark$ridx]
  mark <- mark[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  setkey(mark, motifID, family, motifGene, geneCorrelated)
  # Unmark duplicate motifGene
  mark <- mark[!duplicated(motifGene)]

  ms[ridx %in% mark$ridx, plot := T]
  nrow(ms[plot==T])
  ms[family == 'FOX/4-forkhead-79' & plot==T,]
  ms[short_familyID=="MIES" & plot==T]
  ms[grep(pattern = "MEIS",x = motifID)][plot==T]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'mammary_gland_cluster_refine_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'mammary_gland_cluster_refine_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Top n motif per group
  topn <- 60
  msplot <- ms[plot==T][order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  # These are top
  msplot <- msplot[msplot [, .I[seq_len(topn)], by = .(gex_mammary_gland_refine_clusters)]$V1]
  msplot <- msplot[!is.na(motifID)]

  # Plots the 'best' motif in grid
  msmat <- as.matrix(dcast.data.table(ms[motifID %in% msplot$motifID], motifID ~ gex_mammary_gland_refine_clusters, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(score)], by=list(motifID, gex_mammary_gland_refine_clusters)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, gex_mammary_gland_refine_clusters)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ gex_mammary_gland_refine_clusters, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(rownames(msmat)),]
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                #geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                add.score.to.name = F,
                                labelMotifClusters= '',
                                labelMotifN = F,
                                #mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Top",topn,"enriched TF motifs per cell type\npadj <",pval.cutoff),
                                title.legend = paste0("Scores"),
                                scaling = 'motif',
                                row_name_fontsize = 13,
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = family_palette)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap_mammary_gland_refine-chromvar-spearman.pdf"),
            width = 6, height = 14, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="gex_mammary_gland_refine_clusters",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    show.group.bar.legend=T,
                                    #scale_data = c(0,1),
                                    title = paste("Top",topn,"enriched TF motifs per cell type"),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = mammary_refine_clusters_palette,
                                    column_name_fontsize = 12,
                                    label.col = family_palette) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-mammary-refine-chromvar.pdf")
  message("Plotting to:",f)
  pdf(file = f, width = 10, height = 19, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$gex_mammary_gland_refine_clusters), function(i){
    toplot <- dt[gex_mammary_gland_refine_clusters==i]
    plotdir <- file.path(plot.dir,gsub("/","_",i))
    dir.create(path = plotdir, recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plotdir,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access_chromvar.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- ChromVar nominate by LHor vs LHor-Arid1a- -------------------------------------------
  library(readxl)
  dtVierstraVarDevFilt <- data.table(read_xlsx(file.path(spreadsheets.dir,
                                                         'archr_multiome_atac_motifs_vierstra.xlsx')))

  clu <- unique(proj$gex_mammary_gland_refine_clusters)
  clu <- clu[which(clu %in% c("L-Hor","L-Hor-Arid1a-"))]
  # for each cell type
  ztests <- lapply(clu, function(x){
    notx <- clu[which(!clu %in% x)]
    xcells <- rownames(proj@cellColData[proj@cellColData$gex_mammary_gland_refine_clusters %in% x,])
    othercells <- rownames(proj@cellColData[proj@cellColData$gex_mammary_gland_refine_clusters %in% notx,])
    zx <- assays(vassay)$z[,xcells]
    zothers <- assays(vassay)$z[,othercells]
    ztest <- mclapply(rownames(zx), function(m){
      test <- wilcox.test(zx[m,], zothers[m,], alternative="greater") # we want sig. up and down and filter away later
      data.table(cluster=x, motifID=m, zdev_p_value = test$p.value )
    }, mc.cores = 8) %>% rbindlist()
    # every motif is tests so we adjust pvalue pvalue
    ztest[,zdev_adjust := p.adjust(p = zdev_p_value, method = "fdr")]
  }) %>% rbindlist()
  setnames(x = ztests, old = "cluster", new = "gex_mammary_gland_refine_clusters")
  # each cluster is tested against all others so we adjust for that test
  ztests[,zdev_adjust := p.adjust(p = zdev_adjust, method = "fdr"), by=gex_mammary_gland_refine_clusters]

  ztests[motifID=="FOXA1_480#FOX/4-forkhead-79-0"]
  ztests[grep("MEIS1", motifID, ignore.case=T)]
  ztests[grep("TEAD2", motifID, ignore.case=T)]

  write.table(x = ztests, file = file.path(data.dir,"ztests_LHor_LHorArid1a.txt"), sep = "\t", row.names = F, col.names = T, quote = T)

  ztests <- fread(file = file.path(data.dir,"ztests_LHor_LHorArid1a.txt"), sep = "\t")

  dtLHor_LArid1a_MotifVierstra <- readxl::read_xlsx(path = file.path(spreadsheets.dir,
                                                              'archr_multiome_atac_motifs_vierstra_enriched_in_LHor_vs_LHorArid1a.xlsx')) %>% setDT()

  lHorcells <- rownames(proj@cellColData[proj@cellColData$gex_mammary_gland_refine_clusters=='L-Hor',])
  lHorArid1acells <- rownames(proj@cellColData[proj@cellColData$gex_mammary_gland_refine_clusters=='L-Hor-Arid1a-',])
  dtMotifGeneCorLHorNoImpute <- correlateAccessibilityExpression(build=T,
                                                                    ArchRproj = proj,
                                                                    dtInputChromVar = dtVierstraVarDevFilt,
                                                                    vassay = vassay[,c(lHorcells, lHorArid1acells)],
                                                                    expassay = expassay[,c(lHorcells, lHorArid1acells)],
                                                                    groupby = 'gex_mammary_gland_refine_clusters',
                                                                    ncores = 12,
                                                                    cor.method = 'spearman',
                                                                    xlsx='dtMotifGeneCorCluster_LHor_LHorArid1a-Spearman.xlsx')

  meisz <- assays(vassay["MEIS1_1031#HD/14-homeodomain-19-0",c(lHorcells)])$z
  meise <- assays(expassay[which(rowData(expassay[,c(lHorcells)])$name=="Meis1"),c(lHorcells)])$GeneExpressionMatrix
  identical( colnames(meise), colnames(meisz) )

  gzero <- which(meise == 0)
  meise <- meise[-gzero]
  meisz <- meisz[-gzero]
  identical( colnames(meise), colnames(meisz) )

  plot(meisz, meise)
  plot(as.numeric(meisz), as.numeric(meise))
  cor.test(x = as.numeric(meisz), y = as.numeric(meise), method = "s")

  dtCor <- dtMotifGeneCorLHorNoImpute
  # This simply combines correlation and chromvar statistics with sample enrichment
  dt <- dtLHor_LArid1a_MotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtLHor_LArid1a_MotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'gex_mammary_gland_refine_clusters', value.name = 'mlog10Padj')
  ms <- merge(dtCor, dt, by=c('motifID','gex_mammary_gland_refine_clusters'))
  ms <- merge(ms, ztests, by=c("motifID","gex_mammary_gland_refine_clusters"))
  ms[,mlog10_zdev_adjust := -log10(max(c(zdev_adjust, .Machine$double.xmin))),by=1:nrow(ms)]
  #ms[,score:= mlog10_zdev_adjust * sign(avg_accessibility) * frac_exp]
  ms[,score:= mlog10_zdev_adjust * avg_accessibility]

  ms <- ms[order(abs(score), mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("-.*","",family)]
  ms[,short_familyID := gsub("/.*","",short_familyID)]
  ms[grep("MEIS1", motifID)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  pval.cutoff <- 0.05
  min.cutoff <- round(-log10(pval.cutoff))
  ms[,plot := F]
  msfilt <- ms[cor_padj < 0.05 &
               mlog10_zdev_adjust >= min.cutoff &
               frac_both > 0.1 &
               sign(avg_accessibility) == sign(avg_expression) &
               coefficient > 0 &
               avg_expression > 0.1]
  msfilt <- msfilt[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  mark <- msfilt[msfilt[, .I[head(top_familyN)], by=list(short_familyID, gex_mammary_gland_refine_clusters)]$V1]
  mark <- msfilt[ridx %in% mark$ridx]
  mark <- mark[order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  setkey(mark, motifID, family, motifGene, geneCorrelated)
  # Unmark duplicate motifGene
  mark <- mark[!duplicated(motifGene)]
  ms[ridx %in% mark$ridx, plot := T]
  nrow(ms[plot==T])
  ms[family == 'FOX/4-forkhead-79' & plot==T,]
  ms[family == 'HD/14-homeodomain-19' & plot==T,]
  ms[family == 'HD/14-homeodomain-19']
  ms[short_familyID %in% c("MIES","HD/14-homeodomain-19") & plot==T]
  ms[grep(pattern = "MEIS1",x = motifID)][plot==T]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'LHor_LHorArid1a_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  #------------------------------------------- Variance Rank Filter Grid by C1 C3 -------------------------------------------
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'LHor_LHorArid1a_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Top n motif per group
  topn <- 60
  msplot <- ms[plot==T][order(score, mlog10_zdev_adjust, avg_accessibility, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  # These are top
  msplot <- msplot[msplot [, .I[seq_len(topn)], by = .(gex_mammary_gland_refine_clusters)]$V1]
  msplot <- msplot[!is.na(motifID)]

  # Plots the 'best' motif in grid
  msmat <- as.matrix(dcast.data.table(ms[motifID %in% msplot$motifID], motifID ~ gex_mammary_gland_refine_clusters, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(score)], by=list(motifID, gex_mammary_gland_refine_clusters)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, gex_mammary_gland_refine_clusters)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ gex_mammary_gland_refine_clusters, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(rownames(msmat)),]
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                #geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                add.score.to.name = F,
                                labelMotifClusters= '',
                                labelMotifN = F,
                                #mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Top",topn,"enriched TF motifs\npadj <",pval.cutoff),
                                title.legend = paste0("Scores"),
                                scaling = 'motif',
                                row_name_fontsize = 13,
                                col = pals::brewer.bupu(5), #solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = family_palette)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Marker-Heatmap_LHor_LHorArid1a-chromvar-spearman.pdf"),
            width = 4, height = 12, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),c(lHorcells, lHorArid1acells)])$z
  htV <- plotSingleCellMotifHeatmap(z,
                                    group.by="gex_mammary_gland_refine_clusters",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    show.group.bar.legend=T,
                                    #scale_data = c(0,1),
                                    title = paste("Top",topn,"enriched TF motifs"),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = mammary_refine_clusters_palette,
                                    column_name_fontsize = 12,
                                    label.col = family_palette) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-_LHor_LHorArid1a-chromvar.pdf")
  message("Plotting to:",f)
  pdf(file = f, width = 8, height = 16, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$Clusters), function(i){
    toplot <- dt[Clusters==i]
    dir.create(path = file.path(plot.dir,i), recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plot.dir,i,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- Correlate TF and Gene by C1 vs C3 -------------------------------------------
  dtVierstraVarDevFilt <- data.table(read_xlsx(file.path(spreadsheets.dir,
                                                         'archr_multiome_atac_motifs_vierstra.xlsx')))
  dtC1_C3_MotifVierstra <- readxl::read_xlsx(path = file.path(spreadsheets.dir,
                                                                     'archr_multiome_atac_motifs_vierstra_enriched_in_C1_vs_C3.xlsx')) %>% setDT()
  dtMotifGeneCorClusterNoImpute <- correlateAccessibilityExpression(build=F,
                                                                    ArchRproj = proj,
                                                                    dtInputChromVar = dtVierstraVarDevFilt,
                                                                    vassay = vassay[,c(c1_cells, c3_cells)], #vassayC1_c3, #vassay,
                                                                    expassay = expassay[,c(c1_cells, c3_cells)], #expassay,
                                                                    min.expression = 0.0,
                                                                    groupby = 'Clusters',
                                                                    ncores = 12,
                                                                    cor.method = 'spearman',
                                                                    xlsx='dtMotifGeneCorCluster_C1_C3Spearman.xlsx')
  dtCor <- dtMotifGeneCorClusterNoImpute
  # This simply combines correlation and chromvar statistics with sample enrichment
  dt <- dtC1_C3_MotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtC1_C3_MotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Clusters', value.name = 'mlog10Padj')
  ms <- merge(dtCor, dt, by=c('motifID','Clusters'))
  ms <- merge(ms, wtests, by="motifID")
  # combine p-values (not used anymore)
  #ms[,padj := exp(-mlog10Padj*log(10))]
  #ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  #ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]
  #ms[,score:= (mlog10Padj + cor_mlog10Padj) * sign(avg_accessibility) * frac_exp]
  ms[,mlog10_zdev_adjust := -log10(max(c(zdev_adjust, .Machine$double.xmin))),by=1:nrow(ms)]
  ms[,score:= mlog10_zdev_adjust * sign(avg_accessibility) * frac_exp]

  ms[,score:= mlog10Padj * sign(avg_accessibility) * frac_exp]

  #ms <- ms[order(abs(score), m10log_fisher_p, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms <- ms[order(abs(score), mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("/.*","",family)]

  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  pval.cutoff <- 0.05
  min.cutoff <- round(-log10(pval.cutoff))
  ms[,plot := F]
  mark <- ms[cor_padj < 0.01 &
               score >= min.cutoff &
               frac_both > 0.1 &
               sign(avg_accessibility) == sign(avg_expression) &
               coefficient > 0]
  mark <- mark[order(score, decreasing = T)][mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID, Clusters)][,ridx]
  #mark <- mark[mark[, .I[head(1)] , by = list(short_motifID, short_familyID)]$V1][,head(.SD, top_familyN),by=list(short_familyID)][,ridx]
  mark <- ms[ridx %in% mark]
  setkey(ms, motifID, family, motifGene, geneCorrelated)
  mark <- lapply(mark$ridx, function(i){
    r <- mark[ridx==i, list(motifID, family, motifGene, geneCorrelated)]
    ms[r, unlist(ridx)]
  }) %>% unlist

  ms[ridx %in% mark, plot := T]
  nrow(ms[plot==T])
  ms[plot==T]
  ms[motifID == 'ANDR_18#FOX/4-forkhead-79-1' & plot==T,]
  ms[motifID == 'ANDR_19#NR/20-nuclearreceptor-248-1' & plot==T,]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'C1_C3_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  #------------------------------------------- Variance Rank Filter Grid by C1 C3 -------------------------------------------
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'C1_C3_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))

  # Plots the 'best' motif in grid
  #(msTop <- sort(unique(ms[plot==T, motifID])))

  msmat <- as.matrix(dcast.data.table(ms[plot==T], motifID ~ Clusters, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0

  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, Clusters)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Clusters)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Clusters, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(ms[plot==T, motifID]),]
  msgroup[is.na(msgroup)] <- ""

  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                cluster.rows = T,
                                labelMotifClusters= '',
                                labelMotifN = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste("Enriched TF motifs in C1 C3\npadj <",pval.cutoff),
                                title.legend = "Enrichment score",
                                scaling = 'none',
                                row_name_fontsize = 12,
                                na.col = 'white',
                                col = solarExtraWhite_palette,
                                column_name_fontsize = 12,
                                label.col = godsnot_102) #pals::alphabet(26)
  cairo_pdf(filename = file.path(plot.dir,paste0("Motif-Heatmap-C1-C3-vertical_Vierstra_top",top_familyN,"_family.pdf")),
            width = 4.2, height = 7, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()

  z <- assays(vassay[rownames(msmat),c(c1_cells, c3_cells)])$z
  range(z)
  quantile(z, seq(0.75,1,by=.01))
  htV <- plotSingleCellMotifHeatmap(enrichMotifs = z,
                                    group.by="Clusters",
                                    orientation_horizontal = F,
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    cluster.rows = T,
                                    labelMotifClusters= '',
                                    labelMotifN = F,
                                    labelMotifFrac = F,
                                    show.group.dendrogram = T,
                                    mlog10Padj.cutOff = 0,
                                    labelMotifs = T,
                                    clustering.within = F,
                                    quantile.min = .10,
                                    quantile.max = .90,
                                    show.group.bar.legend=T,
                                    title = paste("Enriched TF motifs per cluster\npadj <",pval.cutoff),
                                    title.legend = paste0("Z-score"),
                                    row_name_fontsize = 13,
                                    na.col = 'gray90',
                                    col = solarExtraWhite_palette, # rev(brewer.ylgnbu(25)),
                                    group.colors = cluster_palette,
                                    column_name_fontsize = 12,
                                    label.col = godsnot_102) #rev(pals::cubehelix())
  f <- file.path(plot.dir,"Motif-SCHeatmap-C1_C3-vertical_Vierstra_family.png")
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 8*res, height = 16*res, res = res, bg = 'transparent')
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
  })
  dev.off()
  image_write(image = fig, path = f, format = 'png')
  pdf(file = file.path(plot.dir,"Motif-SCHeatmap-C1_C3-vertical_Vierstra_family.pdf"), width = 8, height = 16)
  htV <- draw(htV, heatmap_legend_side = "bottom", padding = unit(c(8, 2, 2, 2), "mm")) #bltr
  decorate_heatmap_body(names(htV), { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1)) })
  dev.off()

  #-------------------------------------------  Plot UMAPS of top motifs -------------------------------------------#
  dt <- ms[ms[plot==T  & sign(avg_accessibility) > 0, .I[score == max(score)], by=list(motifID)]$V1]

  mclapply(unique(dt$Clusters), function(i){
    toplot <- dt[Clusters==i]
    dir.create(path = file.path(plot.dir,i), recursive = T, showWarnings = F)
    for(j in 1:nrow(toplot)){
      a <- toplot[j,motifGene]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal = solarExtraWhite_palette,
                                              pal2 = solarHalf_palette, #gist_heat_r,\
                                              pal3 = rocket_palette,
                                              pt.size = 1,
                                              legend.height.factor = 0.5,
                                              imputeValues = F,
                                              plot.title = "",
                                              nrows = nrow)
      png(filename = file.path(plot.dir,i,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png')),
          width = (lengths(gpGene)*4.6)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
    }
  }, mc.cores = 16)

  #------------------------------------------- Correlate TF and Gene by C1 vs C2 -------------------------------------------
  dtMotifGeneCorClusterNoImpute <- correlateAccessibilityExpression(build=F,
                                                                    ArchRproj = proj,
                                                                    dtInputChromVar = dtVierstraVarDevFilt,
                                                                    vassay = vassay,
                                                                    expassay = expassay,
                                                                    #matGE = matGE,
                                                                    min.expression = 0.0,
                                                                    groupby = 'Clusters',
                                                                    ncores = 16,
                                                                    cor.method = 'spearman',
                                                                    xlsx='dtMotifGeneCorClusterNoImpute.xlsx')

  dtKO_vs_LumProgenitorWT_MarkersPeaks <- read_xlsx(path = file.path(spreadsheets.dir,
                                                                 'archr_multiome_atac_motifs_vierstra_enriched_in_C1_vs_C2.xlsx')) %>% setDT

  dt <- dtKO_vs_LumProgenitorWT_MarkersPeaks[,c("motifID",grep("mlog10Padj", colnames(dtKO_vs_LumProgenitorWT_MarkersPeaks), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Clusters', value.name = 'mlog10Padj')
  ms <- merge(dtMotifGeneCorClusterNoImpute, dt, by=c('motifID','Clusters'))
  ms[,padj := exp(-mlog10Padj*log(10))]
  ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]
  ms[,stouffer_p := stouffer(x = c(cor_padj, padj))$p, by=1:nrow(ms) ]
  ms[,m10log_stouffer := -log10(max(c(stouffer_p, .Machine$double.xmin))), by=1:nrow(ms)]
  ms <- ms[order(mlog10Padj, avg_accessibility, rho, decreasing = T)]
  ms <- unique(ms)
  ms[,ridx := 1:nrow(ms)]
  ms[,percent_both := 100*(ncell_both / ncell_motif)]
  ms[,plot := F]

  # set pval threshold
  minlog10pval <- -log10(0.05) # 10^(log10(0.01)), 10
  percent_both_cutoff <- quantile(ms[, percent_both], probs=.25)
  # we want significant motifs, that are significantly positively correlated with expression.
  ms[ridx %in% ms[cor_padj < 0.05 &
                    padj < 0.05 &
                    rho > 0 &
                    avg_accessibility > 0 &
                    percent_both >= percent_both_cutoff &
                    mlog10Padj >= minlog10pval, head(.SD, 1), by=list(family, Clusters)][,ridx], plot := T]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'C1_C2_Motif_Heatmap_Vierstra_Filt.xlsx'), format_headers = T, col_names = T)

  #------------------------------------------- Variance Rank Filter Grid by C1 C2 -------------------------------------------
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'C1_C2_Motif_Heatmap_Vierstra_Filt.xlsx')))
  msmat <- as.matrix(dcast.data.table(ms[!is.na(avg_accessibility)], motifID ~ Clusters, value.var = 'mlog10Padj', fun.aggregate = max), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0
  msmat <- msmat[unique(ms[plot==T,motifID]),]
  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, Clusters)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Clusters)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Clusters, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(ms[plot==T, motifID]),]
  msgroup[is.na(msgroup)] <- ""

  p.val <- 0.05
  mlog10.cutoff <- -log10(p.val)
  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htH <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                labelMotifClusters= '',
                                cutOff = mlog10.cutoff,
                                labelMotifs = T,
                                title = "Enriched transcription factor motifs",
                                title.legend = paste("p-value <",p.val),
                                scaling = 'motif',
                                row_name_fontsize = 10,
                                #col = gist_heat_r,
                                column_name_fontsize = 12)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-C1_vs_C2-vertical-Vierstra-Filt-NoImpute-mlog10.pdf"),
            width = 4.5, height = 9.5, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #tlbr
  dev.off()

  # Plot individual motifs on UMAPS
  dt <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'C1_C2_Motif_Heatmap_Vierstra_Filt.xlsx')))[plot==T]
  for(i in unique(dt$Clusters)) {
    toplot <- dt[Clusters==i]
    for(j in 1:nrow(toplot)){
      a <- toplot[j,gsub("#.*","",motifGene)]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal2 = solarHalf, #gist_heat_r,\
                                              pal3 = solarHalf,
                                              imputeValues = F,
                                              plot.geneScore = F,
                                              plot.title = F,
                                              nrows = nrow)
      out.dir <- ifelse(i=='C1', file.path(plot.dir,"C1_C2"), file.path(plot.dir,"C2_C1"))
      dir.create(path = out.dir, recursive = T, showWarnings = F)
      out.file <- file.path(out.dir, paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png'))
      png(filename = out.file,
          width = (lengths(gpGene)*4)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
      message("Plotted ",out.file)
    }
  }

  #------------------------------------------- Correlate TF and Gene by C2 vs C3 -------------------------------------------
  dtMotifGeneCorClusterNoImpute <- correlateAccessibilityExpression(build=F,
                                                                    ArchRproj = proj,
                                                                    dtInputChromVar = dtVierstraVarDevFilt,
                                                                    vassay = vassay,
                                                                    expassay = expassay,
                                                                    #matGE = matGE,
                                                                    min.expression = 0.0,
                                                                    groupby = 'Clusters',
                                                                    ncores = 16,
                                                                    cor.method = 'spearman',
                                                                    xlsx='dtMotifGeneCorClusterNoImpute.xlsx')

  dtE2MatureWTvsMatureWT_MarkersPeaks <- read_xlsx(path = file.path(spreadsheets.dir,
                                                                 'archr_multiome_atac_motifs_vierstra_enriched_in_C2_vs_C3.xlsx')) %>% setDT

  dt <- dtE2MatureWTvsMatureWT_MarkersPeaks[,c("motifID",grep("mlog10Padj", colnames(dtE2MatureWTvsMatureWT_MarkersPeaks), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Clusters', value.name = 'mlog10Padj')
  ms <- merge(dtMotifGeneCorClusterNoImpute, dt, by=c('motifID','Clusters'))
  ms[,padj := exp(-mlog10Padj*log(10))]
  ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]
  ms[,stouffer_p := stouffer(x = c(cor_padj, padj))$p, by=1:nrow(ms) ]
  ms[,m10log_stouffer := -log10(max(c(stouffer_p, .Machine$double.xmin))), by=1:nrow(ms)]
  ms <- ms[order(mlog10Padj, avg_accessibility, -cor_padj, rho, decreasing = T)]
  ms[,plot := F]
  ms <- unique(ms)
  ms[,percent_both := 100*(ncell_both / ncell_motif)]
  ms[,ridx := 1:nrow(ms)]

  # set pval threshold
  minlog10pval <- -log10(0.05) # 10^(log10(0.01)), 10
  percent_both_cutoff <- quantile(ms[, percent_both], probs=.25)
  # we want significant motifs, that are significantly positively correlated with expression.
  ms[ridx %in%
       ms[cor_padj < 1 &
            padj < 0.05 &
            rho > 0 &
            avg_accessibility > 0 &
            percent_both >= percent_both_cutoff &
            mlog10Padj >= minlog10pval, head(.SD, 1), by=list(family, Clusters)][,ridx], plot := T]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'C2_C3_Motif_Heatmap_Vierstra_Filt.xlsx'), format_headers = T, col_names = T)
  #------------------------------------------- Variance Rank Filter Grid by C2 C3 -------------------------------------------
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'C2_C3_Motif_Heatmap_Vierstra_Filt.xlsx')))
  msmat <- as.matrix(dcast.data.table(ms[!is.na(avg_accessibility)], motifID ~ Clusters, value.var = 'mlog10Padj', fun.aggregate = max), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0
  msmat <- msmat[unique(ms[plot==T,motifID]),]
  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, Clusters)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Clusters)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Clusters, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(ms[plot==T, motifID]),]
  msgroup[is.na(msgroup)] <- ""

  p.val <- 10^-10
  mlog10.cutoff <- -log10(p.val)
  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htH <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                labelMotifClusters= '',
                                cutOff = mlog10.cutoff,
                                labelMotifs = T,
                                title = "Enriched transcription factor motifs",
                                title.legend = paste("p-value <",p.val),
                                scaling = 'motif',
                                row_name_fontsize = 10,
                                #col = gist_heat_r,
                                column_name_fontsize = 12)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-C2_vs_C3-vertical-Vierstra-Filt-NoImpute-mlog10.pdf"),
            width = 4, height = 9.5, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #tlbr
  dev.off()

  # Plot individual motifs on UMAPS
  dt <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'C2_C3_Motif_Heatmap_Vierstra_Filt.xlsx')))[plot==T]
  for(i in unique(dt$Clusters)) {
    toplot <- dt[Clusters==i]
    for(j in 1:nrow(toplot)){
      a <- toplot[j,gsub("#.*","",motifGene)]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal2 = solarHalf, #gist_heat_r,\
                                              pal3 = solarHalf,
                                              imputeValues = F,
                                              plot.geneScore = F,
                                              plot.title = F,
                                              nrows = nrow)
      out.dir <- ifelse(i=='C2', file.path(plot.dir,"C2_C3"), file.path(plot.dir,"C3_C2"))
      dir.create(path = out.dir, recursive = T, showWarnings = F)
      out.file <- file.path(out.dir, paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png'))
      png(filename = out.file,
          width = (lengths(gpGene)*4)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
      message("Plotted ",out.file)
    }
  }

  #------------------------------------------- Correlate TF and Gene by C1 vs C2 C3 -------------------------------------------
  dtMotifGeneCorClusterNoImpute <- correlateAccessibilityExpression(build=F,
                                                                    ArchRproj = proj,
                                                                    dtInputChromVar = dtVierstraVarDevFilt,
                                                                    vassay = vassay,
                                                                    expassay = expassay,
                                                                    #matGE = matGE,
                                                                    min.expression = 0.0,
                                                                    groupby = 'Clusters',
                                                                    ncores = 16,
                                                                    cor.method = 'spearman',
                                                                    xlsx='dtMotifGeneCorClusterNoImpute.xlsx')

  dtKO_vs_LumWT_MarkersPeaks <- read_xlsx(path = file.path(spreadsheets.dir,
                                                                 'archr_multiome_atac_motifs_vierstra_enriched_in_C1_vs_C2C3.xlsx')) %>% setDT

  dt <- dtKO_vs_LumWT_MarkersPeaks[,c("motifID",grep("mlog10Padj", colnames(dtKO_vs_LumWT_MarkersPeaks), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Clusters', value.name = 'mlog10Padj')
  ms <- merge(dtMotifGeneCorClusterNoImpute, dt, by=c('motifID','Clusters'))
  ms[,padj := exp(-mlog10Padj*log(10))]
  ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]
  ms[,stouffer_p := stouffer(x = c(cor_padj, padj))$p, by=1:nrow(ms) ]
  ms[,m10log_stouffer := -log10(max(c(stouffer_p, .Machine$double.xmin))), by=1:nrow(ms)]
  ms <- ms[order(mlog10Padj, avg_accessibility, rho, decreasing = T)]
  ms[,plot := F]
  ms <- unique(ms)
  ms[,ridx := 1:nrow(ms)]

  # set pval threshold
  minlog10pval <- -log10(0.05) # 10^(log10(0.01)), 10
  # we want significant motifs, that are significantly positively correlated with expression.
  ms[ridx %in% ms[cor_padj < 0.05 &
                    padj < 0.05 &
                    rho > 0 &
                    avg_accessibility > 0 &
                    mlog10Padj >= minlog10pval, head(.SD, 1), by=list(family, Clusters)][,ridx], plot := T]

  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'C1_C2C3_Motif_Heatmap_Vierstra_Filt.xlsx'), format_headers = T, col_names = T)
  #------------------------------------------- Variance Rank Filter Grid by C1 C3 -------------------------------------------
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'C1_C2C3_Motif_Heatmap_Vierstra_Filt.xlsx')))
  msmat <- as.matrix(dcast.data.table(ms[!is.na(avg_accessibility)], motifID ~ Clusters, value.var = 'mlog10Padj', fun.aggregate = max), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0
  msmat <- msmat[unique(ms[plot==T,motifID]),]
  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, Clusters)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Clusters)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Clusters, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(ms[plot==T, motifID]),]
  msgroup[is.na(msgroup)] <- ""

  p.val <- 0.05
  mlog10.cutoff <- -log10(p.val)
  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))

  htH <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T,
                                labelMotifClusters= '',
                                cutOff = mlog10.cutoff,
                                labelMotifs = T,
                                title = "Enriched transcription factor motifs",
                                title.legend = paste("p-value <",p.val),
                                scaling = 'motif',
                                row_name_fontsize = 10,
                                #col = gist_heat_r,
                                column_name_fontsize = 12)
  cairo_pdf(filename = file.path(plot.dir,"Motif-Enrichment-Marker-Heatmap-C1_vs_C2C3-vertical-Vierstra-Filt-NoImpute-mlog10.pdf"),
            width = 4.5, height = 9.5, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htH, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #tlbr
  dev.off()

  # Plot individual motifs on UMAPS
  dt <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'C1_C2C3_Motif_Heatmap_Vierstra_Filt.xlsx')))[plot==T]
  for(i in unique(dt$Clusters)) {
    toplot <- dt[Clusters==i]
    for(j in 1:nrow(toplot)){
      a <- toplot[j,gsub("#.*","",motifGene)]
      g <- toplot[j,geneCorrelated]
      nrow <- 1
      gpGene <- PlotMotifExpressionGScoreUmap(proj = proj,
                                              DT = toplot[j, list(motifID, gene=geneCorrelated)],
                                              pal2 = solarHalf, #gist_heat_r,\
                                              pal3 = solarHalf,
                                              imputeValues = F,
                                              plot.geneScore = F,
                                              plot.title = F,
                                              nrows = nrow)
      out.dir <- ifelse(i=='C1', file.path(plot.dir,"C1_C2C3"), file.path(plot.dir,"C2C3_C1"))
      dir.create(path = out.dir, recursive = T, showWarnings = F)
      out.file <- file.path(out.dir, paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.png'))
      png(filename = out.file,
          width = (lengths(gpGene)*4)/nrow,
          height = 4*nrow, bg = 'transparent', units = 'in', res = 200)
      do.call(grid.arrange, c(gpGene, nrow=1))
      dev.off()
      message("Plotted ",out.file)
    }
  }

  ### Scatter plots
  names(cluster_palette) <- paste0('C',names(cluster_palette))

  g <- lapply(unique(ms$Clusters), function(m){
    toplot <- ms[Clusters==m]
    minlog10pval <- head(toplot[order(mlog10Padj, decreasing = T)], 50)[,min(mlog10Padj)]
    toplot[,label:=ifelse(mlog10Padj < minlog10pval, "", paste(gsub("#.*","",motifID), geneCorrelated, sep='/'))]
    ggplot(toplot, aes(x=rho, y=mlog10Padj, label=label)) +
      geom_vline(xintercept = 0, linetype='dashed') +
      geom_text_repel(
        data          = subset(toplot, mlog10Padj >= minlog10pval & rho > 0), #subset(toplot, rho > 0.2 & mlog10Padj > 2),
        nudge_x       = .3 + subset(toplot, mlog10Padj >= minlog10pval & rho > 0)$rho,
        seed = 3,
        segment.size  = 0.2,
        color='blue',
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        size         = 2,
        show.legend = F,
        max.overlaps = 20,
      ) +
      geom_text_repel(
        data          = subset(toplot, mlog10Padj >= minlog10pval & rho <= 0), #subset(toplot, rho < -0.2 & mlog10Padj > 2),
        nudge_x       = -.1 - subset(toplot, mlog10Padj >= minlog10pval & rho <= 0)$rho, # + subset(toplot, rho < 0)$rho,
        seed = 3,
        segment.size  = 0.2,
        color='blue',
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        size          = 2,
        max.overlaps = 20,
        show.legend = F
      ) +
      geom_hline(yintercept = -log(0.01)/log(10), col='red', linetype='dotted') +
      geom_point(aes(color=Clusters), shape=20, size=.8, show.legend = F) +
      #scale_color_manual("", values = cluster_palette) +
      facet_wrap(~Clusters, scales = 'free_y', nrow = 1) +
      scale_y_continuous(breaks = seq(from = 0, to = max(ceiling(toplot$mlog10Padj)), by=100)) +
      ylab("-log10 adjusted pvalues (Accessibility)") +
      xlab("pearson's r") +
      theme_bw() + theme(panel.grid = element_blank())
  })
  png(filename = file.path(plot.dir,'Access_clusters_expression_accessibililty.png'), width = 35, height = 16,
      bg = 'transparent', res = 150, units = 'in')
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  g <- lapply(unique(ms$Clusters), function(m){
    toplot <- ms[Clusters==m]
    toplot[,label:=ifelse(abs(rho) < 0.2, "", paste(gsub("#.*","",motifID),geneCorrelated, sep='/'))]
    ggplot(toplot, aes(x=rho, y=cor_mlog10Padj, label=label)) +
      geom_vline(xintercept = 0, linetype='dashed') +
      geom_text_repel(
        data          = subset(toplot, rho > 0.2 & cor_mlog10Padj > 2),
        nudge_x       = .3 + subset(toplot, rho > 0.2 & cor_mlog10Padj > 2)$rho,
        seed = 3,
        segment.size  = 0.2,
        color='blue',
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        size         = 2,
        show.legend = F,
        max.overlaps = 20,
      ) +
      geom_text_repel(
        data          = subset(toplot, rho < -0.2 & cor_mlog10Padj > 2),
        nudge_x       = -.1 - subset(toplot, rho < -0.2 & cor_mlog10Padj > 2)$rho,
        seed = 3,
        segment.size  = 0.2,
        color='blue',
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        size          = 2,
        max.overlaps = 20,
        show.legend = F
      ) +
      geom_hline(yintercept = -log(0.01)/log(10), col='red', linetype='dotted') +
      geom_point(aes(color=Clusters), shape=20, size=.8, show.legend = F) +
      scale_color_manual("", values = cluster_palette) +
      facet_wrap(~Clusters, scales = 'free_y', nrow = 1) +
      scale_y_continuous(breaks = seq(from = 0, to = max(ceiling(toplot$cor_mlog10Padj)), by=10)) +
      ylab("-log10 adjusted pvalues (Correlation)") +
      xlab("pearson's r") +
      theme_bw() + theme(panel.grid = element_blank())
  })
  png(filename = file.path(plot.dir,'Access_clusters_expression_cor.png'), width = 25, height = 8,
      bg = 'transparent', res = 150, units = 'in')
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  g <- lapply(unique(ms$Clusters), function(m){
    #toplot <- ms[bioNames==m & cor_padj < 0.05]
    toplot <- ms[Clusters==m]
    toplot[,label:=""]
    toplot <- toplot[order(m10log_stouffer, decreasing = T)]
    toplot[,idx:=1:nrow(toplot)]

    tolabel <- toplot[rho > .2 & cor_padj < 0.05, head(.SD, 1), by=list(family)]$idx

    toplot[idx %in% tolabel, label:=paste(gsub("#.*","",motifID),geneCorrelated, sep='/')]

    tolabel <- toplot[rho < -.2 & cor_padj < 0.05, head(.SD, 1), by=list(family)]$idx

    toplot[idx %in% tolabel, label:=paste(gsub("#.*","",motifID),geneCorrelated, sep='/')]

    ggplot(toplot, aes(x=rho, y=m10log_stouffer, label=label)) +
      geom_vline(xintercept = 0, linetype='dashed') +
      geom_vline(xintercept = .2, linetype='dotted', color='grey70') +
      geom_vline(xintercept = -.2, linetype='dotted', color='grey70') +
      geom_hline(yintercept = -log(0.01)/log(10), col='red', linetype='dotted') +
      geom_point(data = toplot[cor_padj>0.05], color='gray80', shape=20, size=.8, show.legend = F) +
      geom_point(data = toplot[cor_padj<0.05], aes(color=Clusters), shape=20, size=.8, show.legend = F) +
      geom_text_repel(
        data          = subset(toplot, rho > 0.2 & m10log_stouffer > 2),
        nudge_x       = .5 - subset(toplot, rho > 0.2 & m10log_stouffer > 2)$rho,
        force        = 0.5,
        seed = 3,
        segment.size  = 0.2,
        color='black',
        segment.color = "grey50",
        segment.curvature = -1e-20,
        direction     = "y",
        hjust         = 0,
        size         = 3,
        show.legend = F,
        max.overlaps = 15,
      ) +
      geom_text_repel(
        data          = subset(toplot, rho < -0.2 & m10log_stouffer > 2),
        nudge_x       = -.5 - subset(toplot, rho < -0.2 & m10log_stouffer > 2)$rho,
        force        = 0.5,
        seed = 3,
        segment.size  = 0.2,
        color='black',
        segment.color = "grey50",
        segment.curvature = -1e-20,
        direction     = "y",
        hjust         = 1,
        size          = 3,
        max.overlaps = 15,
        show.legend = F
      ) +
      scale_color_manual("", values = cluster_palette) +
      facet_wrap(~Clusters, scales = 'free_y', nrow = 1) +
      scale_x_continuous(limits = c(-1,1)) +
      scale_y_continuous(breaks = seq(from = 0, to = max(ceiling(toplot$m10log_stouffer)), by=50)) +
      ylab("-log10 adjusted pvalues (combined)") +
      xlab("pearson's r") +
      theme_bw() + theme(panel.grid = element_blank())
  })
  png(filename = file.path(plot.dir,'Access_clusters_expression_corr_stouffer.png'), width = 25, height = 8,
      #png(filename = file.path(plot.dir,'Access_bionames_expression_corr_s.png'), width = 25, height = 8,
      bg = 'transparent', res = 150, units = 'in')
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  #------------------------------------------- Plot UMAP Accessiblity + Expression -------------------------------------------

  markerRNA <- c("Smarca4")
  dtMotifs <- proj@VierstraMotifToGene
  motifId <- grep("Smarca4|Brg1", dtMotifs$motifID, ignore.case = T, value = T)
  pmg <- PlotMotifExpressionGScoreUmap(proj = proj,
                                       DT=data.table(motifID=motifId, gene=markerRNA),
                                       group = T,
                                       nrows = 1)
  a <- "A"
  g <- 'G'
  nrow <- 4
  pdf(file = file.path(plot.dir,paste0('motif_',a,'_gene_',g,'_Umap_Motif_Exp_Access.pdf')),
      width = (sum(lengths(pmg))*3.4)/nrow,
      height = 3*nrow,
      bg = 'transparent')
      do.call(grid.arrange, c(pmg, nrow=nrow))
  dev.off()

  #------------------------------------------- Correlate accessibility and pseudotime -------------------------------------------
  dtVierstraVarDevFilt <- setDT(read_xlsx(file.path(spreadsheets.dir, 'archr_multiome_atac_motifs_vierstra.xlsx')))
  dtCorPseudotimeMotif <- correlateAccessibilityPseudotime(build=F,
                                                           ArchRproj = proj,
                                                           dtInputChromVar = dtVierstraVarDevFilt,
                                                           vassay = vassay,
                                                           groupby = 'gex_mammary_gland_refine_clusters',
                                                           ncores = 16,
                                                           cor.method = 'spearman',
                                                           xlsx='dtMotifPseudotimeMammaryRefineSpearman.xlsx')
  dtCorPseudotimeMotif <- unique(dtCorPseudotimeMotif[motifGene==geneCorrelated])
  dtCorPseudotimeMotif[,ridx := 1:nrow(dtCorPseudotimeMotif)]
  dtCorPseudotimeMotif[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  dtCorPseudotimeMotif[,short_familyID := gsub("-.*","",family)]
  dtCorPseudotimeMotif[,short_familyID := gsub("/.*","",short_familyID)]

  n_per_family <- 1
  mark <- dtCorPseudotimeMotif[dtCorPseudotimeMotif[cor_padj < 0.05 & avg_accessibility > 0, #& coefficient > 0,
                                                    .I[head( n_per_family )],
                                                    by=list(short_familyID, gex_mammary_gland_refine_clusters)]$V1]
  dtfilt <- dtCorPseudotimeMotif[ridx %in% mark$ridx]
  dtfilt <- dtfilt[order(coefficient, decreasing=T),.SD, by=gex_mammary_gland_refine_clusters]
  dtfilt[,rowid := 1:.N, by=gex_mammary_gland_refine_clusters]

  dtfilt[between(x = rowid, 1, 20)]

  gp <- lapply(unique(dtfilt$gex_mammary_gland_refine_clusters), function(i){
    ggplot(dtfilt[gex_mammary_gland_refine_clusters==i][between(x = rowid, 1, 25)],
           aes(y=reorder(short_motifID, -rowid), x=coefficient,
               #color=gex_mammary_gland_refine_clusters,
               fill=avg_accessibility)) +
      geom_point(shape=1, size=3, color='black') +
      geom_point(size=3, shape=21) +
      facet_wrap( ~ gex_mammary_gland_refine_clusters, scales = "free_y") +
      ylab("TF motifs") +
      scale_fill_gradientn("µ motif \naccessibility", colours = pals::brewer.bupu(20)) +
      theme_bw() + theme(panel.grid = element_line(linewidth = 0.1),
                         legend.title = element_text(size=8),
                         legend.key.width = unit(dev.size()[1] / 80, "inches"),
                         strip.text = element_text(size=16),
                         strip.background = element_rect(fill = mammary_refine_clusters_palette_list[[i]])) +
      guides(fill = guide_colorbar(order=1, draw.ulim = T, draw.llim = T))
  })
  pdf(file = file.path(plot.dir, "correlated_motif_accessibility_pseudotime.pdf"), width = 11, height = 7.5, bg = "transparent")
  do.call(grid.arrange, c(gp, nrow=2))
  dev.off()

  TFfit <- function(signal, pseudotime){
    signal_norm <- signal
    st <- stats::loess(signal_norm ~ pseudotime)
    st$fitted
  }

  TFsmooth <- function(signal, pseudotime, fit="gam"){
    y <- NULL
    pseudotime_equidistant <- seq(from=min(pseudotime), to=max(pseudotime),length.out=length(pseudotime))
    pseudotime_equidistant <- seq(from=0, to=1, length.out=length(pseudotime))
    library(mgcv)
    if(fit=="cellrank"){
      message("Using fit of cellrank")
      # default fit from cellrank
      fit <- gam(signal ~ s(pseudotime, bs="cr", k=5), family=gaussian)
      y <- predict(fit, data.frame(pseudotime=pseudotime_equidistant))
    } else if(fit=="gam"){
      message("Using GAM fit with linear extrapolation")
      # adjusted fit that uses linear model for extrapolation
      # (i.e. ends without sufficient data are linearly increasing or decreasing)
      # fit <- gam(signal ~ ns(pseudotime, df=6), family=gaussian)
      fit <- gam(signal ~ ns(pseudotime, df=5), family=gaussian)
      y <- predict(fit, data.frame(pseudotime=pseudotime_equidistant))
    } else {
      message("Using loess fit with extrapolation")
      st <- stats::loess(signal ~ pseudotime)
      df$cells <- rownames(df)
      # surface = direct "extrapolate values" for points out of range
      y <- predict(object = st, pseudotime_equidistant, surface="direct")
    }
    names(y) <- pseudotime_equidistant
    y
  }

  for(i in unique(dtfilt$gex_mammary_gland_refine_clusters)){
    dtplot <- dtfilt[gex_mammary_gland_refine_clusters==i][between(x = rowid, 1, 25)]
    unique(dtplot$motifID)
    cells <- proj@cellColData[proj@cellColData$gex_mammary_gland_refine_clusters==i,]
    cells <- cells[order(cells$palantir_pseudotime, decreasing = F),]
    pseudotimes <- as.numeric(cells$palantir_pseudotime)
    plot(pseudotimes, main=i)
    cells <- rownames(cells)
    z <- assays(vassay[dtplot$motifID,cells])$z

    ### Test the fits ###
    if(TRUE){
      gp <- lapply(rownames(z), function(motif){
        df <- data.frame(signal=z[motif,])
        df$cells <- rownames(df)
        df$pseudotime <- as.numeric(proj@cellColData[colnames(z), "palantir_pseudotime"])

        library(mgcv)
        library(splines)
        #sp="smoothing_penalty=1.0",
        #fit <- gam(df$signal ~ s(df$pseudotime, bs="cr", k=7, fx =FALSE))
        st <- stats::loess(signal ~ pseudotime, data = df)
        dfpredict <- data.frame(pseudotime = seq(from=0, to=1,length.out=length(df$pseudotime)))
        dfpredict$loess <- predict(object = st, dfpredict$pseudotime, surface="direct")

        # cellrank
        fit <- gam(signal ~ s(pseudotime, bs="cr", k=5), family=gaussian, data = df)
        dfpredict$gam1 <- predict(fit, dfpredict)
        plot(fit)

        # cubic spline
        #fit <- gam(signal ~ ns(pseudotime, df=6), family=gaussian, data = df)
        fit <- gam(signal ~ ns(pseudotime, df=5), family=gaussian, data = df)
        dfpredict$gam2 <- predict(fit, dfpredict)

        gp <- ggplot(df, aes(x=pseudotime, y=signal)) +
          geom_point(size = .1) +
          geom_rug(linewidth = .1) +
          geom_point(data = dfpredict, aes(x=pseudotime, y=loess), color='red', size=1) +
          geom_point(data = dfpredict, aes(x=pseudotime, y=gam1), color='green4', size=1) +
          geom_point(data = dfpredict, aes(x=pseudotime, y=gam2), color='darkblue', size=1) +
          geom_smooth() +
          ggtitle(motif) + theme_bw(16) + theme(panel.grid = element_line() )
        gp
      }) %>% set_names(gsub("#.*","",rownames(z)))
      i <- gsub("\\/","_",i)
      png(file=file.path(plot.dir,paste0("correlated_motif_accessibility_pseudotime_curves_",i,".png")),
          width=35, height=35, units = "in", res = 200)
      do.call(grid.arrange, c(gp, nrow=5))
      dev.off()
    }

    ##

    ### apply the fits to the matrix
    pseudotimes <- as.numeric(proj@cellColData[colnames(z), "palantir_pseudotime"])

    pseudotime_equidistant <- seq(from=0, to=1, length.out=length(pseudotimes))
    zs <- apply(z, 1, function(x) {
      vals <- TFsmooth(signal = x, pseudotime = pseudotimes, fit = "gam")
      #vals <- TFsmooth(signal = x, pseudotime = pseudotimes, fit = "cellrank")
    }) %>% t()

    colnames(zs) <- colnames(z)[seq(ncol(zs))]

    htV <- plotSingleCellMotifHeatmap(zs,
                                      cellset=colnames(zs),
                                      cellset_ordered = T,
                                      pseudotimes = pseudotime_equidistant[1:dim(zs)[2]],
                                      pseuodtime_cutoff = c(0, 1),
                                      pseudotime_palette = pals::brewer.bupu,
                                      scale_data = c(0,1),
                                      group.by="gex_mammary_gland_refine_clusters",
                                      orientation_horizontal = F,
                                      remove.redundant = T,
                                      cluster.columns = F,
                                      cluster.rows = T,
                                      labelMotifClusters= '',
                                      labelMotifN = F,
                                      labelMotifFrac = F,
                                      show.group.dendrogram = T,
                                      mlog10Padj.cutOff = 0,
                                      labelMotifs = T,
                                      clustering.within = F,
                                      show.group.bar.legend=T,
                                      title = paste("Top enriched TF motifs in",i),
                                      title.legend = paste0("Z-score"),
                                      row_name_fontsize = 13,
                                      na.col = 'gray90',
                                      col = brewer.blues(100), # rev(brewer.ylgnbu(25)),
                                      group.colors = mammary_refine_clusters_palette,
                                      column_name_fontsize = 12,
                                      use.raster = T,
                                      label.col = family_palette) #rev(pals::cubehelix())
    i <- gsub("\\/","_",i)
    width <- 9; height <- 5.5
    cairo_pdf(filename = file.path(plot.dir,paste0("correlated_motif_accessibility_pseudotime_heatmap_",i,".pdf")),
              width = width, height = height, bg = 'transparent', onefile = T, family = 'Helvetica')
    draw(htV, heatmap_legend_side = "right", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
    dev.off()

  }

  #------------------------------------------- Plot Motifs -------------------------------------------
  ### order the motifs by largest deviation z-score (no quantile cut)

  (motifIDs <- grep( "FOXC1_492" , names(proj@peakAnnotation$Vierstra$motifs), value=T, ignore.case = T))

  gp <- lapply(motifIDs, function(motifId){
    (smotifId <- gsub("#.*","",motifId))
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = NULL,
                    feature = paste0('z:',motifId),
                    colorBy = 'VierstraMatrix',
                    pt.size = .5,
                    addLabels = F,
                    panel.border.size = 0,
                    embedding = 'UMAP',
                    title = smotifId,
                    label.size = 14,
                    title.size = 18,
                    break.unit = NA,
                    panel_width_in = 3,
                    scale_color_bar = T,
                    pal = solarExtraWhite_palette,
                    legend.title = 'Motif deviation zscore',
                    legend.title.position = "right",
                    legend.height.factor = .6,
                    filename.prefix = paste0('Umap_zDeviation_',smotifId))
  })

  gp[[1]]
  gp[[2]]

  x <- grep( "Arid3a_32#HD/18-homeodomain-24-1" , names(proj@peakAnnotation$Vierstra$motifs), value=T, ignore.case = T)

  gp <- lapply(x, function(motifId){
    smotifId <- gsub("#.*","",motifId)
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = NULL,
                    feature = paste0('z:',motifId),
                    colorBy = 'VierstraMatrix',
                    pt.size = 1,
                    addLabels = F,
                    panel.border.size = 0,
                    embedding = 'UMAP',
                    title = smotifId,
                    label.size = 14,
                    title.size = 18,
                    break.unit = NA,
                    panel_width_in = 3,
                    scale_color_bar = T,
                    pal = solarExtraWhite_palette,
                    legend.title = 'Motif deviation zscore',
                    legend.title.position = "right",
                    legend.height.factor = .6)
  })
  length(gp)
  do.call(grid.arrange, c(gp, nrow=3))

  #------------------------------------------- Variance Rank Filter Plot -------------------------------------------
  # After viewing top deviations and grid of top motifs by family we filter "retained" motifs for
  # 1) best in family and 2) most correlated with expression

  dtMotifGeneCorTopCorr <- read_xlsx(path = file.path(spreadsheets.dir,'Clusters_Motif_Heatmap_Vierstra_Filt.xlsx')) %>% setDT
  dtMotifGeneCorTopCorr <- dtMotifGeneCorTopCorr[plot==T]
  dtMotifGeneCorTopCorr[,rank:=seq_along(1:.N)]
  #dtMotifGeneCorTopCorr <- dtMotifGeneCorTopCorr[order(rank)][, head(.SD, 1), by=motifID] # top by rank
  g <- PlotRankCorrelations(dt = dtMotifGeneCorTopCorr, n = 50, text.adjust = T)

  plotVierstraVarDev <- getVarDeviations(proj, name = "VierstraMatrix", plot = TRUE, n = 100)
  dtVierstraVarDevFilt <- data.table(plotVierstraVarDev$data) %>% setnames(., 'name','motifID')
  # Plots Top deviations
  (PlotDeviations(dt = dtVierstraVarDevFilt, n=75, text.adjust = T))

  ntop <- 50
  cairo_pdf(filename = file.path(plot.dir, "Variable-Vierstra-Motif-Deviation-Scores-Filter.pdf"), width = 10, height = 12, bg = 'transparent')
  grid.arrange(PlotDeviations(dt = dtVierstraVarDevFilt, n = ntop, text.adjust = T),
               PlotRankCorrelations(dt = dtMotifGeneCorTopCorr, n = ntop, text.adjust = T),
               nrow=1)
  dev.off()

  ## ChromVar variations
  dtMotifGeneCorTopCorSample <- dtMotifGeneCorSample[order(-rho, pval, decreasing = F)][, head(.SD, 1), by=list(gene,Sample)] # top by corr
  dtMotifGeneCorTopCorSample <- dtMotifGeneCorTopCorSample[,head(.SD, 1), by=list(family, Sample)]
  dtMotifGeneCorTopCorSample[,corRank:=seq_along(1:.N), by='Sample']
  ntop <- 30

  pl <- lapply(sort(unique(dtMotifGeneCorTopCorSample$Sample)), function(i){
    PlotRankCorrelations(dt = dtMotifGeneCorTopCorSample[Sample==i],
                         n = ntop,
                         text.adjust = T,
                         title = paste(i,"top",ntop,"correlated deviation/gene"))
    })

  cairo_pdf(filename = file.path(plot.dir, "Variable-Vierstra-Motif-Deviation-Scores-Filter-bySample.pdf"), width = 30, height = 12, bg = 'transparent')
  do.call(grid.arrange, c(pl, nrow=1))
  dev.off()

  #### Plot UMAP AR modalities ####

  writexl::write_xlsx(dtTopVarDev,
                      path = file.path(spreadsheets.dir,'multiome_chromvar_deviations_Vierstra.xlsx'),
                      col_names = T, format_headers = T)

  saveArchRProject(ArchRProj = proj)

  #--------------------------------------- Plots TF Locus in Browser ---------------------------------------
  #### Browse peaks ####
  ArchRBrowser(proj)

  p <- plotBrowserTrack(
    ArchRProj = proj,
    #region = GRanges(markerList$`C1-Luminal Mature`[3,]),
    groupBy = "Clusters",
    #geneSymbol = c("Foxa1","Trp63","Cd55","Ctcf","Krt19","Krt14","Krt18","Pgr","Gata3","Esr1"),
    geneSymbol = sort(unique(unname(unlist(markers_list)))),
    title = '(+- 20Kb)',
    plotSummary = c("bulkTrack","scTrack","featureTrack","geneTrack"),
    sizes = c(10, 10, 2, 2.5),
    # features =  getMarkers(clusterMarkersPeaks, cutOff = cutoff, returnGR = TRUE), # We can restrict to significant peaks or use all peaks
    upstream = 20000,
    downstream = 20000
  )
  grid.draw(p$Foxa1)
  grid.draw(p$Krt5)
  grid.draw(p$Krt15)

  plotPDF(p, name = "Plot-Tracks-With-Features_Groups", width = 6, height = 12, ArchRProj = proj, addDOC = FALSE)

  p <- plotBrowserTrack(
    ArchRProj = proj,
    #region = GRanges(markerList$`C1-Luminal Mature`[3,]),
    groupBy = "Clusters",
    geneSymbol = c("Hif3a","Pou2f1","Pou2f2","Pou2f3","Ar","Lef1",""),
    title = '(+- 20Kb)',
    plotSummary = c("bulkTrack","scTrack","featureTrack","geneTrack"),
    sizes = c(10, 10, 2, 2.5),
    # features =  getMarkers(clusterMarkersPeaks, cutOff = cutoff, returnGR = TRUE), # We can restrict to significant peaks or use all peaks
    upstream = 20000,
    downstream = 20000
  )
  grid.draw(p$Hif3a)

  plotPDF(p, name = "Plot-Tracks-With-Features_Groups", width = 6, height = 12, ArchRProj = proj, addDOC = FALSE)

  saveArchRProject(ArchRProj = proj, load = TRUE)

  # bionames de tests
  nset <- c("Basal","Luminal") #sort(unique(proj@cellColData$luminal_basal))
  (combinations <- combn(x = nset, m = 2))
  seDELBTests <- lapply(seq_len(ncol(combinations)), function(j){
    pair <- combinations[,j]
    getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = "luminal_basal",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = pair[1],
      bgdGroups = pair[2]
    )
  })
  (names(seDELBTests) <- apply(combinations, 2, function(x) { paste(x,collapse = "_vs_")}))

  # Takes a long time
  # message("Plotting pariwise with cutoff: ",cutoff)
  # seDEPairPlots <- lapply(names(seDETests), function(i){
  #   x <- seDETests[[i]]
  #   if(is.null(x))
  #     return(x)
  #   fg <- x@metadata$Params$useGroups
  #   bg <- x@metadata$Params$bgdGroups
  #   g <- plotMarkers(seMarker = x, name = fg, cutOff = cutoff, plotAs = "MA")
  #   g <- g + labs(subtitle = i) + theme(plot.margin=unit(c(0,0,0,0), "cm"))
  # })
  # plotPDF(seDEPairPlots, name = "MA_Peaks_in_pairwise_de_clusters.pdf", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

  PlotBarDEPeaks(ArchRProj = proj, absLog2FC = 0.5, fdr = 0.1, seDETests = seDEbioNamesTests, width = 8, height = 6, outfile.name = "Peaks_in_pairwise_de_bioNames_Barplot.pdf")

  PlotBarDEPeaks(ArchRProj = proj, seDETests = seDEClusterTests, outfile.name = "Peaks_in_pairwise_de_clusters_Barplot.pdf")

  PlotBarDEPeaks(ArchRProj = proj, absLog2FC = 1, fdr = 0.01, seDETests = seDELBTests, width = 3, height = 6, outfile.name = "Peaks_in_pairwise_de_luminal_basal_Barplot.pdf")

  nset <- sort(unique(proj@cellColData$Clusters2))
  (combinations <- combn(x = nset, m = 2))
  deTests <- apply(combinations, 2, function(x){ paste0(x[1],'_vs_',x[2]) } )

  library(doParallel)
  registerDoParallel(cores = 4)

  # Finds motifs enriched in first vs second
  motifsUp <- foreach(i = names(seDETests)) %dopar% {
    peakAnnoEnrichment(
      seMarker = seDETests[[i]],
      ArchRProj = proj,
      peakAnnotation = "Motif",
      cutOff = peak.cutoff)
  }
  names(motifsUp) <- names(seDETests)

  # Finds motifs enriched in second vs first
  motifsDown <- foreach(i = names(seDETests)) %dopar% {
    peakAnnoEnrichment(
      seMarker = seDETests[[i]],
      ArchRProj = proj,
      peakAnnotation = "Motif",
      cutOff = peak.cutoff)
  }
  names(motifsDown) <- names(seDETests)

  dtup <- lapply(names(motifsUp), function(i){
    se <- motifsUp[[i]]
    dt <- data.table(do.call(cbind, se@assays@data)) %>% setnames(names(se@assays@data))
    #dt <- data.table(TF = rownames(se), mlog10Padj = assay(se)[,1])
    dt <- dt[order(mlog10Padj, decreasing = TRUE),]
    dt$rank <- seq_len(nrow(dt))
    dt$group <- i
    dt
  }) %>% rbindlist

  dtdown <- lapply(names(motifsDown), function(i){
    se <- motifsDown[[i]]
    dt <- data.table(do.call(cbind, se@assays@data)) %>% setnames(names(se@assays@data))
    #dt <- data.table(TF = rownames(se), mlog10Padj = assay(se)[,1])
    dt <- dt[order(mlog10Padj, decreasing = TRUE),]
    dt[,mlog10Padj := -1*mlog10Padj]
    dt$rank <- seq_len(nrow(dt))
    dt$group <- i
    dt
  }) %>% rbindlist

  ## Plot motif grid
  pdf(file = file.path(getOutputDirectory(proj),'Plots',"Motif-Enrichment-in-DE-peaks.pdf"), width = 10, height = 5, bg = 'transparent', onefile = T)
  for(i in names(seDETests)){
    up <- motifsUp[[i]]
    down <- motifsDown[[i]]
    pUp <- PlotMotifs(up, topn = 10) + labs(title = paste(gsub("_vs_"," vs\n",i),"Up")) + theme(plot.title = element_text(size = 12))
    pDown <- PlotMotifs(down, topn = 10) + labs(title = paste(gsub("_vs_"," vs\n",i),"Down")) + theme(plot.title = element_text(size = 12))
    grid.arrange(pUp, pDown, nrow=1)
  }
  dev.off()

  (g1 <- PlotMotifsTileGridSymmetric(seListUp = motifsUp, seListDown = motifsDown, "Motifs", scale = 'group', topn=10))
  ggsave(plot = g1, filename = file.path(getOutputDirectory(proj),'Plots',"Motif-Enrichment-in-DE-peaks-grid_groupscale.pdf"), width=22, height=12, units = 'in',
         device = cairo_pdf, dpi=320, bg = 'transparent')

  (g1 <- PlotMotifsTileGridSymmetric(seListUp = motifsUp, seListDown = motifsDown, "Motifs", scale = 'motif', topn=10))
  ggsave(plot = g1, filename = file.path(getOutputDirectory(proj),'Plots',"Motif-Enrichment-in-DE-peaks-grid_motifscale.pdf"), width=22, height=12, units = 'in',
         device = cairo_pdf, dpi=320, bg = 'transparent')

  saveArchRProject(ArchRProj = proj)

  #### Plot scRNA Integration UMAP ####
  markerRNA <- getFeatures(proj, select = paste( gsub("_.*","",c(motifs,'Tfap2')), collapse="|"), useMatrix = "GeneIntegrationMatrix")
  markerRNA <- markerRNA[markerRNA %ni% c("Cd55b","Fosb","Fosl1","Fosl2","Wt1os","Ctcflos","Ctcfl","Celf5","Esr1",
                                          "Ppifos","Pgr15l","Pgrmc1","Pgrmc2","Rpgr","Rpgrip1","Rpgrip1l",
                                          "Tfap2b","Tfap2c","Tfap2e","Swt1")]
  sort(markerRNA)
  p <- plotEmbedding(
    ArchRProj = proj2,
    colorBy = "GeneIntegrationMatrix",
    name = sort(markerRNA),
    embedding = "UMAP",
    plotAs = 'points',
    size = 0.5,
    imputeWeights = getImputeWeights(proj2)
  )
  p1 <- lapply(p, function(x){
    UmapPlotEnhance(p1 = x, addLabels = F, title = x$labels$title, legend.title = 'scRNA linked expression')
  })

  png(filename = file.path(getOutputDirectory(proj2),'Plots','UMAP-GeneScores-scRNA.png'), width = 20, height = 8,
      units = 'in', bg = 'transparent', res = 150)
  do.call(grid.arrange, c(p1, nrow=3, top="Inferred gene expression per cell (scRNA)"))
  dev.off()

  #### Deviation z scores plot ####
  motifmat <- getMatrixFromProject(ArchRProj = proj2, useMatrix = "MotifMatrix")
  z <- assays(motifmat)

  png(filename = file.path(getOutputDirectory(proj2),'Plots','EnrichedMotifs-devZscore-density.png'),
      width = 12, height = 5, units = 'in', bg = 'transparent', res = 150)
  par(mfrow=c(2,7))
  for(i in motifs){
    zn <- grep(i, rownames(z$z), value = T)[1]
    zd <- density(z$z[zn,])
    hist(z$z[zn,], freq=FALSE, breaks = seq(from=min(z$z[zn,]), to=max(z$z[zn,]), length.out = 50),
         col = 'lightblue', angle=75, border='white', main=zn, ylab="Density", xlab="dev. zscore")
    polygon(zd, col=alpha('steelblue', 0.1), border='steelblue', lwd=2)
  }
  dev.off()

  #### Deviation z scores per cluster ####
  p1 <- PlotDeviationZScoreDensity(motifmat, motifs)

  png(filename = file.path(getOutputDirectory(proj2),'Plots','EnrichedMotifs-devZscore-percluster-density.png'),
      width = 20, height = 10, units = 'in', bg = 'transparent', res = 150)
  l <- lapply(p1, function(x){x + theme(legend.position = 'none')})
  l$legend <- ggpubr::get_legend(p1[[1]])
  do.call(grid.arrange, c(l, nrow=3))
  dev.off()

  ## Scores distribution
  p <- plotGroups(ArchRProj = proj2,
                  groupBy = "Clusters2",
                  colorBy = "MotifMatrix",
                  name = markerMotifs,
                  imputeWeights = NULL #getImputeWeights(proj2)
  )
  png(filename = file.path(getOutputDirectory(proj2),'Plots','Plot-Groups-Deviations-w-Imputation.png'),
      width = 20, height = 7, units = 'in', bg = 'transparent', res = 150)
  l <- lapply(p, function(x){x + ylab("") + theme_classic()+ theme(plot.background = element_blank(),
                                                                   plot.margin = unit(c(0,0.1,0,0), units = 'cm'),
                                                                   text = element_text(size=10))})
  do.call(grid.arrange, c(l, nrow=2))
  dev.off()

  #### accessible regions per cluster ####
  # Identify the differential peaks per cluster and then you can subset those peaks from the motifmatches.
  # Try getMatches() to get a peak x motif occurrence matrix and then subset it with GRanges manipulations or by pasting the ranges i.e paste0(seqnames(gr),"",start(gr),"",end(gr)) and doing string intersection.

  #### Motif footprinting ####

  # We need large number of reads
  # cells are grouped to create pseudo-bulk profiles and stored as group coverage files
  proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters2")

  motifPositions <- getPositions(proj2) # gets motif positions from first motif matrix = CIS-BP

  motifPositions

  ##motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
  markerMotifs <- sort(unlist(lapply(c(motifs,'Klf5','Smad1'), function(x) grep(x, names(motifPositions), value = TRUE))))
  #markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
  markerMotifs

  seFoot <- getFootprints(
    ArchRProj = proj2,
    positions = motifPositions[markerMotifs],
    groupBy = "Clusters2"
  )

  pf <- plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj2,
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias",
    addDOC = FALSE,
    smoothWindow = 5,
    plot = T,
    height=7,
    width=5,
    baseSize = 5
  )

  # Number of motifs found
  sapply(motifPositions, length)

  #### Co-accessibility ####
  proj2 <- addCoAccessibility(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI"
  )

  # get DF of results
  cA <- getCoAccessibility(
    ArchRProj = proj2,
    corCutOff = 0.5,
    resolution = 1000,
    returnLoops = T
  )

  p <- plotBrowserTrack(
    ArchRProj = proj2,
    groupBy = "Clusters2",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000,
    loops = cA
  )

  grid.draw(p$Esr1)
  grid.draw(p$Krt8)
  grid.draw(p$Sparc)
  plotPDF(plotList = p,
          name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
          ArchRProj = proj2,
          addDOC = FALSE, width = 5, height = 5)

  gr <- head(cA$CoAccessibility[order(cA$CoAccessibility$value, decreasing = T),], n=20)

  pGR <- plotBrowserTrack(
    ArchRProj = proj2,
    groupBy = "Clusters2",
    region = gr,
    #geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000,
    loops = cA
  )

  grid.draw(pGR[[20]])

  #### Peak2Gene links ####
  ## Adopts approach of Cicero to create low-overlapping aggregates of single cells.
  ##
  proj2 <- addPeak2GeneLinks(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI"
  )

  p2g <- getPeak2GeneLinks(
    ArchRProj = proj2,
    corCutOff = 0.45,
    DRCutOff = 1e-04,
    resolution = 1,
    returnLoops = FALSE
  )

  p2g
  metadata(p2g)[[1]]

  p <- plotBrowserTrack(
    ArchRProj = proj2,
    groupBy = "Clusters2",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(projHeme5)
  )

  ## Heatmap
  p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "Clusters2", k = 7,
                            palGroup = rna_group_palette,
                            palATAC =  ArchRPalettes$solarExtra,
                            palRNA = ArchRPalettes$blueYellow)

  for(i in 1:length(p@ht_list)){
    p@ht_list[[i]]@top_annotation@anno_list[[1]]@label <- ""
    p@ht_list[[i]]@top_annotation@anno_list[[1]]@color_mapping@name <- "clusters"
    p@ht_list[[i]]@top_annotation@anno_list[[1]]@name_param$label <- ""
  }

  png(filename = file.path(getOutputDirectory(proj2),'Plots','ATAC-RNA-links-heatmap.png'),width = 10, height = 12, units='in', res = 150, bg = 'transparent')
  draw(p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()

  pdf(file = file.path(getOutputDirectory(proj2),'Plots','ATAC-RNA-links-heatmap.pdf'),width = 10, height = 12, bg = 'transparent')
  draw(p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()

}

#' Create the ArchR project from Cell Ranger ARC output
#'
#' Creates Arrow files for every multiome sample found in `raw.dir`, computes
#' doublet scores, builds the ArchRProject, filters cells by TSS enrichment,
#' doublets and the cells retained in the RNA analysis, then runs iterative LSI,
#' clustering and UMAP and saves the QC / cluster plots. Sample MMO2 is excluded.
#' Run section by section.
BuildProj <- function(){

  multiome.data <- list.files(path = raw.dir, 'multiome', full.names = T)
  multiome.sample <- names(multiome.data) <- basename(multiome.data)
  message("Removing MMO2")
  multiome.sample <- grep("^MMO2", multiome.sample, invert = T, value = T)
  multiome.data <- multiome.data[multiome.sample]

  (rna.counts <- list.files(path = multiome.data, pattern = '^filtered_feature_bc_matrix.h5$', recursive = T, full.names = T))

  (atac.frags <- list.files(path = multiome.data, pattern = '^atac_fragments.tsv.gz$', recursive = T, full.names = T))

  sample.list <- lapply(multiome.sample, function(x){
    list('rna'=grep(x, rna.counts, value = T),
         'atac'=grep(x, atac.frags, value = T))
  }) %>% set_names(multiome.sample)

  (fragment.files <- sapply(names(sample.list), function(x) {  sample.list[[x]]$atac }))

  #### Creating Arrow Files ####
  minTSS <- 1

  ArrowFiles <- createArrowFiles(
    inputFiles = fragment.files,
    sampleNames = names(fragment.files),
    minTSS = minTSS, #Minimum TSS score required Dont set this too high because you can always increase later
    minFrags = nfragMin, #1000, #Minimum number of mapped ATAC-seq fragments required
    maxFrags = 1e+06,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    force = TRUE,
    threads = 4,
    verbose = TRUE
  )
  ArrowFiles

  #### Inferring Doublets ####
  doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
    LSIMethod = 1,
    threads = 8
  )
  names(doubScores) <- ArrowFiles

  #### Creating an ArchRProject ####
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = archr.dir,
    copyArrows = F #T
  )

  paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), " MB")

  getAvailableMatrices(proj)

  quantile(proj$TSSEnrichment)

  table(proj$Sample)
  proj$bioNames <- unlist(sapply(proj$Sample, function(x) { INPUT_LABELS[[x]] }))
  proj$bioNames2 <- unlist(sapply(proj$Sample, function(x) { INPUT_LABELS[[x]] }))
  #proj@cellColData$bioNames <- factor(x = proj@cellColData$bioNames, levels = names(bionames_palette_list), ordered = T)
  proj$Condition <- unlist(sapply(proj$Sample, function(x) { CONDITION_LABELS[[x]] }))
  proj$Treatment <- unlist(sapply(proj$Sample, function(x) { TREATMENT_LABELS[[x]] }))

  getCellColData(proj, select = c("bioNames"))
  getCellColData(proj, select = c("Condition"))
  getCellColData(proj, select = c("Treatment"))

  ## Can subset here, if desired
  library(SingleCellExperiment)

  # pbmc[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA.
  # To assist with visualization and interpretation. we also convert Pearson residuals back to ‘corrected’ UMI counts.
  # You can interpret these as the UMI counts we would expect to observe if all cells were sequenced to the same depth. If you want to see exactly how we do this, please look at the correct function here.
  # The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@counts.
  # We store log-normalized versions of these corrected counts in pbmc[["SCT"]]@data, which are very helpful for visualization.

  scRNA.multiome <- readRDS(file = multiRNA.file)
  a <- SingleCellExperiment(assays = list(counts = GetAssayData(object = scRNA.multiome, slot = "counts"),
                                          data = GetAssayData(object = scRNA.multiome, slot = "data"),
                                          normcounts = GetAssayData(object = scRNA.multiome, slot = "scale.data"),
                                          rnarawcounts = GetAssayData(object = scRNA.multiome, assay='RNA', slot = "counts"),
                                          rnadata = GetAssayData(object = scRNA.multiome, assay='RNA', slot = "data")),
                            colData=DataFrame(scRNA.multiome@meta.data),
                            rowData=DataFrame(scRNA.multiome@assays$SCT@meta.features),
                            reducedDims = scRNA.multiome@reductions)
  head(colnames(a))
  head(a@colData)
  message("# cells in common:", length(intersect(colnames(a), rownames(proj@cellColData))))
  message("# cells not in ATAC:", length(setdiff(colnames(a), rownames(proj@cellColData))))
  message("# cells not in RNA:", length(setdiff(rownames(proj@cellColData), colnames(a))))
  adiff <- lapply(setdiff(colnames(a), rownames(proj@cellColData)), function(x){
    as.data.table(a@colData[x,], keep.rownames = 'cells')
  }) %>% rbindlist()
  table(adiff$seurat_clusters)

  ## Plot and place in separate directory
  QCPlots(proj)

  # Already done during ArrowFile creation
  # idxPass <- which(proj$TSSEnrichment >= 8)
  # idxPass <- which(proj$nFrags >= nfragMin)
  # cellsPass <- proj$cellNames[idxPass]
  # proj <- proj[cellsPass, ]
  ## Inspected output and decided to filter for larger TSS enrichment
  proj

  range(proj$TSSEnrichment)
  range(proj$nFrags)

  proj <- proj[proj$nFrags >= 1000]
  proj <- proj[proj$TSSEnrichment >= 4]

  #### Save object ####
  proj <- saveArchRProject(ArchRProj = proj, load = T)

  #### Load object ####
  proj <- loadArchRProject(path = archr.dir)

  QCPlots(proj)

  #### Before doublet removal ####
  #### LSI Dimensionality Reduction and Clustering with Latent Semantic Indexing ####
  # https://www.archrproject.com/bookdown/iterative-latent-semantic-indexing-lsi.html
  # Creates a reducedDims object.
  # components correlated with depth are filtered out based on corCutoff.
  ## Identify doublets
  dfDoublets <- filterDoubletsFix(proj, filterRatio = .5, dryrun = T)
  write.table(file = file.path(data.dir, "scATAC_doublets.txt"), x = dfDoublets, quote = F, sep = "\t", col.names = T)

  ## Remove doublets ##
  proj <- filterDoubletsFix(proj, filterRatio = .5, dryrun = F)

  proj <- saveArchRProject(ArchRProj = proj, outputDirectory = archr.dir, load = T)

  ### Keep seurat cells ###
  arid1a.multiome <- readRDS(multiRNA.file)

  table(gsub("#.*","",rownames(proj@cellColData)))
  table(gsub("#.*","",Cells(arid1a.multiome)))
  #### Remove cells not in ATAC ####
  message("Filtering cells removed in RNA")
  proj <- subsetArchRProject(proj, intersect(ArchR::getCellNames(proj),
                                             Cells(arid1a.multiome)),
                             outputDirectory = archr.dir,
                             force = TRUE)
  table(gsub("#.*","",rownames(proj@cellColData)))
  table(gsub("#.*","",Cells(arid1a.multiome)))

  setequal(rownames(proj@cellColData), Cells(arid1a.multiome))

  proj <- saveArchRProject(ArchRProj = proj, outputDirectory = archr.dir, load = T)

  message("Consider the use of scaleDims")
  #### We will not remove doublets as ArchR warned about little heterogeneity, dimension reduction (LSI) and evaluate clusters ####
  #### Make final determinaton of clustering ####
  #### Should by now have a good idea which resolution and clustering to apply
  proj
  #### LSI Dimensionality Reduction and Clustering with Latent Semantic Indexing ####
  # https://www.archrproject.com/bookdown/iterative-latent-semantic-indexing-lsi.html
  # Creates a reducedDims object.
  # components correlated with depth are filtered out based on corCutoff.
  # choose preferred resolution
  proj <- addIterativeLSI(ArchRProj = proj,
                           useMatrix = "TileMatrix",
                           force = T,
                           corCutOff = 0.75, # removes dimensions having more that this correlation with sequencing depth
                           iterations = 3,
                           clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start=10), # See Seurate:FindClusters, uses KNN->SNN
                           varFeatures = 25000,
                           scaleDims = F,  # Set to False so that we remove highly correlated dimensions
                           dimsToUse = 1:30)

  message("Adding umaps for bioNames")
  proj <- addUMAP(ArchRProj = proj,
                  force = T,
                  reducedDims = "IterativeLSI",
                  name = "UMAP",
                  nNeighbors = 30,
                  seed = 42,
                  minDist = 0.5,
                  spread = 1,
                  metric = "cosine")

  ## test
  proj1 <- addUMAP(ArchRProj = proj,
                  force = T,
                  reducedDims = "IterativeLSI",
                  name = "UMAP",
                  nNeighbors = 30,
                  seed = 42,
                  minDist = .5,
                  spread = .8,
                  metric = "cosine")
  p <- plotEmbedding(ArchRProj = proj1,
                     colorBy = "cellColData",
                     name = "gex_mammary_gland_refine_clusters",
                     embedding = "UMAP")
  p <- UmapPlotEnhance(ArchRProj = proj1,
                       p1 = p,
                       addLabels = F,
                       title.size = 16,
                       pal =  mammary_refine_clusters_palette,
                       legend.title = paste("cell type"),
                       legend.title.rotation = 0)
  p

  #--
  pTruePromoter <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                                 name = "ReadsInPromoter", embedding = "UMAP",
                                 size=0.4, baseSize = 14, plotAs = 'points') %>%
    UmapPlotEnhance(ArchRProj = proj, p1 = ., addLabels = F,
                    pt.size = 2,
                    title.size = 16,
                    break.unit = 10000,
                    legend.title = "Reads in Promoter")
  pTruePromoter

  p <- plotEmbedding(ArchRProj = proj,
                       colorBy = "cellColData",
                       name = "bioNames",
                       embedding = "UMAP")
  p <- UmapPlotEnhance(ArchRProj = proj,
                       p1 = p,
                       addLabels = F,
                       title.size = 16,
                       pal =  bionames_palette_list,
                       legend.title = paste("Samples"),
                       legend.title.rotation = 0,
                       filename.prefix = "UMAP_bionames")
  p

  scaleDimsFalseClustering <- EvaluateUMAPClustering(ArchRProj = proj,
                                                     ncores = 6,
                                                     file.prefix = 'scaleDimsFalse',
                                                     scaleClusterDims = F,
                                                     resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5,
                                                                     0.6, 0.7, 0.8, 0.9, 1.0,
                                                                     1.1, 1.2, 1.3, 1.4, 1.5,
                                                                     1.6, 1.7, 1.8, 1.9, 2))

  #--

  #projscaleDimsFalseClustering <- EvaluateUMAPClustering(ArchRProj = proj, file.prefix = 'scaleDimsFalse', minDist = 0.8, spread = 1.2)

  ## Plot bionames to see samples
  p1 <- plotEmbedding(ArchRProj = projscaleDimsFalseClustering$rlist[[1]], colorBy = "cellColData", name = "bioNames", embedding = "UMAP")
  p1 <- UmapPlotEnhance(ArchRProj = projscaleDimsFalseClustering$rlist[[1]],
                        p1 = p1,
                        addLabels = F,
                        outline.color = 'black',
                        pt.size = .8,
                        title.size = 16,
                        pal = bionames_palette_list)
  p1

  library(MASS)
  library(plotly)
  df <- projscaleDimsFalseClustering$rlist[[1]]@embeddings$UMAP$df
  df$mammary_gland_cluster_classification <- proj1$mammary_gland_cluster_classification
  colnames(df) <- c('x','y','z','mammary_gland_cluster_classification')

  p <- plot_ly(df, x=df$x, y=df$y, z=df$z,
              color=~mammary_gland_cluster_classification) %>% add_markers(size=.5)
  p

  ## calculate silhouette scores and plot clusters
  resm <- names(projscaleDimsFalseClustering$rlist)
  sils <- lapply(resm, function(i){
    x <- as.integer(gsub("^C","",projscaleDimsFalseClustering$rlist[[i]]$Clusters))
    d <- dist(projscaleDimsFalseClustering$rlist[[i]]@embeddings$UMAP$df)
    sil <- silhouette(x = as.numeric(x = as.factor(x = x)), dist=d)
    #(sshi <- summary(sil))
  }) %>% set_names(resm)

  ## Umaps for different cluster resolutions
  um <- lapply(resm, function(r){
    score <- summary(sils[[r]])$avg.width
    p1 <- projscaleDimsFalseClustering$rplots[[r]]
    p1 <- p1 + ggtitle(paste("Cluster res.:", gsub("SCT_snn_res.","",r),
                             "\nmean silhouette ", formatC(score, digits = 3)))
  }) %>% set_names(resm)

  png(filename = file.path(plot.dir,paste0('scaleDimsFalse_clusters_assess_spread_1.png')),
      width = 22, height = 8, units = 'in', res = 200)
  do.call(grid.arrange, c(um, nrow=2))
  dev.off()

  library(cluster)
  pdf(file = file.path(plot.dir,'silhouette_scaleDimsFalseUMAP_final.pdf'), width = 10, height = 8, onefile = T)
  for(r in resm){
    sil <- sils[[r]]
    p1 <- ~plot(sil, main=paste("resolution:",r), cex.names = par("cex.axis")) # RStudio sometimes does not display silhouette plots correctly
    grid.arrange(cowplot::as_grob(p1), projscaleDimsFalseClustering$rplots[[r]], nrow=1, widths=c(.4,.6))
  }
  dev.off()

  r <- projscaleDimsFalseClustering$rlist[['0.1']]

  #### Plot final UMAPS ####
  table(proj@cellColData$bioNames)
  #BioNames
  pclust <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                          name = "Clusters", embedding = "UMAP") %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    addLabels = T,
                    title.size = 16,
                    pal = NA,
                    label.size = 4)

  pBio <- plotEmbedding(ArchRProj = proj,
                        colorBy = "cellColData",
                        name = "bioNames",
                        embedding = "UMAP",
                        size=0.2, baseSize = 14) %>%
    UmapPlotEnhance(ArchRProj = proj, p1 = ., addLabels = F,
                    label.size = 4,
                    pal=bionames_palette_list,
                    title.size = 18,
                    panel.border.size = 0,
                    filename.prefix = "UMAP-ATAC-BioNames")

  png(filename = file.path(plot.dir,'UMAP-ATAC-BioNamesClusters.png'), width = 5.5, height = 3.5, units = 'in', bg = 'transparent', res = 250)
  pBio+pclust$layers
  dev.off()

  #Cluster #
  plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                name = "Clusters", embedding = "UMAP") %>%
    UmapPlotEnhance(ArchRProj = proj,
                    p1 = .,
                    addLabels = T,
                    title.size = 16,
                    pal = cluster_palette_list,
                    label.size = 6,
                    panel.border.size = 0,
                    filename.prefix = "UMAP-ATAC-ClusterNumber")

  #nFrag
  dtl <- list(data=data.table(proj@embeddings$UMAP$df, log10(proj@cellColData$nFrags), keep.rownames = 'cells') %>% setnames(.,c('cells',"x","y","color")))
  dtl$labels$x <- 'UMAP 1'
  dtl$labels$y <- 'UMAP 2'
  p <- UmapPlotEnhance(ArchRProj = proj,
                       break.unit = 1,
                       legend.height.factor = .3,
                       p1 = dtl,
                       addLabels = F,
                       panel.border.size = 0,
                       title.size = 16,
                       scale_color_bar = F,
                       legend.title = "log10(nFrags)",
                       filename.prefix = "UMAP-ATAC-nFRAG")

  #### Plot bar abundance ####
  PlotBarVariableInClusters(proj, variable = 'bioNames', height_scale_perc = .8,
                            palette = bionames_palette_list,
                            orientation = 'vertical',
                            filename = "bionames_bar")
  PlotBarVariableInClusters(proj, variable = 'Condition', height_scale_perc = .8,
                            orientation = 'vertical', filename = "condition_bar")
  PlotBarVariableInClusters(proj, variable = 'gex_mammary_gland_cluster_classification', height_scale_perc = .8,
                            orientation = 'vertical', palette = mammary_cluster_class_palette, filename = "mammary_bar")

  pdf(file = file.path(plot.dir,'DepthDimensionCor.pdf'), width = 8, height = 4)
  PlotCorDimensionDepth(proj)
  dev.off()

  #### Confusion Matrix Sample and Cluster ####
  table(proj$Clusters)
  cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$bioNames))
  cM
  library(pheatmap)
  # apply(cM, 1, function(x){ x/sum(x)})
  cM <- t(cM) / rowSums(t(cM))
  cM
  p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = ArchRPalettes$solarExtra, #
    clustering_method = 'ward.D',
    treeheight_row = 20,
    treeheight_col = 20,
    cellwidth = 20,
    cellheight =20,
    #color = RColorBrewer::brewer.pal(n=8, name = "YlOrRd"), # paletteContinuous(set = "Set1"),
    border_color = "black",
    number_color = 'white',
    display_numbers = T,
    number_format="%.2f"
  )

  plotPDF(p, name = "confusion-matrix-wardmethod-cluster-sample.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 3)

  #### Save object ####
  proj <- saveArchRProject(ArchRProj = proj, load = T)

}

#### Helpers ####
#' Value with the largest absolute magnitude in a vector
absmax <- function(x) { x[which(abs(x)==max(abs(x)))][1] }

#' Convert an ArchR matrix (default: motif deviations) to a Seurat/SCE object
#'
#' @param build If FALSE, read the previously saved object from `output.file`.
#' @param ArchRProj ArchR project.
#' @param matrix Name of the ArchR matrix to export.
#' @param output.file .rds file to read or write.
ArchRProjectToSeurat <- function(build=T,
                                 ArchRProj,
                                 matrix = 'VierstraMatrix',
                                 output.file=file.path(robj.dir,'arid1a_archR_atac_peaks_multiome_sce.rds')){

  if(!build){
    return(readRDS(output.file))
  }

  library(SingleCellExperiment)
  library(Seurat)

  message("Using ",matrix)

  rse <- getMatrixFromProject(ArchRProj, matrix, threads = 8)

  if(matrix=='GeneScoreMatrix'){
    names(assays(rse)) <- "counts"
    assays(rse) <- list(counts = assays(rse)$counts,
                        logcounts = log1p(assays(rse)$counts) )
    rownames(rse) <- rowData(rse)$name
  }
  else if(matrix=='PeakMatrix'){
    names(assays(rse)) <- "counts"
    assays(rse) <- list(counts = assays(rse)$counts,
                        logcounts = log1p(assays(rse)$counts) )
    rownames(rse) <- paste0("Peak",seq_along(rse))
  }
  else if(matrix=='VierstraMatrix'){
    i <- names(assays(rse))
    i <- sapply(i, function(j){
      ifelse(j=="deviations","logcounts",ifelse(j=="z","counts",j)) })
    names(assays(rse)) <- i
    assays(rse) <- list(counts = assays(rse)$counts,
                        logcounts = assays(rse)$logcounts )
  }

  ## create single cell experiment
  sce <- as(rse, "SingleCellExperiment")

  ## create seurat object
  message("Creating Seurat object")
  obj <- CreateSeuratObject(counts = counts(sce), data = logcounts(sce), meta.data = as.data.frame(colData(sce)))
  obj@assays
  obj <- RenameAssays(obj, RNA = 'ATAC')

  lsi <- ArchRProj@reducedDims$IterativeLSI$matSVD[,-1]
  umap <- as.matrix(ArchRProj@embeddings$UMAP$df)

  obj[['LSI']] <- Seurat::CreateDimReducObject(embeddings = lsi,
                                               key = "LSI_",
                                               assay = "ATAC")
  obj[['UMAP']] <- Seurat::CreateDimReducObject(embeddings = umap,
                                                key = "UMAP_",
                                                assay = "ATAC")

  if(length(output.file) > 0)
    if(!is.na(output.file)){
      message("Saving output to ",output.file)
      saveRDS(object = f, file=output.file)
    }

  f

}

# spread can therefore be used to control the inter-cluster distances to some extent,
# where as min_dist controls the size of the clusters.
# https://jlmelville.github.io/uwot/abparams.html
#' Evaluate clustering resolutions on an ArchR project
#'
#' Runs clustering at several resolutions on the IterativeLSI reduction,
#' computes UMAPs and silhouette scores for each, and writes diagnostic plots.
#'
#' @param ArchRProj ArchR project.
#' @param file.prefix Prefix of the output files.
#' @param minDist,spread UMAP parameters.
#' @param resolutions Clustering resolutions to evaluate.
EvaluateUMAPClustering <- function(ArchRProj, file.prefix="", minDist=0.5, spread=1,
                                   scaleClusterDims=F, ncores=8,
                                   resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5,
                                                   0.6, 0.7, 0.8, 0.9, 1.0,
                                                   1.1, 1.2, 1.3, 1.4, 1.5)) {
  library(cluster)

  dt <- data.table(ArchRProj@reducedDims$IterativeLSI$matSVD[,c(1,2,3)])
  dt[,ReadsInPromoter := ArchRProj@cellColData$ReadsInPromoter]
  dt[,labels := as.character(ArchRProj@cellColData$bioNames)]
  colors <-  sapply(dt$labels, function(x) { bionames_palette_list[[x]] })

  png(filename = file.path(plot.dir,paste0(file.prefix,'_LSI_Correlation.png')),
      width = 6, height = 6, units = 'in', res = 200)
  plot(dt[,list(LSI1, LSI2, LSI3, ReadsInPromoter)], col=colors, pch=1, lwd=1, cex=.5)
  dev.off()

  png(filename = file.path(plot.dir,paste0(file.prefix,'_clusters_pca.png')),
      width = 6, height = 6, units = 'in', res = 200)
  plot(dt$LSI2, dt$LSI3, col=colors, pch=1, lwd=1, cex=.5)
  # Add a legend
  legend("bottomleft", inset=.02, title="Sample", bg = NULL,
         names(bionames_palette_list), col=unlist(bionames_palette_list), horiz=FALSE, cex=.8, pch=1)
  dev.off()

  #### Clustering ####
  # quick clustering approach
  # Uses Seurat::FindClusters(), RAM hungry because "Seurat" casts data as dense 'matrix' rather than keeping it sparse
  # must use my own method for clustering, it uses Seurat sparse matrix, but it is slow
  rclust <- mclapply(resolutions, function(r){
    message("computing for resolution: ",r)
    #addClustersIgraph(
    addClusters(input = ArchRProj,
                reducedDims = "IterativeLSI",
                method = "Seurat",
                name = "Clusters",
                algorithm = 4, #leiden not louvain
                resolution = r,
                ncores = ncores,
                scaleDims = scaleClusterDims, # Uses value from LSI ArchRProj@reducedDims$IterativeLSI$scaleDims
                force = T)
  }, mc.cores = ncores ) %>% set_names(resolutions)

  if(!dir.exists(robj.dir))
    dir.create(robj.dir)
  saveRDS(object = rclust, file = file.path(robj.dir,'archrR_clustering_resolutions.RDS'))

  message("Adding umaps for resolutions")
  rlist <- mclapply(rclust, function(r){
    addUMAP(ArchRProj = r,
            force = T,
            reducedDims = "IterativeLSI",
            name = "UMAP",
            nNeighbors = 30,
            seed = 42,
            minDist = minDist,
            spread = spread,
            metric = "cosine")
  }, mc.cores = ncores) %>% set_names( names(rclust) ) # efficiently operates on sparse matrix

  rplots <- lapply(names(rlist), function(r){
    p <- plotEmbedding(ArchRProj = rlist[[r]],
                         colorBy = "cellColData",
                         name = "Clusters",
                         embedding = "UMAP")
    p <- UmapPlotEnhance(ArchRProj = rlist[[r]],
                         p1 = p,
                         addLabels = T,
                         title.size = 16,
                         pal =  cluster_palette_list,
                         title = paste("resolution",r))

  })
  names(rplots) <- names(rlist)

  png(filename = file.path(plot.dir,paste0(file.prefix,'_clusters_assess_spread_',spread,'.png')),
      width = 28, height = 18, units = 'in', res = 200)
  do.call(grid.arrange, c(rplots, nrow=4))
  dev.off()

  # spread changed from default(=1) to 0.6
  # minDist default(=0.5)
  # spread can therefore be used to control the inter-cluster distances to some extent,
  # where as min_dist controls the size of the clusters.
  # https://jlmelville.github.io/uwot/abparams.html
  psample <- plotEmbedding(ArchRProj = rlist[[1]],
                             colorBy = "cellColData",
                             name = "bioNames",
                             embedding = "UMAP")
  psample <- UmapPlotEnhance(ArchRProj = rlist[[1]],
                             p1 = psample,
                             addLabels =  F,
                             title.size = 16,
                             pal =  bionames_palette_list,
                             legend.title = "",
                             filename.prefix = "Umap_bioNames")
  psample

  ## Decide best resolution
  return(list('rlist'=rlist, 'rplots'=rplots))

}

#### Integrate scATAC-seq with scRNA-seq ####

#' Patched copy of ArchR::addGeneExpressionMatrix
#'
#' Local modification of the ArchR function (MIT license, GreenleafLab/ArchR).
addGeneExpressionMatrixUpdate <- function(input = NULL,
                                          seRNA = NULL,
                                          chromSizes = getChromSizes(input),
                                          excludeChr = c("chrM", "chrY", "dsRED"),
                                          scaleTo = 10000,
                                          normalize = T,
                                          verbose = TRUE,
                                          threads = getArchRThreads(),
                                          parallelParam = NULL,
                                          force = TRUE,
                                          gtfFile = NULL,
                                          matrixName = "GeneExpressionMatrix",
                                          logFile = createLogFile("addGeneExpressionMatrix"))
{
  if (inherits(input, "ArchRProject")) {
    ArrowFiles <- getArrowFiles(input)
    allCells <- rownames(getCellColData(input))
    outDir <- getOutputDirectory(input)
    if (is.null(chromSizes)) {
      chromSizes <- getChromSizes(input)
    }
  }
  else if (inherits(input, "character")) {
    outDir <- ""
    ArrowFiles <- input
    allCells <- NULL
    if (is.null(chromSizes)) {
      chromSizes <- getChromSizes()
    }
  }
  else {
    stop("Error Unrecognized Input!")
  }
  if (!all(file.exists(ArrowFiles))) {
    stop("Error Input Arrow Files do not all exist!")
  }
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), "addGeneExpressionMatrix Input-Parameters", logFile = logFile)
  cellsInArrows <- unlist(lapply(ArrowFiles, ArchR:::.availableCells),use.names = FALSE)
  if (!is.null(allCells)) {
    cellsInArrows <- allCells
  }
  overlap <- sum(cellsInArrows %in% colnames(seRNA))/length(cellsInArrows)
  ArchR:::.logMessage("Overlap w/ scATAC (ratio) = ", round(overlap, 3), logFile = logFile, verbose = TRUE)
  if (overlap == 0) {
    stop("No overlap found with scATAC!")
  }
  splitCells <- split(cellsInArrows, stringr::str_split(cellsInArrows,
                                                        pattern = "#", simplify = TRUE)[, 1])
  overlapPerSample <- unlist(lapply(splitCells, function(x) sum(x %in% colnames(seRNA))))
  ArchR:::.logMessage("Overlap Per Sample w/ scATAC : ", paste(paste(names(overlapPerSample), round(overlapPerSample, 3), sep = "="), collapse = ","),
              logFile = logFile, verbose = TRUE)
  assay(seRNA) <- Matrix::Matrix(assay(seRNA), sparse = TRUE)
  nUMI <- Matrix::colSums(assay(seRNA))
  mb <- assay(seRNA)
  mb@x[mb@x > 0] <- 1
  nGenes <- Matrix::colSums(mb)
  rm(mb)
  ribo_genes <- fread(file.path(annot.dir, 'mouse_ribosomal_genes.txt'), header = F)
  ribo_genes <- paste(paste0("^",unname(unlist(ribo_genes))), collapse = "|")
  RiboRatio <- Matrix::colSums(assay(seRNA)[grep(ribo_genes, rownames(assay(seRNA))), ])/nUMI

  mito_genes <- fread(file.path(annot.dir, 'mouse_mitochondrial_genes.txt'), header = F)
  mito_genes <- paste(paste0("^",unname(unlist(mito_genes))), collapse = "|")
  MitoRatio <- Matrix::colSums(assay(seRNA)[grep(mito_genes, rownames(assay(seRNA))), ])/nUMI
  qcInfo <- DataFrame(nUMI = nUMI, nGenes = nGenes, MitoRatio = MitoRatio,
                      RiboRatio = RiboRatio)
  colnames(qcInfo) <- paste0("Gex_", colnames(qcInfo))

  ## A big fix to accomodate seRNAs that do not contain genomic range information
  if(sum(lengths(seqnames(seRNA@rowRanges)))==0){
    ArchR:::.logMessage(paste("seRNA does not contain gene coordinate information.  Trying to import with gtf..",gtfFile ),
                        logFile = logFile, verbose = TRUE)
    gtf <- importGTF(gtfFile, typeof = 'gene')

    # keeping chromosomes defined above
    gtf$seqnames <- as.character(gtf$seqnames)
    gtf <- gtf[which(as.character(gtf$seqnames) %in% as.character(seqnames(chromSizes)))]
    gtf <- gtf[which(as.character(gtf$seqnames) %ni% excludeChr)]

    gtfmatch <- match(names(seRNA@rowRanges), gtf$gene_name)
    gtfabsent <- which(is.na(gtfmatch))
    if( length(gtfmatch)>0 ){
      message("A total of ",length(gtfabsent)," out of ",length(seRNA@rowRanges)," genes found in seRNA were not found in input")
      message("Removing the following genes from scRNA")
      # missing or duplicate genes (looks like pseudogenes)
      rem <- names(seRNA@rowRanges[gtfabsent])
      message(paste(rem, collapse = ","))

      gtfmatch <- gtfmatch[which(!is.na(gtfmatch))]
      gtf <- gtf[gtfmatch]

      seRNA <- seRNA[!names(seRNA) %in% rem]
      gr <- GRanges(seqnames = gtf$seqnames, ranges = IRanges(start = gtf$start, end = gtf$end), gene_name=gtf$gene_name)
      names(gr) <- gr$gene_name
      gr <- keepSeqlevels(gr, unique(seqnames(gr)))

      if(identical(names(seRNA@rowRanges), names(gr))){
        rowRanges(seRNA) <- gr
      } else {
        stop("There was an error when trying to copy ranges to seRNA. Names did not match!")
      }
    } else{
      stop("There was an error when trying to copy ranges to seRNA. Names did not match!")
    }
    message("Successfully added gene ranges")
  }else {
    seRNA <- seRNA[BiocGenerics::which(seqnames(seRNA) %bcin% seqnames(chromSizes))]
    seRNA <- seRNA[BiocGenerics::which(seqnames(seRNA) %bcni% excludeChr)]
  }

  ## It seems there is a sorting issue below.  So we do this beforehand
  features <- sort(sortSeqlevels(rowRanges(seRNA)), ignore.strand = TRUE)
  seRNA <- seRNA[names(features)]

  idxDup <- which(rownames(seRNA) %in% rownames(seRNA[duplicated(rownames(seRNA))]))
  names(idxDup) <- rownames(seRNA)[idxDup]
  if (length(idxDup) > 0) {
    dupOrder <- idxDup[order(Matrix::rowSums(assay(seRNA[idxDup])),
                             decreasing = TRUE)]
    dupOrder <- dupOrder[!duplicated(names(dupOrder))]
    seRNA <- seRNA[-as.vector(idxDup[idxDup %ni% dupOrder])]
  }
  features <- rowRanges(seRNA)
  features <- sort(sortSeqlevels(features), ignore.strand = TRUE)
  features <- split(features, seqnames(features))
  features <- lapply(features, function(x) {
    mcols(x)$idx <- seq_along(x)
    return(x)
  })
  features <- Reduce("c", features)
  rowData(seRNA)$idx <- features[rownames(seRNA)]$idx
  ArchR:::.logThis(qcInfo, "qcInfo", logFile = logFile)
  args <- mget(names(formals()), sys.frame(sys.nframe()))
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  #args$FUN <- ArchR:::.addGeneExpressionMat
  args$FUN <- addGeneExpressionMatUpdate
  args$registryDir <- file.path(outDir, "addGeneExpressionMatRegistry")
  args$qcInfo <- qcInfo
  args$seRNA <- seRNA
  args$input <- NULL
  args$chromSizes <- NULL
  args$gtfFile <- NULL
  outList <- ArchR:::.batchlapply(args)
  ArchR:::.endLogging(logFile = logFile)
  message("Adding qc info")
  if (inherits(input, "ArchRProject")) {
    qcInfo <- qcInfo[rownames(qcInfo) %in% input$cellNames, ]
    message("Set qc info")

    #return(qcInfo)

    for (i in seq_len(ncol(qcInfo))) {
      colName <- colnames(qcInfo)[i]
      message("Adding column:",colName)
      input <- addCellColData(ArchRProj = input, data = as.vector(qcInfo[, i]), name = colName, cells = paste0(rownames(qcInfo)), force = force)
    }
    return(input)
  }
  else {
    return(unlist(outList))
  }
}

#' Per-sample worker of addGeneExpressionMatrixUpdate()
#'
#' Local modification of an ArchR internal (MIT license, GreenleafLab/ArchR).
addGeneExpressionMatUpdate <- function (i = NULL, ArrowFiles = NULL, seRNA = NULL, qcInfo = NULL,
          excludeChr = NULL, scaleTo = NULL, normalize = FALSE, cellNames = NULL, allCells = NULL,
          matrixName = "GeneExpressionMatrix",
          tstart = NULL, subThreads = 1, force = FALSE, verbose = TRUE,
          logFile = NULL)
{
  ArrowFile <- ArrowFiles[i]
  sampleName <- ArchR:::.sampleName(ArrowFile)
  o <- h5closeAll()
  o <- ArchR:::.createArrowGroup(ArrowFile = ArrowFile, group = matrixName,
                         force = force, logFile = logFile)
  if (is.null(tstart)) {
    tstart <- Sys.time()
  }
  if (is.null(cellNames)) {
    cellNames <- ArchR:::.availableCells(ArrowFile)
  }
  if (!is.null(allCells)) {
    cellNames <- cellNames[cellNames %in% allCells]
  }
  cellNames <- cellNames[cellNames %in% colnames(seRNA)]
  seRNA <- seRNA[, cellNames]
  dfParams <- data.frame(scaleTo = ifelse(is.null(scaleTo),'Null',scaleTo), exclude = excludeChr,
                         stringsAsFactors = FALSE)
  featureDF <- data.frame(seqnames = paste0(seqnames(seRNA)),
                          idx = mcols(seRNA)$idx, start = start(seRNA), end = end(seRNA),
                          name = rownames(seRNA), strand = as.integer(strand(seRNA)),
                          stringsAsFactors = FALSE)
  ArchR:::.logThis(featureDF, "featureDF", logFile = logFile)
  o <- ArchR:::.initializeMat(ArrowFile = ArrowFile, Group = matrixName,
                      Class = "double", Units = "NormCounts", cellNames = cellNames,
                      params = dfParams, featureDF = featureDF, force = TRUE)
  if(normalize){
    ArchR:::.logMessage("normalizing expression data", logFile = logFile, verbose = TRUE)
    assay(seRNA) <- ArchR:::.normalizeCols(assay(seRNA), scaleTo = scaleTo)
  } else{
    ArchR:::.logMessage("Skipping expression data normalization", logFile = logFile, verbose = TRUE)
  }

  uniqueChr <- unique(featureDF$seqnames)
  for (z in seq_along(uniqueChr)) {
    o <- tryCatch({
      o <- h5closeAll()
      chr <- uniqueChr[z]
      matz <- assay(seRNA[BiocGenerics::which(seqnames(seRNA) == chr), ])
      ArchR:::.logDiffTime(sprintf("Adding %s to %s for Chr (%s of %s)!",
                           sampleName, matrixName, z, length(uniqueChr)),
                   tstart, verbose = verbose, logFile = logFile)
      o <- ArchR:::.addMatToArrow(mat = matz, ArrowFile = ArrowFile,
                          Group = paste0(matrixName, "/", chr), binarize = FALSE,
                          addColSums = TRUE, addRowSums = TRUE, addRowVarsLog2 = FALSE) #addRowVarsLog2 = TRUE (default)
      gc()
      if (z%%3 == 0 | z == length(uniqueChr)) {
        gc()
      }
    }, error = function(e) {
      errorList <- list(ArrowFile = ArrowFile, chr = chr,
                        mat = if (exists("matz", inherits = FALSE)) matz else "matz")
      ArchR:::.logError(e, fn = ".addGeneExpressionMatUpdate AddToArrow",
                info = sampleName, errorList = errorList, logFile = logFile)
    })
  }
  allCells <- ArchR:::.availableCells(ArrowFile, passQC = FALSE)
  qcInfoi <- qcInfo[rownames(qcInfo) %in% colnames(seRNA),
  ]
  for (i in seq_len(ncol(qcInfo))) {
    infoi <- rep(-1, length(allCells))
    names(infoi) <- allCells
    infoi[rownames(qcInfoi)] <- qcInfoi[, i]
    o <- h5closeAll()
    h5write(infoi, file = ArrowFile, paste0("Metadata/",
                                            colnames(qcInfoi)[i]))
    o <- h5closeAll()
  }
  ArrowFile
}

#### Bug in filter doublets
#' Patched copy of ArchR::filterDoublets
#'
#' Local modification of the ArchR function (MIT license, GreenleafLab/ArchR).
#' With `dryrun = TRUE` only reports the number of cells that would be removed.
filterDoubletsFix <- function(ArchRProj = NULL, cutEnrich = 1, cutScore = -Inf, filterRatio = 1, dryrun=FALSE){

  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for (i in seq_along(fn)) {
    tryCatch({
      eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
  }

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = cutEnrich, name = "cutEnrich", valid = c("numeric"))
  .validInput(input = cutScore, name = "cutScore", valid = c("numeric"))
  .validInput(input = filterRatio, name = "filterRatio", valid = c("numeric"))

  if(any(grepl("filterDoublets", names(ArchRProj@projectSummary)))){
    stop("Already ran filterDoublets on ArchRProject! Cannot be re-ran on an ArchRProject!")
  }

  df <- getCellColData(ArchRProj, c("Sample", "DoubletEnrichment", "DoubletScore"))
  splitDF <- split(seq_len(nrow(df)), as.character(df$Sample))

  cellsFilter <- lapply(splitDF, function(y){

    x <- df[y, ,drop = FALSE]

    n <- nrow(x)

    x <- x[order(x$DoubletEnrichment, decreasing = TRUE), ]

    if(!is.null(cutEnrich)){
      x <- x[which(x$DoubletEnrichment >= cutEnrich), ]
    }

    if(!is.null(cutScore)){
      x <- x[which(x$DoubletScore >= cutScore), ]
    }

    if(nrow(x) > 0){
      head(rownames(x), filterRatio * n * (n / 100000))
    }else{
      NULL
    }

  }) %>% unlist(use.names=FALSE)

  message("Total ", length(cellsFilter), " cells to be filtered from ArchRProject!")
  tabRemove <- table(df[cellsFilter,]$Sample)
  tabAll <- table(df$Sample)
  samples <- unique(df$Sample)
  for(i in seq_along(samples)){
    if(!is.na(tabRemove[samples[i]])){
      message("\t", samples[i], " : ", tabRemove[samples[i]], " of ", tabAll[samples[i]], " (", round(100 * tabRemove[samples[i]] / tabAll[samples[i]], 1),"%)")
    }else{
      message("\t", samples[i], " : ", 0, " of ", tabAll[samples[i]], " (0%)")
    }
  }

  if(length(cellsFilter) > 0){

    if(dryrun){
      message("Dry run, not filtering. Returning cells to filter")
      return(df[cellsFilter,])
    }

    else
      ArchRProj@cellColData <- ArchRProj@cellColData[rownames(ArchRProj@cellColData) %ni% cellsFilter,,drop=FALSE]

  }

  ArchRProj

}

#' Collapse highly similar motifs into non-redundant groups
#'
#' Clusters the motif PWMs by similarity and keeps one representative per
#' cluster.
#'
#' @param ArchRProj ArchR project.
#' @param motifs PWM list (default: chromVARmotifs mouse_pwms_v2).
#' @param similarityMethod Similarity metric for the PWM comparison.
#' @param combine.dist.thresh Distance below which motifs are merged
#'   (1 = perfect correlation for Pearson, 0 = identity for Euclidean).
NonRedundantMotifs <- function(ArchRProj,
                               motifs=chromVARmotifs::mouse_pwms_v2,
                               similarityMethod='Pearson',
                               combine.dist.thresh=1,
                               plot.cor=T){
  # combine.dist.thresh (i.e. 1 with similarityMethod='Pearson' means perfect correlation,
  # 0 with method='Euclidean' means perfect identity)
  library(doMC)
  registerDoMC(4)
  library(foreach)
  library(TFBSTools)

  # Make a list of motif matrices
  iMotifList <- lapply(names(motifs), function(x) motifs[[x]]@profileMatrix)

  # Rename based on ArchR
  motifnames <- rownames(ArchRProj@peakAnnotation$Motif$motifSummary)[
    match(names(motifs), ArchRProj@peakAnnotation$Motif$motifSummary$ID)
  ]

  if(length(motifnames) == length(iMotifList)){
    message("Using motif names from peakAnnotation")
  } else{
    message("Using motif names from motif argument")
    motifnames <- names(motifs)
  }
  names(iMotifList) <- motifnames
  iMotifList <- iMotifList[sort(names(iMotifList))]

  .scoreMatrix <- function(m, method){
    start_time <- Sys.time()

    x <- matrix(0, length(m), length(m))
    x <-
      foreach(b=m, .combine='cbind') %:%
      foreach(a=m, .combine='c') %dopar% {
        # frobenius norm = norm(iMotifList[[1]], type='f')  = sqrt(sum(iMotifList[[1]]^2))
        if(method=='f')
          norm(x = a, type = 'f') - norm(x = b, type = 'f') # frobenius norm difference
        else
          PWMSimilarity(a, b, method = method)
      }
    rownames(x) <- colnames(x) <- names(m)
    for(i in colnames(x)) x[i,i] <- NA

    message("score matrix execution time: ", round(Sys.time() - start_time, digits = 3)," mins.")

    x
  }

  message("Scoring PWM similarity using ",similarityMethod)
  x <- .scoreMatrix(m = iMotifList, method = similarityMethod)

  # These are close matches of another
  ## idx <- which( abs(x) >= combine.dist.thresh, arr.ind = T)
  idx <- which( !is.na(x), arr.ind = T)
  dt <- data.table('TF1'=rownames(x)[idx[,1]], "TF2"=colnames(x)[idx[,2]])
  dt[, dist:= x[TF1, TF2], by=list(1:nrow(dt))]
  dt <- dt[!is.na(dist)]
  if(similarityMethod %in% c('Pearson','KL'))
    dt <- dt[dist >= combine.dist.thresh]
  else
    dt <- dt[ abs(dist) <= combine.dist.thresh]

  dt <- dt[order(TF1)]

  message("Combining entries")

  str <- foreach(i=1:nrow(dt)) %dopar% {
    tf <- dt[i,TF1]
    c(tf, dt[TF1==tf, TF2]) %>% sort %>% unique()
  }

  message("Removing duplicates entries")
  # remove duplicate pairings
  dups <- which(duplicated(str))
  if(length(dups))
    str <- str[-dups]

  # remove pairings that are subset of another
  message("Retaining largest subset entries")
  rem <- sapply(seq_along(str), function(i){
    sapply(seq_along(str), function(j){
      if(i==j)
        return()
      else if(all(str[[i]] %in% str[[j]]))
        i
    }) %>% unlist
  }) %>% unlist %>% unique

  if(length(rem))
    str <- str[-rem]

  if(any(duplicated(unlist(str))))
    stop("Problem with removing duplicates\nThis is due to combine cutoff not returning unique entries\n")

  # package this into a new list for use in addMotifAnnotations
  indep <- setdiff(rownames(x), unlist(str))
  Itf <- iMotifList[indep]

  # representative PWM for combined motifs
  repmotif <- sapply(str, function(i) i[1])
  Dtf <- iMotifList[repmotif]
  names(Dtf) <- sapply(str, function(i) { paste(i,collapse='|')})

  message("# unique motif entries ",length(Itf))
  message("# combined motif entries ",length(Dtf))

  oMotifList <- c(Itf, Dtf)

  y <- x[c(indep,repmotif), c(indep,repmotif)]
  rownames(y) <- colnames(y) <- names(oMotifList)

  if(plot.cor){
    message("Plotting input motif correlations")
    cairo_pdf(filename = file.path(plot.dir,paste0('MotifInput_Cor_Plot_',similarityMethod,'_',combine.dist.thresh,'.pdf')), width = 40, height = 40)
    corrplot::corrplot(x, type='upper',
                       is.corr = F,
                       na.label = "-",
                       diag = T,
                       tl.cex = .2,
                       cl.cex = 5,
                       method = 'square',
                       outline = F,
                       addgrid.col = NA)
    dev.off()

    message("Plotting output motif correlations")
    cairo_pdf(filename = file.path(plot.dir,paste0('MotifOutput_Cor_Plot_',similarityMethod,'_',combine.dist.thresh,'.pdf')), width = 40, height = 40)
    corrplot::corrplot(y, type='upper',
                       is.corr = F,
                       na.label = "-",
                       diag = T,
                       tl.cex = .2,
                       cl.cex = 5,
                       method = 'square',
                       outline = F,
                       addgrid.col = NA)
    dev.off()
  }

  ## Finalize List as PWM list
  pwmList <- lapply(names(oMotifList), function(x){
    PWMatrix(ID = x, name = x, profileMatrix = oMotifList[[x]])
  })
  pwmListNR <- do.call(PWMatrixList, pwmList)
  names(pwmListNR) <- names(oMotifList)
  pwmListNR
}

#' Correlate motif accessibility (chromVAR) with TF expression
#'
#' @param ArchRproj ArchR project.
#' @param dtInputChromVar Table of chromVAR motif deviations.
#' @param vassay chromVAR (deviations) assay.
#' @param expassay Gene expression assay.
#' @param matGE Optional gene expression matrix (otherwise taken from `expassay`).
#' @param min.dev.zscore,min.expression,min.accessibililty Filters on the inputs.
correlateAccessibilityExpression <- function(ArchRproj,
                                             dtInputChromVar,
                                             vassay,
                                             expassay,
                                             matGE,
                                             min.dev.zscore=-Inf,
                                             min.expression=0,
                                             min.accessibililty=-Inf,
                                             groupby=NULL,
                                             imputeValues=F,
                                             ncores=16,
                                             build=T,
                                             cor.method='pearson',
                                             xicor.symmetric=T, # If xicor should we perform symmetric test?
                                             xlsx=NULL){
  stopifnot("cor.method must be one of spearman, pearson, xicor"=cor.method %in% c('spearman','pearson','xicor'))

  require(readxl)
  require(XICOR)

  if(!build)
    return( setDT(read_xlsx(path = file.path(spreadsheets.dir, xlsx))) )

  message("#---------------------------------- Correlatating TF and Gene -----------------------------------------")
  if(missing(vassay))
    vassay <- getMatrixFromProject(ArchRProj = ArchRproj, useMatrix = 'VierstraMatrix', binarize = F, threads = ncores)

  if(missing(expassay))
    expassay <- getMatrixFromProject(ArchRProj = ArchRproj, useMatrix = 'GeneExpressionMatrix', binarize = F, threads = ncores)

  ## make function
  DT <- copy(dtInputChromVar)
  DT[,gene:=gsub(",.*","",gene)] # pick one from list , separated
  setkey(DT, idx)

  DTcellColData <- data.table(data.frame(ArchRproj@cellColData), keep.rownames = 'cell')

  dt <- NULL

  if(missing(matGE)){
    if(imputeValues){
      # Imputed values
      message("Getting imputed values")
      matGE <- imputeMatrix(assay(expassay), getImputeWeights(ArchRproj))
      rownames(matGE) <- rowData(expassay)$name
    } else {
      message("Using Unimputed expression matrix values")
      matGE <- expassay@assays@data$GeneExpressionMatrix
      rownames(matGE) <- rowData(expassay)$name
    }
  } else {
    message("Using provided gene expression matrix values")
  }

  #1. “deviations” - A deviation is a bias-corrected measurement of how far the per-cell accessibility of a given feature (i.e motif)
  # deviates from the expected accessibility based on the average of all cells or samples.

  #2. “z-score” - The z-score, also known as a “deviation score” is the z-score for each bias-corrected
  # deviation across all cells. The absolute value of the deviation score is correlated with the per-cell read depth.
  # This is because, with more reads, you have higher confidence that the difference in per-cell accessibility of the given feature (i.e. motif) from the expectation is greater than would occur by chance.

  # deviations
  X <- vassay@assays@data$deviations[DT$motifID,]
  # z-score
  Z <- vassay@assays@data$z[DT$motifID,]
  # expression
  Y <- as.matrix(matGE[which(rownames(matGE) %in% DT$gene), ,drop=T])

  DT <- DT[DT$motifID %in% rownames(X)]
  DT <- DT[DT$gene %in% rownames(Y)]

  if( !is.null(groupby) ){
    if(!groupby %in% colnames(DTcellColData))
      stop("Error ",groupby," is not a valid value for grouping")
    message("Grouping correlations by ",groupby)

    # we will scan a motif/gene combination based on the motif 'family'
    # i.e. Tcf7l1, Tcfl2, and Tcf7 are all found within the TCF/LEF-HMG-29-0 family.
    # correlate each of these genes with each motif member in the family
    #dt <- foreach(i=DT$idx) %dopar% {
    families <- unique(DT$family)
    dt <- mclapply(families, function(f){
      #dt <- lapply(families, function(f){
      dtf <- DT[family == f]
      # for each motif
      dtm <- lapply( unique(dtf$motifID), function(m){
        # for each gene in Family
        dg <- lapply(unique( dtf$gene ), function(g){
          message("Correlating deviations of ",m," with gene ",g)
          motif.assigned.gene <- dtf[motifID==m, unique(gene)]

          x <- X[m,]
          x <- x[x > min.accessibililty]   # insist that cell accessibility score is above some threshold
          y <- Y[g,]
          y <- y[y >= min.expression]       # insist that cell expresses gene above the min
          z <- Z[m,]
          z <- z[z > min.accessibililty]   # insist that cell accessibility score is above some threshold

          # for each groupby
          result <- lapply( unique(DTcellColData[,get(groupby)]), function(j){
            sc <- DTcellColData[get(groupby) == j, cell]
            gn <- length(sc) # number of cells in group
            yn <- length(intersect(sc, names(y)))
            xn <- length(intersect(sc, names(x)))
            # cells with access. and express gene
            common <- intersect(sc, intersect(names(y), names(x)))
            sy <- y[which(names(y) %in% common)] # cell expression
            sx <- x[which(names(x) %in% common)] # cell deviations or accessibility
            sz <- z[which(names(z) %in% common)] # cell z-score deviations
            expr <- tryCatch({
              if(cor.method=='xicor'){
                # is y ~ f(x), i.e. is expression a function of accessibility?
                #t1 <- xicor(x = sx, y = sy, pvalue = T, ties = T)
                t1 <- xicor(x = sz, y = sy, pvalue = T, ties = T)
                names(t1) <- c("estimate","sd","p.value")
                out <- t1
                if(xicor.symmetric){
                  #t2 <- xicor(x = sy, y = sx, pvalue = T, ties = T)
                  t2 <- xicor(x = sy, y = sz, pvalue = T, ties = T)
                  names(t2) <- c("estimate","sd","p.value")
                  if(t2$estimate > out$estimate)
                    out <- t2
                }
                out
              }
              else # correlating cell accessibility with expression
                #cor.test(sx, sy, method = cor.method, exact=F)
                cor.test(sz, sy, method = cor.method, exact=F)
            }, warning = function(warning_condition) {
              # message(warning_condition)
              list(estimate=expr$estimate, p.value=ifelse(is.na(expr$p.value), 1, expr$p.value))
            }, error = function(error_condition) {
              # message(error_condition)
              list(estimate=NA, p.value=1)
              #stop(error_condition)
            }, finally = {
              # cleanup-code
            })
            d <- data.table(j,
                            family=f,
                            motifID=m,
                            motifGene=motif.assigned.gene,
                            geneCorrelated=g,
                            coefficient=expr$estimate,
                            cor_pval=expr$p.value,
                            avg_accessibility=mean(sx),
                            avg_expression=mean(sy),
                            avg_zscore_deviation=mean(sz),
                            ncell_group=gn,
                            ncell_both=length(common),
                            ncell_motif=xn,
                            ncell_exp=yn,
                            frac_both=length(common)/gn,
                            frac_motif=xn/gn,
                            frac_exp=yn/gn
            )
            setnames(d, 'j', eval(groupby))
            d
          }) %>% rbindlist #groupby
        }) %>% rbindlist #gene
        dg[,cor_padj := p.adjust(cor_pval), by=geneCorrelated]
      }) %>% rbindlist #motif
    }, mc.cores = ncores) # family
    #})
    dt <- rbindlist(dt)
    dt[,cor_mlog10Padj := -log(cor_padj+.Machine$double.xmin)/log(10)]
  } else {
    message("Calculating correlations for all entries")
    dt <- mclapply( DT$idx, function(i){
      f <- DT[i, family] %>% unique
      m <- DT[i, motifID] %>% unique
      g <- DT[i, gene] %>% unique
      gn <- length(DT$idx)
      # Imputed values
      # seGS <- getMatrixFromProject(ArchRProj)
      # matGS <- imputeMatrix(assay(seGS), getImputeWeights(ArchRProj))
      # Raw values
      #x <- vassay@assays@data$deviations[m,]
      #x <- vassay@assays@data$z[m,]
      #x <- x[which(x > min.dev.zscore)] # z-scores cutoff, otherwise may correlate with many 0s
      #x <- vassay@assays@data$z[m,]
      #y <- expassay@assays@data$GeneExpressionMatrix[which(rowData(expassay)$name==g),]
      #y <- matGE[g,]
      x <- X[m,]
      z <- Z[m,]
      y <- Y[g,]
      y <- y[y > min.expression]
      common <- intersect(names(y), names(x))
      sy <- y[names(y) %in% common]
      sx <- x[names(x) %in% common]
      sz <- x[names(z) %in% common]
      expr <- tryCatch({
        if(cor.method=='xicor'){
          #out <- xicor(sx, sy, pvalue = T)
          out <- xicor(sz, sy, pvalue = T)
          names(out) <- c("estimate","sd","p.value")
          out
        }
        else # correlating cell accessibility with expression
          #cor.test(sx, sy, method = cor.method, exact=F)
          cor.test(sz, sy, method = cor.method, exact=F)

      }, warning = function(warning_condition) {
        # message(warning_condition)
        list(estimate=expr$estimate, p.value=ifelse(is.na(expr$p.value), 1, expr$p.value))
      }, error = function(error_condition) {
        # message(warning_condition)
        list(estimate=NA, p.value=1)
      }, finally={
        # cleanup-code
      })
      d <- data.table(family=f,
                      motifID=m,
                      geneCorrelated=g,
                      coefficient=expr$estimate,
                      cor_pval=expr$p.value,
                      avg_accessibility=mean(sx),
                      avg_expression=mean(sy),
                      avg_zscore_deviation=mean(sz),
                      ncell_both=length(common),
                      ncell_motif=length(x),
                      ncell_exp=length(y),
                      frac_both=length(common)/gn,
                      frac_motif=length(x)/gn,
                      frac_exp=length(y)/gn)
    }, mc.cores = ncores)
    dt <- rbindlist(dt)
    dt[,cor_padj := p.adjust(cor_pval), by=geneCorrelated]
    dt[,cor_mlog10Padj := -log(cor_padj+.Machine$double.xmin)/log(10)]
  }
  dt <- dt[order(coefficient, decreasing = T)]

  message("#---------------------------------- Finished correlations -----------------------------------------")

  if(!is.null(xlsx)){
    out.file <- file.path(spreadsheets.dir, xlsx)
    message("Writing to : ",out.file)
    writexl::write_xlsx(x = dt, path = out.file)
  }

  dt
}

#' Correlate motif accessibility (chromVAR) with pseudotime
#'
#' Same as correlateAccessibilityExpression() but against a pseudotime column
#' (default `palantir_pseudotime`, computed in scanpy).
#'
#' @inheritParams correlateAccessibilityExpression
#' @param pseudotime_key Name of the pseudotime column.
correlateAccessibilityPseudotime <- function(ArchRproj,
                                             dtInputChromVar,
                                             vassay,
                                             pseudotime_key="palantir_pseudotime",
                                             matGE,
                                             min.dev.zscore=-Inf,
                                             min.expression=0,
                                             min.accessibililty=-Inf,
                                             groupby=NULL,
                                             imputeValues=F,
                                             ncores=16,
                                             build=T,
                                             cor.method='pearson',
                                             xicor.symmetric=T, # If xicor should we perform symmetric test?
                                             xlsx=NULL){
  stopifnot("cor.method must be one of spearman, pearson, xicor"=cor.method %in% c('spearman','pearson','xicor'))

  require(readxl)
  require(XICOR)

  if(!build)
    return( setDT(read_xlsx(path = file.path(spreadsheets.dir, xlsx))) )

  message("#---------------------------------- Correlating TF and Pseudotime -----------------------------------------")
  if(missing(vassay))
    vassay <- getMatrixFromProject(ArchRProj = ArchRproj, useMatrix = 'VierstraMatrix', binarize = F, threads = ncores)

  DT <- copy(dtInputChromVar)
  DT[,gene:=gsub(",.*","",gene)] # pick one from list , separated
  setkey(DT, idx)

  DTcellColData <- data.table(data.frame(ArchRproj@cellColData), keep.rownames = 'cell')
  DTcellColData <- DTcellColData[order(match(cell, colnames(vassay)))]
  if(!identical(DTcellColData$cell, colnames(vassay)))
    stop("cell names not same in motif and expression matrices")

  dt <- NULL

  #1. “deviations” - A deviation is a bias-corrected measurement of how far the per-cell accessibility of a given feature (i.e motif)
  # deviates from the expected accessibility based on the average of all cells or samples.

  #2. “z-score” - The z-score, also known as a “deviation score” is the z-score for each bias-corrected
  # deviation across all cells. The absolute value of the deviation score is correlated with the per-cell read depth.
  # This is because, with more reads, you have higher confidence that the difference in per-cell accessibility of the given feature (i.e. motif) from the expectation is greater than would occur by chance.

  # deviations
  X <- vassay@assays@data$deviations[DT$motifID,]
  # z-score
  Z <- vassay@assays@data$z[DT$motifID,]
  Y <- DTcellColData[,get(pseudotime_key)]
  names(Y) <- DTcellColData$cell

  DT <- DT[DT$motifID %in% rownames(X)]
  #DT <- DT[DT$gene %in% rownames(Y)]

  if( !is.null(groupby) ){
    if(!groupby %in% colnames(DTcellColData))
      stop("Error ",groupby," is not a valid value for grouping")
    message("Grouping correlations by ",groupby)

    # we will scan a motif/gene combination based on the motif 'family'
    # i.e. Tcf7l1, Tcfl2, and Tcf7 are all found within the TCF/LEF-HMG-29-0 family.
    # correlate each of these genes with each motif member in the family
    #dt <- foreach(i=DT$idx) %dopar% {
    families <- unique(DT$family)
    dt <- mclapply(families, function(f){
      dtf <- DT[family == f]
      # for each motif
      dtm <- lapply( unique(dtf$motifID), function(m){
        # for each gene in Family
        dg <- lapply(unique( dtf$gene ), function(g){
          message("Correlating deviations of ",m," with pseudotime")
          motif.assigned.gene <- dtf[motifID==m, unique(gene)]
          x <- X[m,]
          x <- x[x > min.accessibililty]   # insist that cell accessibility score is above some threshold
          y <- Y
          z <- Z[m,]
          z <- z[z > min.accessibililty]   # insist that cell accessibility score is above some threshold

          # for each groupby
          result <- lapply( unique(DTcellColData[,get(groupby)]), function(j){
            sc <- DTcellColData[get(groupby) == j, cell]
            gn <- length(sc) # number of cells in group
            yn <- length(intersect(sc, names(y)))
            xn <- length(intersect(sc, names(x)))
            # cells with access. and express gene
            common <- intersect(sc, intersect(names(y), names(x)))
            sy <- y[which(names(y) %in% common)] # cell expression
            sx <- x[which(names(x) %in% common)] # cell deviations or accessibility
            sz <- z[which(names(z) %in% common)] # cell z-score deviations
            expr <- tryCatch({
              if(cor.method=='xicor'){
                # is y ~ f(x), i.e. is expression a function of accessibility?
                #t1 <- xicor(x = sx, y = sy, pvalue = T, ties = T)
                t1 <- xicor(x = sz, y = sy, pvalue = T, ties = T)
                names(t1) <- c("estimate","sd","p.value")
                out <- t1
                if(xicor.symmetric){
                  #t2 <- xicor(x = sy, y = sx, pvalue = T, ties = T)
                  t2 <- xicor(x = sy, y = sz, pvalue = T, ties = T)
                  names(t2) <- c("estimate","sd","p.value")
                  if(t2$estimate > out$estimate)
                    out <- t2
                }
                out
              }
              else # correlating cell accessibility with expression
                #cor.test(sx, sy, method = cor.method, exact=F)
                cor.test(sz, sy, method = cor.method, exact=F)
            }, warning = function(warning_condition) {
              # message(warning_condition)
              list(estimate=expr$estimate, p.value=ifelse(is.na(expr$p.value), 1, expr$p.value))
            }, error = function(error_condition) {
              # message(error_condition)
              list(estimate=NA, p.value=1)
              #stop(error_condition)
            }, finally = {
              # cleanup-code
            })
            d <- data.table(j,
                            family=f,
                            motifID=m,
                            motifGene=motif.assigned.gene,
                            geneCorrelated=g,
                            coefficient=expr$estimate,
                            cor_pval=expr$p.value,
                            avg_accessibility=mean(sx),
                            avg_expression=mean(sy),
                            avg_zscore_deviation=mean(sz),
                            ncell_group=gn,
                            ncell_both=length(common),
                            ncell_motif=xn,
                            ncell_exp=yn,
                            frac_both=length(common)/gn,
                            frac_motif=xn/gn,
                            frac_exp=yn/gn
            )
            setnames(d, 'j', eval(groupby))
            d
          }) %>% rbindlist #groupby
        }) %>% rbindlist #gene
        dg[,cor_padj := p.adjust(cor_pval), by=geneCorrelated]
      }) %>% rbindlist #motif
    }, mc.cores = ncores) # family
    #})
    dt <- rbindlist(dt)
    dt[,cor_mlog10Padj := -log(cor_padj+.Machine$double.xmin)/log(10)]
  } else {
    message("Calculating correlations for all entries")
    dt <- mclapply( DT$idx, function(i){
      f <- DT[i, family] %>% unique
      m <- DT[i, motifID] %>% unique
      g <- DT[i, gene] %>% unique
      gn <- length(DT$idx)

      x <- X[m,]
      z <- Z[m,]
      y <- Y
      common <- intersect(names(y), names(x))
      sy <- y[names(y) %in% common]
      sx <- x[names(x) %in% common]
      sz <- x[names(z) %in% common]
      expr <- tryCatch({
        if(cor.method=='xicor'){
          #out <- xicor(sx, sy, pvalue = T)
          out <- xicor(sz, sy, pvalue = T)
          names(out) <- c("estimate","sd","p.value")
          out
        }
        else # correlating cell accessibility with expression
          #cor.test(sx, sy, method = cor.method, exact=F)
          cor.test(sz, sy, method = cor.method, exact=F)

      }, warning = function(warning_condition) {
        # message(warning_condition)
        list(estimate=expr$estimate, p.value=ifelse(is.na(expr$p.value), 1, expr$p.value))
      }, error = function(error_condition) {
        # message(warning_condition)
        list(estimate=NA, p.value=1)
      }, finally={
        # cleanup-code
      })
      d <- data.table(family=f,
                      motifID=m,
                      geneCorrelated=g,
                      coefficient=expr$estimate,
                      cor_pval=expr$p.value,
                      avg_accessibility=mean(sx),
                      avg_expression=mean(sy),
                      avg_zscore_deviation=mean(sz),
                      ncell_both=length(common),
                      ncell_motif=length(x),
                      ncell_exp=length(y),
                      frac_both=length(common)/gn,
                      frac_motif=length(x)/gn,
                      frac_exp=length(y)/gn)
    }, mc.cores = ncores)
    dt <- rbindlist(dt)
    dt[,cor_padj := p.adjust(cor_pval), by=geneCorrelated]
    dt[,cor_mlog10Padj := -log(cor_padj+.Machine$double.xmin)/log(10)]
  }
  dt <- dt[order(coefficient, decreasing = T)]

  message("#---------------------------------- Finished correlations -----------------------------------------")

  if(!is.null(xlsx)){
    out.file <- file.path(spreadsheets.dir, xlsx)
    message("Writing to : ",out.file)
    writexl::write_xlsx(x = dt, path = out.file)
  }

  dt
}

#' Table linking motifs to the genes encoding their TFs
#'
#' Splits composite motif names into gene symbols, maps human symbols to mouse
#' orthologs, and flags TFs that are expressed in the RNA data.
#'
#' @param proj ArchR project.
#' @param motifpwm Named list of motif PWMs.
Motif2GeneLookupDT <- function(proj, motifpwm){

  ### Gene relateing to these motifs NOT expressed ##
  dtmotif <- lapply(names(motifpwm), function(i){
    motif.names <- gsub("_.*","", toupper(i))
    x <- unlist(strsplit(motif.names, '\\+'), use.names = F)
    data.table(motifID=i, motif=x)
  }) %>% rbindlist()

  dtmotif <- dtmotif[order(motifID)]
  ## Update a few known ones
  dtmotif[grep("^ZN[0-9]+", motif, value = F),
          motif := gsub("ZN","ZNF", motif)]

  # replace NKX21 with NKX2-1
  dtmotif[grep("^NKX[0-9]+$",motif, value = F),
          motif := gsub('^(NKX[0-9]+)([0-9]+)$', '\\1-\\2', motif)]

  dtmotif[grep("^HX[A-C][0-9]+$", motif, value = F),
          motif := gsub('HX', 'HOX', motif)]

  dtmotif[grep("^PRD[0-9]+$", motif, value = F),
          motif := gsub('PRD', 'PRDM', motif)]

  dtmotif[grep('ZSCAN4', motif, value=F), motif := 'ZSCAN4D']

  # First assume motifs are human and translate human to mouse homologs
  #require(biomaRt)
  #human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  #mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
  human <- biomaRt::useMart(host="https://useast.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  results <- getBM(c("ensembl_gene_id","hgnc_symbol"), "hgnc_symbol", dtmotif$motif, human)

  h2m <- humanGeneList2MouseOrthologs(humangenes = dtmotif$motif)

  setDT(h2m)
  # Create a look up table from motif to gene
  lookup <- cbind(dtmotif, h2m)

  setnames(lookup, old = c("Mouse_symbol"), new=c("gene"))
  # Several genes were not found
  missing <- lookup[is.na(gene)]

  ## make manual adjustments to names
  # Try uniprot for lookup
  # https://www.uniprot.org/uniprot/?query=suh&sort=score
  missing[motif=='ANDR', gene:='Ar']
  missing[motif=='AP2A', gene:='TFAP2A']
  missing[motif=='AP2B', gene:='TFAP2B']
  missing[motif=='AP2C', gene:='TFAP2C']
  missing[motif=='TCFAP2A', gene:='Ap2a1']
  missing[motif=='TCFAP2B', gene:='Ap2b1']
  missing[motif=='TCFAP2C', gene:='Ap2c1']
  missing[motif=='ARI5B', gene:='Arid5b']
  missing[motif=='ATF6A', gene:='Atf6']
  missing[motif=='ATF7', gene:='Atf7']

  missing[motif=='BC11A', gene:='Bcl11a']
  missing[motif=='BHA15', gene:='Bhlha15']
  missing[motif=='BHLHB2', gene:='Bhlhe40']
  missing[motif=='BHLHB3', gene:='Bhlhe41']
  missing[motif=='BMAL1', gene:='Arntl']
  missing[motif=='BRAC', gene:='T']

  missing[motif=='CART1', gene:='Alx1']
  missing[motif=='COT1', gene:='Map3k8']
  missing[motif=='COT2', gene:='Nr2f2']

  missing[motif=='ERR1', gene:='Esrra']
  missing[motif=='ERR2', gene:='Esrrb']
  missing[motif=='ERR3', gene:='Esrrg']
  missing[motif=='EVI1', gene:='Mecom']
  missing[motif=='EWSR1-FLI1', gene:='Ewsr1']
  missing[motif=='FOXP1', gene:='Foxp1']

  missing[motif=='GCR', gene:='Nr3c1']
  missing[motif=='Hen1', gene:='Henmt1']
  missing[motif=='HINFP1', gene:='Hinfp']
  missing[motif=='HNF6', gene:='Onecut1']
  missing[motif=='HTF4', gene:='Tcf12']
  missing[motif=='ITF2', gene:='Tcf4']

  missing[motif=='MYBA', gene:='Mybl1']
  missing[motif=='NR1A4', gene:='Nr4a1']
  missing[motif=='NDF1', gene:='Neurod1']
  missing[motif=='NDF2', gene:='Neurod2']

  missing[motif=='NF2L1', gene:='Nfe2l1']
  missing[motif=='NF2L2', gene:='Nfe2l2']

  missing[motif=='NFAC1', gene:='Nfatc1']
  missing[motif=='NFAC2', gene:='Nfatc2']
  missing[motif=='NFAC3', gene:='Nfatc3']
  missing[motif=='NFAC4', gene:='Nfatc4']
  missing[motif=='NGN2', gene:='Neurog2']

  missing[motif=='P53', gene:='Trp53']
  missing[motif=='P63', gene:='Trp63']
  missing[motif=='P73', gene:='Trp73']
  missing[motif=='PO2F1', gene:='Pou2f1']
  missing[motif=='PO2F1', gene:='Pou2f1']
  missing[motif=='PO2F2', gene:='Pou2f2']
  missing[motif=='PO3F1', gene:='Pou3f1']
  missing[motif=='PO3F2', gene:='Pou3f2']
  missing[motif=='PO5F1', gene:='Pou5f1']
  missing[motif=='POU5F1P1', gene:='Pou5f1']

  missing[motif=='PRGR', gene:='Pgr']

  missing[motif=='RHOX11', gene:='Rhox11']
  missing[motif=='SMCA1', gene:='Smc1a']
  missing[motif=='SMCA5', gene:='Smc5a']
  missing[motif=='SRBP1', gene:='Srebf1']
  missing[motif=='SRBP2', gene:='Srebf2']
  missing[motif=='STA5A', gene:='Stat5a']
  missing[motif=='STA5B', gene:='Stat5b']
  missing[motif=='STF1', gene:='Pdx1']

  missing[motif=='T', gene:='T']
  missing[motif=='TAF1', gene:='Taf1']
  missing[motif=='TF2L1', gene:='Tfcp2l1']
  missing[motif=='TF7L1', gene:='Tcf7l1']
  missing[motif=='TF7L2', gene:='Tcf7l2']
  missing[motif=='TFE2', gene:='E2f1']
  missing[motif=='TRP63', gene:='Trp63']
  missing[motif=='TRP73', gene:='Trp73']
  missing[motif=='TWST1', gene:='Twist1']
  missing[motif=='TYY1', gene:='Yy1']
  missing[motif=='YY2', gene:='Mbtps2']
  missing[motif=='ZKSC1', gene:='Zkscan1']
  missing[motif=='ZSCAN4', gene:='Zscan4d']
  missing[motif=='ZSCAN4D', gene:='Zscan4d']

  ##gmex <- getMatrixFromProject(ArchRProj = proj, useMatrix = "geneExpressionMatrix")
  gmex <- ArchR:::.getFeatureDF(ArrowFiles = getArrowFiles(proj), subGroup = "GeneExpressionMatrix")
  y <- sapply(missing[is.na(gene), motif], function(x){
    grep(paste0("^",x), gmex$name, ignore.case = T, value = T) %>% paste(., collapse = ',')
  })
  missing[is.na(gene)]$gene <- y

  lookup <- rbindlist(list(lookup[!is.na(gene)] , missing))[order(motif)]
  lookup <- unique(lookup)

  expsed <- sapply(lookup$gene, function(i){
    x <- NULL
    if(grepl(",", i, fixed = T))
      x <- sapply( unlist(strsplit(i,',')), function(j){
        grep(paste0("^",j,"$"), gmex$name, ignore.case = T, value = T) %>% paste(., collapse = ',')
      })
    else
      x <- grep(paste0("^",i,"$"), gmex$name, ignore.case = T, value = T) %>% paste(., collapse = ',')
    any(nchar(x)) > 0
  })

  lookup[,expressed := unlist(expsed)]
  lookup
}

# Given -log10pvalues will return combined statistic
#' Fisher's method for combining p-values given as -log10(p)
fisherCpvalLog10 <- function(x){
  chi <- 2*sum(x*log(10))
  # (2*2.811076*log(10)) + (2*255.9578*log(10))
  df <- 2*length(x) # df = 2k
  p <- pchisq(q = chi, df = df, lower.tail = F)
  list(chisq=chi, df=df, p=p)
}

#' Stouffer's method (metap::sumz) for combining p-values
stouffer <- function(x){
  ge <- which(x >= 1)
  if(length(ge))
    x[ge] <- .99999999
  zero <- which(x==0)
  if(length(zero))
    x[zero] <- .Machine$double.xmin
  metap::sumz(x)
}

#### Plots ####
#' Per-sample QC plots (TSS enrichment vs. fragments)
#'
#' @param ArchRproj ArchR project with `bioNames`, `nFrags` and `TSSEnrichment`.
QCPlots <- function(ArchRproj){
  df <- getCellColData(ArchRproj, select = c("log10(nFrags)", "TSSEnrichment","bioNames","nFrags"))

  p <- lapply(unique(df$bioNames), function(s){
    data <- df[df$bioNames==s,]
    fdata <- data[data$nFrags >= nfragMin,]
    ncell <- nrow(fdata)
    median_frag <- median(fdata$nFrags)
    median_tss <- median(fdata$TSSEnrichment)
    ggPoint(
      x = data[,1],
      y = data[,2],
      title = paste0(s,"\n","no. cells pass filter:",ncell,
                     "\n","median fragments:",median_frag,
                     "\nmedian TSS enrichment:",median_tss),
      colorDensity = TRUE,
      continuousSet = "sambaNight",
      xlabel = "Log10 (Unique Fragments)",
      ylabel = "TSS Enrichment",
      #xlim = c(log10(500), quantile(data[,1], probs = 1)),
      xlim = c(log10(500), max(data[,1])+1),
      ylim = c(0, quantile(data[,2], probs = 1)),
      baseSize = 9,
    ) + geom_hline(yintercept = 4, lty = "dashed") +
      geom_vline(xintercept = round(log10(nfragMin), digits = 1), lty = "dashed")
  })
  names(p) <- unique(df$bioNames)
  p$Arid1a
  plotPDF(plotList = p, name='TSS-vs-Frags.pdf', ArchRProj = ArchRproj, addDOC = T)
  cairo_pdf(filename = file.path(plot.dir,'TSS-vs-Fragments.pdf'), width = 9, height = 8.5, bg = 'transparent')
  do.call(grid.arrange, c(p, nrow=2))
  dev.off()

  #### QC plots ####
  p1 <- plotGroups(
    ArchRProj = ArchRproj,
    groupBy = "bioNames",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
  )
  p2 <- plotGroups(
    ArchRProj = ArchRproj,
    groupBy = "bioNames",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
  )
  p3 <- plotGroups(
    ArchRProj = ArchRproj,
    groupBy = "bioNames",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "ridges"
  )
  p4 <- plotGroups(
    ArchRProj = ArchRproj,
    groupBy = "bioNames",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
  )

  plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics-postfilter.pdf", ArchRProj = ArchRproj, addDOC = FALSE, width = 4, height = 4)

  png(filename = file.path(plot.dir,'QC-Sample-Statistics-postfilter.png'), width = 8, height = 8, units = 'in', bg = 'transparent', res = 150)
  grid.arrange(p1,p2,p3,p4, nrow=2)
  dev.off()

  p1 <- plotFragmentSizes(ArchRProj = ArchRproj, groupBy = 'bioNames')
  p1 <- p1 + geom_vline(xintercept = 147, linetype='dotted') + annotate(geom = 'text', x = 147, y=0.05, label="147nt")
  p2 <- plotTSSEnrichment(ArchRProj = ArchRproj, groupBy = 'bioNames')
  plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = ArchRproj, addDOC = FALSE, width = 5, height = 5)

}

#' Heatmap of marker peaks with ComplexHeatmap
#'
#' @param marker.peaks Marker peak SummarizedExperiment from getMarkerFeatures().
#' @param marker.cutoff Filter expression applied to the markers.
#' @param heatmap.colors Colour ramp of the heatmap body.
DoHeatmapEnhance <- function(marker.peaks,
                             marker.cutoff="!is.na(FDR)",
                             filt.heatmap.peak = NULL,
                             label.features,
                             subset.features,
                             label.sigpeak.counts = T,
                             orientation.vertical = TRUE,
                             heatmap.colors = heatmap_palette,
                             heatmap.limits = c(-2,2),
                             plotLog2FC=F, #If true plots Log2FC, otherwise Mean (default)
                             sample.feature.ratio = 1,
                             title=NULL,
                             legend.title="Peak accessibility",
                             spreadsheets.path=NULL,
                             group.bar = TRUE,
                             group.colors = NULL,
                             group.bar.count = TRUE,
                             disp.min = -2.5,
                             show.feature.names = FALSE,
                             show_column_names = TRUE,
                             show.group.dendrogram = TRUE,
                             show.cell.dendrogram = FALSE,
                             show.feature.dendrogram = FALSE,
                             group.fontsize = 10,
                             feature.fontsize = 10,
                             legend.fontsize = 9,
                             raster = TRUE,
                             draw.lines = TRUE,
                             lines.width = NULL,
                             disp.identity.legend = TRUE,
                             clustering.within = T, # this only matters if plotting individual single cells
                             archRSort = TRUE,
                             archRSortDescending = TRUE,
                             group.bar.height = 1)
{

  if(!missing(subset.features)){
    message("Subsetting peaks for features")
    idx <- which(rowData(marker.peaks)$name %in% unlist(subset.features))
    marker.peaks <- marker.peaks[idx,]
    message("Recommend setting both archRsort and clustering.within to FALSE")
  }

  ## returns a z-scaled filtered list ##
  if(is.null(filt.heatmap.peak)){
    message("Generating matrix of peaks")
    filt.heatmap.peak <- plotMarkerHeatmap(
      seMarker = marker.peaks,
      plotLog2FC = plotLog2FC, #If true plots Log2FC, otherwise Mean (default)
      #log2Norm = !plotLog2FC, # default (log2Norm=TRUE & plotLog2FC=FALSE)
      returnMatrix = T,
      cutOff = marker.cutoff,
      limits = heatmap.limits, #default (-2,2)
      transpose = TRUE
    )
  }
  if(plotLog2FC)
    legend.title <- paste(legend.title,"log2FC")

  # Write peaks and FDR to spreadsheets #
  dtheatmap_peaks <- data.table(as.data.table(rowData(marker.peaks)),
                                "FDR"=do.call(cbind, assays(marker.peaks)$FDR),
                                "Log2FC"=do.call(cbind, assays(marker.peaks)$Log2FC),
                                "MeanDiff"=do.call(cbind, assays(marker.peaks)$MeanDiff))
  # Write sample peaks and FDR to spreadsheets #
  if(!is.null(spreadsheets.path)){
    message("Writing to spreadsheet ",spreadsheets.path)
    writexl::write_xlsx(x = dtheatmap_peaks, path = file.path(spreadsheets.path), col_names = T, format_headers = T)
  }

  message("Filtering peaks with ",marker.cutoff)
  sigpeaks <- lapply(colnames(marker.peaks), function(x){
    sig.str <- gsub("Log2FC",paste0("`Log2FC.",x,"`"),marker.cutoff) %>%
      gsub("FDR",paste0("`FDR.",x,"`"),.) %>%
      gsub("MeanDiff",paste0("`MeanDiff.",x,"`"), .)
    dtheatmap_peaks[eval(parse(text = sig.str))]
  }) %>% set_names(colnames(marker.peaks))

  sigpeakcounts <- sapply(sigpeaks, nrow)
  sigpeakcounts <- sigpeakcounts[order(sigpeakcounts, decreasing = T)]
  print("sig peak counts:")
  print(sigpeakcounts)

  # order largest to smallest
  nlabel <- sapply(names(sigpeakcounts), function(i){ i })
  if(label.sigpeak.counts){
    nlabel <- sapply(names(sigpeakcounts), function(i){ paste0(i," (n=",sigpeakcounts[i],")") })
  }
  print("Label counts:")
  print(nlabel)
  message("")

  # peaks are rows and clusters on columns
  mat <- t(filt.heatmap.peak[names(sigpeakcounts),])

  print("dim(mat):")
  print(dim(mat))
  message("")

  print("Mat data range")
  print(range(mat))
  message("")

  cluster_anno <- factor(rep(colnames(mat), nrow(mat)))
  print("cluster anno:")
  levels(cluster_anno) <- names(sigpeakcounts)

  print(paste("Levels:", levels(cluster_anno)))
  message("")

  ngroup <- length(unique(cluster_anno))

  if(sample.feature.ratio<1){
    message("Sampling ",100*sample.feature.ratio,"% of features in each group")
    d <- dim(mat)
    # r <- sample(x = rownames(mat), size = d[1]*sample.cell.ratio)
    r <- sapply( levels(cluster_anno), function(i){
      idx <- sample(x = seq_len(length.out = nrow(sigpeaks[[i]])), size = nrow(sigpeaks[[i]])*sample.feature.ratio)
      sigpeaks[[i]][idx][,paste0(seqnames,":",start,"-",end)]
    })
    print(lengths(r))
    idx <- unlist(r, use.names = F)
    # subset matrix to peaks
    mat <- mat[idx,]
    # subset cluster annotation to subsampled cells
    cluster_anno <- factor(rep(colnames(mat), nrow(mat)))
    # re assign group size
    ngroup <- length(unique(cluster_anno))
  } else {
    message("Using all cells")
  }

  # reorder matrix columns to match cluster_anno
  mat <- mat[,levels(cluster_anno)]

  if(archRSort){
    bS <- ArchR:::.binarySort(m = mat, clusterCols = TRUE)
    mat <- bS[[1]][, colnames(mat)]
    clusterCols <- bS[[2]][['order']]
    if(archRSortDescending){
      #clusterCols <- rev(clusterCols)
      clusterCols <- levels(cluster_anno)
    }
    mat <- mat[, clusterCols, drop = FALSE]
    cluster_anno <- factor(colnames(mat), levels = colnames(mat))
  }

  print("cluster Anno:")
  print(levels(cluster_anno))
  #Don't scale data as we get log2FC and a certain threshold see: plotMarkerHeatmap
  #message("Scaling data")
  #mat <- scale(mat)
  # what's the value range in the matrix
  q <- quantile(mat, c(0.1, 0.95))
  q <- c(-max(abs(q)), max(abs(q)))

  if(is.null(heatmap.colors)){
    heatmap.colors <- colorRampPalette(c("royalblue4","royalblue2","deepskyblue",
                                         "grey90","orangered3","red3","darkred"), space="rgb")(1024);
  }
  col_fun = circlize::colorRamp2(seq(floor(q[1]), ceiling(q[2]), length = length(heatmap.colors)), heatmap.colors)
  library(grid)
  library(ComplexHeatmap)
  top_anno <- NULL
  if(group.bar){

    lvls <- levels(cluster_anno)

    if(is.null(group.colors)){
      group.colors <- pals::alphabet(length(lvls))
      names(group.colors) <- lvls
    }
    else if(length(names(group.colors))==length(group.colors)){
      group.colors <- group.colors[lvls]
      cat(group.colors,"\n")
    }

    df <- data.frame(sample = lvls)
    anno.counts <- sapply(df$sample, function(x){ sigpeakcounts[x] })
    df$sample <- sapply(df$sample, function(x){ nlabel[x] })
    names(group.colors) <- sapply(names(group.colors), function(x){ nlabel[x] })

    group.bar.dist <- NULL
    if(group.bar.count){
      if(orientation.vertical){
        group.bar.dist <- anno_barplot(
          x = anno.counts,
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = "gray40"), #group.colors
          width = unit(x = group.bar.height, units = 'cm'),
        )
      } else {
        group.bar.dist <- anno_barplot(
          x = anno.counts,
          which = 'row',
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = "gray40"), #group.colors
          height = unit(x = group.bar.height, units = 'cm'),
        )
      }
    }

    group_anno <- NULL

    if(orientation.vertical){
      group_anno = HeatmapAnnotation(
        show_annotation_name = F,
        dist = group.bar.dist,
        df = df,
        col = list(sample = c(group.colors)),
        border = TRUE,
        show_legend = F)
      colnames(mat) <- sapply(colnames(mat), function(x){ nlabel[x] })
      cluster_anno <- factor(rep(colnames(mat), nrow(mat)))
    } else {
      group_anno = rowAnnotation(
        df = df,
        col = list(sample = c(group.colors)),
        gp = gpar(fontsize = group.fontsize),
        #height = unit(x = group.bar.height, units = 'cm'),
        border = TRUE,
        show_legend = F,
        show_annotation_name = F,
        dist = group.bar.dist)
      colnames(mat) <- sapply(colnames(mat), function(x){ nlabel[x] })
      cluster_anno <- factor(rep(colnames(mat), nrow(mat)))
    }

  }

  # prefer group
  dendc <- F
  column_split <- NULL
  show_column_dend <- FALSE
  cluster_column_slices <- FALSE
  cluster_rows <- FALSE

  ## We don't input single cell so clustering within is same as group clustering
  if(archRSort | nrow(mat) >= 65000){
    msg <- ifelse(nrow(mat) >= 65000,
                  "Complexheatmap cannot handle dendrogram for large matrices, turning off dendrogram.",
                  "Using archR Sort")
    message(msg)
    show.group.dendrogram <- show.cell.dendrogram <- F
  }
  else if(clustering.within){
    message("Performing within & between group clustering")
    dendc <- cluster_within_group(mat, levels(cluster_anno)) # both groups and within groups are reordered, better for discerning clusters, but looks rougher

    # For large heatmaps clustering rows with dendrogram works best
    cluster_rows <- as.dendrogram(hclust(dist(mat, method = 'binary')))
    cluster_column_slices <- TRUE
  }
  else if(archRSort==F){
    message("Performing between group clustering")
    dendc <- cluster_between_groups(mat, levels(cluster_anno)) # only reordered between groups, looks more smooth
    cluster_column_slices <- TRUE
    # For smaller heatmaps (with labels) clustering rows based on dendc works better
    cluster_rows <- TRUE
  }

  if(show.group.dendrogram){
    message("Show group dendrogram")
    column_split <-  ngroup
    show_column_dend <- T
  } else if (show.cell.dendrogram) {
    message("Show cell dendrogram")
    column_split <-  cluster_anno
    show_column_dend <- T
  }

  heatmap_legend_param <- NULL
  if(orientation.vertical){
    heatmap_legend_param <- list(title = legend.title,
                                 title_position = 'topleft',
                                 direction = 'vertical',
                                 legend_position="right",
                                 legend_gp = gpar(fontsize = legend.fontsize))
  } else {
    heatmap_legend_param = list(#at = range(mat),
      title = legend.title,
      title_position = 'topcenter',
      legend_gp = gpar(fontsize = legend.fontsize),
      legend_width = unit(3, "cm"),
      grid_height = unit(.3, "cm"),
      legend_position="bottom",
      legend_direction = 'horizontal')
  }

  feature_anno <- NULL
  if(!missing(label.features)){
    if(orientation.vertical)
      feature_anno <- rowAnnotation(foo = anno_mark(at = match(x = label.features, table = rownames(mat)), labels = label.features))
    else
      feature_anno <- HeatmapAnnotation(foo = anno_mark(at = match(x = label.features, table = rownames(mat)), labels = label.features))
  }

  message("Generating heatmap")
  ht <- NULL
  if(orientation.vertical)
    ht <- Heatmap(mat,
                  name = title,
                  row_title = marker.cutoff,
                  cluster_columns = dendc,
                  column_split = column_split,
                  show_column_dend = show_column_dend,
                  cluster_column_slices = cluster_column_slices,
                  column_title_gp = gpar(fontsize = feature.fontsize),
                  column_gap = unit(0.5, "mm"),
                  cluster_rows = cluster_rows,
                  show_row_dend = show.feature.dendrogram,
                  col = col_fun,
                  row_names_gp = gpar(fontsize = feature.fontsize),
                  column_title_rot = 90,
                  top_annotation = group_anno,
                  right_annotation = feature_anno,
                  heatmap_legend_param = heatmap_legend_param,
                  show_column_names = show_column_names,
                  show_row_names = show.feature.names,
                  use_raster = raster,
                  border = 'black',
                  raster_quality = 5)
  else
    ht <- Heatmap(t(mat),
                  name = title,
                  column_title = title,
                  row_title = marker.cutoff,
                  cluster_rows = dendc,
                  row_split = column_split,
                  show_row_dend = show_column_dend,
                  cluster_row_slices = cluster_column_slices,
                  row_title_gp = gpar(fontsize = feature.fontsize),
                  row_gap = unit(0.5, "mm"),
                  cluster_columns = cluster_rows,
                  show_column_dend = show.feature.dendrogram,
                  col = col_fun,
                  column_names_gp = gpar(fontsize = feature.fontsize),
                  row_title_rot = 90,
                  top_annotation = feature_anno,
                  right_annotation = group_anno,
                  heatmap_legend_param = heatmap_legend_param,
                  show_row_names = show_column_names,
                  show_column_names = show.feature.names,
                  use_raster = raster,
                  border = 'black',
                  raster_quality = 5)
  ht

}

#' Publication-style UMAP of an ArchR project
#'
#' Plots any cellColData column or matrix feature (e.g. gene score, motif
#' deviation) on the chosen embedding, with custom palettes, rasterised points,
#' tidy legends and a fixed panel size. Saves png/pdf files to the project's
#' Plots directory.
#'
#' @param ArchRProj ArchR project.
#' @param colorBy "cellColData" or the name of a matrix (e.g. "GeneScoreMatrix").
#' @param feature Column or feature name to colour by.
#' @param embedding Embedding name.
#' @param pal Colour palette.
UmapPlotEnhance <- function(ArchRProj,
                            p1=NULL,
                            useMatrix="VierstraMatrix",
                            colorBy = "cellColData",
                            feature,
                            embedding='UMAP',
                            pal=as.vector(ArchRPalettes$solarExtra),
                            panel.border.size=0,
                            na.color="gray80",
                            axis.label=NA,
                            addLabels=NULL,
                            outline.size=.2,
                            outline.color='black', #'grey80',
                            label.size=3,
                            panel_width_in=2.5,
                            order.random=FALSE,
                            pt.size=1,
                            pt.alpha=1,
                            break.unit=1,
                            vmax=NULL,
                            vmin=NULL,
                            title,
                            title.size=16,
                            show.legend=T,
                            legend.title="",
                            legend.size=11,
                            legend.title.size=10,
                            legend.title.rotation=90,
                            legend.height.factor=.5,
                            legend.position="right",
                            legend.justification="top",
                            legend.direction="vertical",
                            use.continuous=T,
                            scale_color_bar=T,
                            fixed_coord=T,
                            small_axis_title=T,
                            filename.prefix,
                            nrow=1,
                            ncol=1,
                            as.pdf=F){
  # Plot the integration UMAP
  require(shadowtext)
  require(ggrastr)
  require(ggh4x)
  library(gtable)

  panel_width <- unit(panel_width_in, "in")

  if(is.null(p1)){
    p1 <- plotEmbedding(ArchRProj = ArchRProj,
                        colorBy = colorBy,
                        name = feature,
                        embedding = embedding,
                        plotAs = 'points',
                        size = 0.5,
                        imputeWeights = NULL)
    if(use.continuous)
      p1$data$color <- as.numeric(p1$data$color)
    else
      p1$data$color <- as.character(p1$data$color)
  }

  dt <- data.table(copy(p1$data))
  ylabel <- gsub("Dimension ","",p1$labels$y)
  xlabel <- gsub("Dimension ","",p1$labels$x)
  ## This creates 0 valued points not in faceting window
  if(!is.na(axis.label)){
    if(!axis.label %in% colnames(dt)){
      dt[,(axis.label) := color]
    }
    if(axis.label %in% colnames(dt)){
      message("Creating 0 points")
      dt[,axis.label := dt[,get(axis.label)]]
      dt <- lapply(unique(dt$axis.label), function(i){
        ds <- copy(dt[axis.label!=i])
        dt_zero <- copy(ds)
        dt_zero$color <- NA
        dt_zero$axis.label <- i
        rbindlist(list(dt_zero, ds))
      }) %>% rbindlist
    }
  }

  if(order.random)
    dt <- dt[order(sample(color), na.last = F)]
  else if(!is.numeric(dt$color) & length(levels(dt$color)) > 0){
    dt$color <- gsub("[0-9]+-","", dt$color)
    message("Setting color levels")
    dt <- dt[order(color, na.last = F)]
    dt$color <- gsub("[0-9]+-","", dt$color)
    dt$color <- as.factor(dt$color)
    lvls <- levels(dt$color)
    cat(lvls)
    dt[,color:=factor(color, levels = lvls, ordered = T)]
  }
  else {
    message("Sorting colors as numeric")
    p1$data$color <- as.numeric(p1$data$color)
    dt <- dt[order(color, decreasing = F, na.last = F)]
  }

  median.labels <- as.data.table(dt)
  median.labels <- median.labels[,.(x=median(x),y=median(y)), by=color]
  median.labels <- median.labels[order(color)]

  p1 <- NULL

  if(all(is.na(pal))){
    p1 <- ggplot(dt, aes(x=x, y=y)) + theme_ArchR()
    outline.size <- 0
  }
  else{
    p1 <- ggplot(dt, aes(x=x, y=y, fill=color)) + theme_ArchR()
  }

  if(outline.size > 0){
    stroke <- outline.size # deintensified stroke
    p1 <- p1 + geom_point_rast(size=pt.size, shape=1, stroke=outline.size,
                               col=outline.color, alpha=pt.alpha) # standard stroke #'grey40'
  }

  if(small_axis_title){
    xmin <- min(dt[,x])
    xmax <- xmin + abs(xmin)*.2
    ymin <- min(dt[,y])
    ymax <- ymin + abs(ymin)*.3
    message("xmax",xmax)
    message("xmin",xmin)

    p1 <- p1 +
      # vertical arrow
      annotate("segment", x=xmin, y=ymin, xend=xmin, yend=ymax, linewidth=.5,
               arrow = arrow(type = "closed", length = unit(0.15, "cm")), color="black") +
      annotate("segment", x=xmin, y=ymin, xend=xmax, yend=ymin, linewidth=.5,
               arrow = arrow(type = "closed", length = unit(0.15, "cm")), color="black") +
      theme(axis.title.x =  element_text(hjust = 0, vjust = 2.4, size = 12),
            axis.title.y = element_text(hjust = 0, vjust = -1.6, size = 12))

  }

  p1 <- p1 + geom_point_rast(size=pt.size, shape=21, stroke=0, alpha=pt.alpha) +
    xlab(xlabel) + ylab(ylabel) +
    theme(text = element_text(size = 15),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=title.size, hjust = 0.5),
          plot.background = element_blank(),
          panel.border = element_rect(linewidth = panel.border.size, color='black'),
          panel.spacing = unit(0, "lines"), # in case facet
          strip.background = element_blank(),
          legend.direction = legend.direction,
          legend.position = legend.position,
          legend.justification = legend.justification,
          legend.title = element_text(angle = legend.title.rotation, size = legend.title.size),
          legend.text = element_text(size = legend.size),
          legend.background = element_blank(), # get rid of legend bg element_rect(fill = 'transparent')
          legend.box.background = element_blank(), # get rid o
          legend.key=element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(t=0,r=0,b=0,l=0)) #trbl
  if(fixed_coord)
    p1 <- p1 + coord_fixed(clip = 'off')

  if(!is.null(addLabels)){
    if(is.logical(addLabels)){
      dt[,label := color]
      median.labels <- dt[,.(x=median(x),y=median(y), color), by=label]
      median.labels <- median.labels[order(label)]
      median.labels <- unique(median.labels)
      p1 <- p1 + geom_shadowtext(data = median.labels,
                                 aes(x=x, y=y, label=label, fontface='bold', family = 'Helvetica'),
                                 color='black',
                                 bg.color='white',
                                 bg.r=0.1,
                                 label.padding=0.01,
                                 size=label.size)
    }
    else if(addLabels %in% colnames(ArchRProj@cellColData)){
      ds <- plotEmbedding(ArchRProj = ArchRProj,
                          colorBy = colorBy,
                          name = addLabels,
                          embedding = embedding,
                          plotAs = 'points',
                          size = 0.5,
                          imputeWeights = NULL)$data
      ds <- setDT(ds)
      setnames(ds, c("x","y","label"))
      ds$label <- gsub("-C[0-9]+","", ds$label)
      median.labels <- ds[,.(x=median(x),y=median(y), color="Black"), by=label]
      median.labels <- median.labels[order(label)]
      median.labels <- unique(median.labels)
      p1 <- p1 + geom_shadowtext(data = median.labels,
                                 aes(x=x, y=y, label=label, fontface='bold', family = 'Helvetica'),
                                 color='black',
                                 bg.color='white',
                                 bg.r=0.05,
                                 label.padding=0.01,
                                 size=label.size)
    }
    else
      warning(addLabels, "not found in meta.data")
  }

  ## Set dimensions of plot
  gt <- ggplotGrob(p1)
  panel <- gtable_filter(gt, "panel")
  # height is scaled compared to width
  panel_height_in <- round(panel_width_in*as.numeric(gsub("null","",panel$heights)), 1)
  # set the panel height based on width
  panel_height <- unit(panel_height_in, "in") # reduce for title

  if(is.numeric(dt$color)){
    message("continuous: limits")

    lim <- range(dt$color)
    message("Data limits: [",lim[1],":", lim[2],']')
    ## round the limit up to nearest nth
    lim <- c(floor(100*lim[1])/100, ceiling(100*lim[2])/100)

    if(!is.null(vmax)){
      lim[2] <- vmax
      dt[color > vmax, color := vmax]
    }

    my.breaks <- NULL

    if(scale_color_bar){
      message("scaling color bar")
      lim <- c(-max(abs(lim)), max(abs(lim)))
      numbreaks <- 4
      if(is.null(break.unit))
        break.unit <- ceiling(lim[2]/numbreaks)
      if(lim[2] > 1){
        lim <- c(floor(-max(abs(lim))), ceiling(max(abs(lim))))
        break.unit <- ceiling(lim[2]/numbreaks)
      }
      my.breaks <- seq(from=-break.unit*numbreaks,
                       to=break.unit*numbreaks,
                       by = break.unit)
    }

    if(!missing(pal) & !is.null(pal)){
      #names(pal) <- levels(dt$color)
      names(pal) <- levels(lim)
    } else {
      pal <- ArchRPalettes$solarExtra
    }

    message("Limits: [",lim[1],":", lim[2],']')
    if(is.null(break.unit) | is.na(break.unit)){
      numbreaks <- 4
      break.unit <- round((lim[2] - lim[1]) / numbreaks, digits = 0)
    }
    if(lim[2] < 1 & break.unit==1){
      break.unit <- 0.2
    }

    message("break unit:",break.unit)

    # if limit is smaller than break
    if(break.unit > lim[2]){ # & scale_color_bar){
      # center at 0
      my.breaks <- seq(from=floor(lim[1]/break.unit)*break.unit,
                       to=ceiling(lim[2]/break.unit)*break.unit,
                       by=break.unit)
    } else if(is.null(my.breaks)) {
      my.breaks <- seq(from=lim[1],
                       to=lim[2],
                       by = break.unit)
    }

    if(!is.null(vmin)){
      my.breaks[1] <- lim[1] <- vmin
    }

    message("Breaks: [",paste(my.breaks,collapse = ","),']')

    p1 <- p1 + scale_fill_gradientn(legend.title,
                                    guide = "colourbar",
                                    colours = pal,
                                    na.value = na.color,
                                    limits=lim,
                                    breaks=my.breaks)

    p1 <- p1 + theme(legend.title = element_text(angle = 90, size=legend.title.size))

    p1 <- p1 + guides(fill = guide_colourbar(title.position = "left",
                                             ticks = T,
                                             ticks.linewidth = .2,
                                             ticks.colour = 'black',
                                             frame.colour = 'black',
                                             frame.linewidth = 0.4,
                                             nbin = 1000, # set high to put ticks at extreme ends
                                             direction=legend.direction,
                                             barwidth = .6,
                                             barheight = panel_height*legend.height.factor))
    panel_height_in <- panel_height_in+.5
  }
  #discrete
  else{
    message("discrete")
    n <- length(unique(dt$color))
    message("discrete color : ", n)
    if(missing(pal))
      pal <- pals::brewer.accent(n = n)
    else if(is.null(pal)){
      message("No colors defined")
    }
    else{
      p1 <- p1 + scale_fill_manual(legend.title, values = pal[1:n], na.value = "gray80")
      ncol.legend <- ifelse(n>10,ifelse(n>20,3,2),1)
      p1 <- p1 + guides(fill = guide_legend(override.aes = list(size=5),
                                            ncol=ncol.legend))
      p1 <- p1 + theme(legend.box.margin=margin(-10,10,-10,-10)) #trbl

      if(ncol.legend>1)
        panel_width_in <- panel_width_in+.5
    }
  }

  message(paste(colnames(dt), collapse = "," ))
  if(axis.label %in% colnames(dt)){
    ncol=1
    if(!is.null(nrow)){
      p1 <- p1 + facet_wrap(~ axis.label, nrow = nrow)
      message("Facet on ",axis.label)
    }
    else{
      p1 <- p1 + facet_wrap(~ color, nrow = nrow)
      message("Faceting on ",color)
    }
    message("ncol=",ncol)
  }

  if(!show.legend){
    p1 <- p1 + theme(legend.position = 'none')
    panel_width_in <- panel_width_in+1
    panel_height_in <- panel_height_in+.5
  }
  else{
    panel_width_in <- panel_width_in+2
    panel_height_in <- panel_height_in+.5
  }
  if(legend.title!=""){
    p1 <- p1 + theme(legend.box.margin=margin(-10,10,-10,-10), legend.title = element_text(angle = legend.title.rotation)) #trbl
  }
  if(!missing(title)){
    p1 <- p1 + ggtitle(label = title) + theme(plot.title = element_text(size = title.size))
    panel_height_in <- panel_height_in+.5 # bit more
  }

  message("Dimensions: [",panel_width,",",panel_height,"]")

  p1 <- p1 + force_panelsizes(cols = panel_width, rows = panel_height)

  if(!missing(filename.prefix)){
    if(as.pdf){
      pdf(file = file.path(getOutputDirectory(ArchRProj),'Plots',paste0(filename.prefix,".pdf")),
          width = panel_width_in, height = panel_height_in, bg = 'transparent')
      print(p1)
      dev.off()
    } else{
      png(file = file.path(getOutputDirectory(ArchRProj),'Plots',paste0(filename.prefix,".png")),
          width = panel_width_in, height = panel_height_in, units = 'in', bg = 'transparent', res = 300)
      print(p1)
      dev.off()
    }
  }
  p1
}

#' UMAP coloured by local point density, split by `feature`
#'
#' @param ArchrProj ArchR project.
#' @param feature cellColData column to facet by.
#' @param pal Colour palette for the density.
UmapPlotDensity <- function(ArchrProj,
                            feature='Phase',
                            title.size=16,
                            pt.size=1,
                            panel_width_in=3.5,
                            pal=pals::kovesi.linear_kryw_5_100_c64(10), #pal=pals::viridis(10), #pal=pals::parula(10)
                            na.color="gray85",
                            contour.color = 'blue',
                            ngridpoints=500,
                            filename.prefix,
                            fixed_coord=T,
                            embedding='UMAP',
                            small_axis_title=T,
                            coord1=NULL,
                            coord2=NULL,
                            panel.border.size=0.2,
                            legend.height.factor=1,
                            legend.title="",
                            legend.size=12,
                            legend.title.size=10,
                            legend.title.position="right",
                            legend.title.justification=1,
                            legend.title.rotation=90,
                            legend.position="right",
                            legend.justification="top",
                            legend.direction="vertical",
                            nrow=1,
                            as.pdf=F){

  require(ggrastr)
  require(ggrepel)

  .get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }

  dt <- merge( data.table( ArchrProj@embeddings[[embedding]][['df']], keep.rownames = 'cellid' ),
               data.table( as.data.frame(ArchrProj@cellColData), keep.rownames = 'cellid' ), by='cellid')
  if(is.null(coord1) | is.null(coord2)){
    coord1 <- colnames(dt)[2]
    coord2 <- colnames(dt)[3]
  }
  dt[,temp_features := get(feature)]

  dt <- dt[!is.na(get(coord1))]

  require(ggh4x)
  require(gtable)

  panel_height_in <<- 0
  panel_height <<- 0
  panel_width <<- panel_width_in

  dt[,x:=get(coord1)]
  dt[,y:=get(coord2)]
  fsort <- sort(unique(dt[,temp_features]))
  p <- lapply(fsort, function(s){
    # select data of interest
    toplot <- dt[temp_features == s]
    if(nrow(toplot) > 1){
      # increases space between plots but also allows for contours not being cutoff
      #xrange <- range(dt[,get(coord1)])
      #yrange <- range(dt[,get(coord2)])

      toplot[,density := .get_density(x = get(coord1), y= get(coord2), n=ngridpoints)]
      toplot[,density := density/max(density)]
      toplot <- toplot[order(density, decreasing = F)]
      palx <- pal
      if(is.list(pal)){
        palx <- pal[[s]]
      }

      g <- ggplot(toplot, aes(x=get(coord1), y=get(coord2), color=density)) +
        geom_point_rast(data = dt[temp_features != s], size = pt.size, fill = na.color, shape=21, color='black', stroke=0.5) + # include all points and color gray
        geom_point_rast(data = dt[temp_features != s], size = pt.size, color = na.color, shape=1) + # include all points and color gray
        geom_point_rast(size = pt.size, stroke=0.5) +
        #scale_x_continuous(limits=xrange) +
        #scale_y_continuous(limits=yrange) +
        scale_color_gradientn(legend.title, guide = "colourbar", colours = palx,
                              limits=c(0,1),
                              breaks=seq(from=0, to=1, by=.2)) +
        #ggtitle(paste("Density", gsub("_"," ",s) )) +
        ggtitle(gsub("_"," ",s)) +
        xlab( gsub(".*#|_Dimension","",coord1) ) + ylab( gsub(".*#|_Dimension","",coord2) )

      if(!is.null(contour.color))
        if(!is.na(contour.color))
          g <- g + geom_density_2d(data=toplot, linewidth=.2, color=contour.color) #pal[length(pal)])

      if(small_axis_title){
        xmin <- dt[,min(x)]
        #xmin <- xmin -(abs(xmin)*.25)
        xmax <- xmin + abs(xmin)*.3
        xlen <- abs(xmax-xmin)

        ymin <- dt[,min(y)]
        #ymin <- ymin -(abs(ymin)*.25)
        ymax <- ymin + abs(ymin)*.3
        ylen <- abs(ymax-ymin)

        minlength <- min(xlen,ylen)
        xmax <- xmin + minlength
        ymax <- ymin + minlength

        g <- g +
          # vertical arrow
          annotate("segment", x=xmin, y=ymin, xend=xmin, yend=ymax, linewidth=.5,
                   arrow = arrow(type = "closed", length = unit(0.15, "cm")), color="black") +
          annotate("segment", x=xmin, y=ymin, xend=xmax, yend=ymin, linewidth=.5,
                   arrow = arrow(type = "closed", length = unit(0.15, "cm")), color="black") +
          theme(axis.title.x =  element_text(hjust = 0, size = 12),
                axis.title.y = element_text(hjust = 0, size = 12))

      } else {
        g <- g + xlab("") + ylab("")
      }

      g <- g + theme(text = element_text(size = 15),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(size=title.size, hjust = 0.5),
            plot.background = element_blank(),
            panel.border = element_rect(linewidth = panel.border.size, color='black'),
            panel.background = element_blank(),
            panel.spacing = unit(0, "lines"), # in case facet
            legend.direction = legend.direction,
            legend.position = legend.position,
            legend.justification = legend.justification,
            legend.title = element_text(angle = legend.title.rotation, size = legend.title.size, vjust = legend.title.justification),
            legend.text = element_text(size = legend.size),
            legend.background = element_blank(), # get rid of legend bg element_rect(fill = 'transparent')
            legend.box.background = element_blank(), # get rid o
            legend.key=element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(t=0,r=0,b=0,l=0)) #trbl

      if(fixed_coord)
        g <- g + coord_fixed(clip = 'off')

      #  get dimensions of plot
      gt <- ggplotGrob(g)
      panel <- gtable_filter(gt, "panel")
      # height is scaled compared to width
      panel_height_in <<- panel_width_in*as.numeric(gsub("null","",panel$heights))
      # set the panel height based on width
      panel_height <- unit(panel_height_in, "in") # expa
      panel_width <- unique(panel_width_in, "in")

      g <- g + guides(color = guide_colourbar(title.position = "left",
                                              ticks = T,
                                              frame.colour = 'black',
                                              frame.linewidth = 0.2,
                                              nbin = 1000, # high to set extreme ends
                                              direction='vertical',
                                              barwidth = .5,
                                              barheight = (panel_height / ifelse(is.null(nrow), 1, nrow))*legend.height.factor))

      g <- g + force_panelsizes(cols = panel_width, rows = panel_height)
    }
  })

  p <- p[lengths(p) > 0]

  panel_width_in <- round((panel_width_in+.5) * (length(p)/nrow),2)
  panel_height_in <- round((panel_height_in+.5) * nrow, 2)

  message("dimensions:[",panel_width_in,",",panel_height_in,"]")

  if(!missing(filename.prefix)){
    if(as.pdf){
      f <- file.path(plot.dir,paste0(filename.prefix,".pdf"))
      message("Plotting to:",f)
      pdf(file = f, width = panel_width_in, height = panel_height_in, bg = 'transparent')
      grid_arrange_list_shared_legend(p, nrow = nrow, legend.position = "topright")
      dev.off()
    } else{
      f <- file.path(plot.dir,paste0(filename.prefix,".png"))
      message("Plotting [",panel_width_in,"x",panel_height_in,"] to:",f)
      png(filename = f, width = panel_width_in, height = panel_height_in, bg = 'transparent', units = 'in', res=200)
      grid_arrange_list_shared_legend(p, nrow = nrow, legend.position = "topright")
      dev.off()
    }
  }

  p
}

#' UMAPs of motif deviation, TF expression and gene score side by side
#'
#' @param proj ArchR project.
#' @param DT Table of motifs/genes to plot.
#' @param motifMatrix Name of the motif matrix.
#' @param plot.motif,plot.expression,plot.geneScore Which panels to include.
PlotMotifExpressionGScoreUmap <- function(proj,
                                          DT,
                                          motifMatrix="VierstraMatrix",
                                          plot.motif=T,
                                          plot.expression=T,
                                          plot.geneScore=F,
                                          embedding="UMAP",
                                          imputeValues=F,
                                          scale_color_bar=T,
                                          scale_color_bar_expression=F,
                                          legend.height.factor=1,
                                          group=T,
                                          plot.title=T,
                                          pt.size=0.9,
                                          nrows=NULL,
                                          pal1=solarExtraWhite_palette,
                                          pal2=solarHalf,
                                          pal3=solarHalf)
{
  if(plot.motif & !'motifID' %in% colnames(DT))
    stop("DT columns must contain column 'motifID'")
  if(plot.expression & !'gene' %in% colnames(DT))
    stop("DT columns must contain column 'gene'")
  if(plot.geneScore & !'gene' %in% colnames(DT))
    stop("DT columns must contain column 'gene'")

  imputeWeights <- NULL
  if(imputeValues){
    imputeWeights <- getImputeWeights(proj)
  }

  dt <- unique(DT)
  plotList <- list()
  genelabel <- ""
  # Plots the deviation z score (i.e.)
  # vassay <- getMatrixFromProject(ArchRProj = proj, useMatrix = 'VierstraMatrix', binarize = F, threads = 2)
  # vassay@assays@data$z
  if(plot.motif){
    motifId <- dt[, motifID] %>% unique
    genelabel <- dt[, motifID]
    p <- plotEmbedding(
      ArchRProj = proj,
      colorBy = eval(motifMatrix),
      name = paste0('z:',motifId),
      embedding = embedding,
      plotAs = 'points', size = 0.5,
      imputeWeights = imputeWeights
    )
    if(is.ggplot(p))
      p <- list(p)
    plotList[['1.motif']] <-
      lapply(p, function(i){
        UmapPlotEnhance(ArchRProj = proj,
                        p1 = i,
                        pt.size = pt.size,
                        panel.border.size = 0,
                        legend.justification = "bottom",
                        title = gsub('#.*','', gsub(".*:","",i$labels$title)),
                        label.size = 12,
                        title.size = 18,
                        legend.size = 14,
                        legend.height.factor = legend.height.factor,
                        use.continuous = T,
                        break.unit = 1,
                        scale_color_bar = scale_color_bar,
                        pal = pal1,
                        legend.title = 'Motif deviation zscore')
      })
    names(plotList[['1.motif']]) <- motifId
  }

  if(plot.expression){
    genelabel <- dt[, as.character(gene)] %>% unique()
    cat("Plotting expression\n")
    p <- plotEmbedding(
      ArchRProj = proj,
      colorBy = "GeneExpressionMatrix",
      name = genelabel,
      embedding = embedding,
      plotAs = 'points',
      log2Norm = F,
      size = 0.9,
      imputeWeights = imputeWeights
    )
    if(is.ggplot(p))
      p <- list(p)
    plotList[['2.expression']] <- lapply(p, function(i){
      UmapPlotEnhance(ArchRProj = proj,
                      p1 = i,
                      pt.size = pt.size,
                      panel.border.size = 0,
                      legend.justification = "bottom",
                      title = gsub('#.*','', gsub(".*:","",i$labels$title)),
                      label.size = 12,
                      title.size = 18,
                      legend.size = 14,
                      use.continuous = T,
                      legend.height.factor = legend.height.factor,
                      scale_color_bar = scale_color_bar_expression,
                      pal = pal2,
                      legend.title = 'scRNA normalized data')
    })
    names(plotList[['2.expression']]) <- genelabel
  }

  if(plot.geneScore){
    cat("Plotting gene score\n")
    genelabel <- dt[, as.character(gene)] %>% unique()
    p <- plotEmbedding(
      ArchRProj = proj,
      colorBy = "GeneScoreMatrix",
      name = genelabel,
      embedding = embedding,
      plotAs = 'points', size = 0.5,
      imputeWeights = imputeWeights
    )
    if(is.ggplot(p))
      p <- list(p)

    plotList[['3.genescore']] <- lapply(p, function(i){
      UmapPlotEnhance(ArchRProj = proj,
                      p1 = i,
                      pt.size = pt.size,
                      title = gsub('#.*','', gsub(".*:","",i$labels$title)),
                      label.size = 12,
                      title.size = 18,
                      legend.size = 14,
                      use.continuous = T,
                      scale_color_bar = scale_color_bar_expression,
                      legend.height.factor = legend.height.factor,
                      pal = pal3,
                      legend.title = 'Locus accessibility')
    })
    names(plotList[['3.genescore']]) <- genelabel
  }

  # group by gene
  if(group & (plot.expression | plot.geneScore)){
    cat("Grouping by gene\n")
    genelabel <- dt[, as.character(gene)] %>% unique()
    plotList <- lapply(genelabel, function(i){
      retL <- list()
      if('1.motif' %in% names(plotList)){
        retL <- plotList[['1.motif']][dt[gene==i, motifID]]
      }
      if('2.expression' %in% names(plotList)){
        retL <- c(retL, plotList[['2.expression']][i])
      }
      if('3.genescore' %in% names(plotList)){
        retL <- c(retL, plotList[['3.genescore']][i])
      }
      if(is.ggplot(retL))
        retL <- list(retL)
      do.call(grid.arrange, c(retL, nrow=nrows))
    })
    names(plotList) <- genelabel
  } else {
    plotList <- lapply(1:nrow(dt), function(i){
      cat("Plotting\n")
      x <- ""
      retL <- list()
      if('1.motif' %in% names(plotList)){
        x <- dt[i, motifID]
        retL <- plotList[['1.motif']][[x]]
      }
      if('2.expression' %in% names(plotList)){
        x <- dt[i, gene]
        retL <- c(retL, plotList[['2.expression']][[x]])
      }
      if('3.genescore' %in% names(plotList)){
        x <- dt[i, gene]
        retL <- c(retL, plotList[['3.genescore']][[x]])
      }
      if(is.ggplot(retL))
        retL <- list(retL)
      if(plot.title)
        retL[['top']] <- textGrob(x, gp=gpar(fontsize=14,font=2))
      retL
    })
  }

  message("Returning ",length(plotList)," records")

  plotList
}

#' Stacked bar chart of a cellColData variable per cluster
#'
#' @param ArchRProj ArchR project.
#' @param variable cellColData column plotted as stacked fill.
#' @param palette Named colours for `variable`.
#' @param filename Output file name.
PlotBarVariableInClusters <- function(ArchRProj,
                                      variable='Condition',
                                      palette=condition_palette,
                                      addcounts=T, filename="",
                                      orientation='horizontal',
                                      xlabel = "Cluster number",
                                      width_scale_perc=1,
                                      height_scale_perc=1){
  dt <- as.data.table(ArchRProj@cellColData[,c("Clusters",eval(variable))])

  dt[,nCluster := .N, by=Clusters]
  dt[,nClusterExperiment := .N, by=list(Clusters, get(variable))]
  dt <- unique(dt)
  g <- NULL
  if(orientation=='horizontal'){
    dt <- dt[order(nCluster, decreasing = T),]
    dt[,Clusters:=factor(x = Clusters, levels = unique(Clusters), ordered = T)]
    dt[,percClusterExperiment := 100 * (nClusterExperiment / nCluster) ]
    dtlab <- unique(dt[,list(Clusters, nCluster)])
    g <- ggplot(dt, aes(x=Clusters, y=percClusterExperiment)) +
      geom_bar(stat='identity', aes(fill=get(variable)), width = .95) +
      scale_fill_manual(eval(variable),values = palette) +
      ylab("Percentage of cells") + xlab(xlabel) +
      theme_ArchR() +
      coord_cartesian(ylim = c(0,100), expand = F, clip = 'off') +
      theme(panel.ontop = T,
            panel.grid = element_blank(),
            # panel.grid.major = element_blank(),
            # panel.grid.minor = element_blank(),
            # panel.grid.major.x = element_blank(),

            panel.grid.major.y = element_line(colour = 'black', linewidth = .2, linetype = 'longdash'),
            legend.direction = "horizontal",
            legend.position = "bottom",
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 12),
            text = element_text(size = 14),
            legend.text = element_text(size = 11),
            legend.background = element_blank(),
            plot.background = element_blank())
    if(addcounts)
      g <- g + geom_text(data = dtlab, aes(label=paste0("n=",nCluster)), y=100, color='black',
                         vjust=0, hjust=0, angle=45, size=3.5)
    if(filename !=""){
      fout <- file.path(getOutputDirectory(ArchRProj),'Plots',paste0(filename,'.pdf'))
      message("Plotting to:",fout)
      pdf(file = fout, width = 5*width_scale_perc, height = 3*height_scale_perc, bg = 'transparent')
      print(g)
      dev.off()
    }
  } else {
    dt <- dt[order(nCluster, decreasing = F),]
    dt[,Clusters:=factor(x = Clusters, levels = unique(Clusters), ordered = T)]
    dt[,percClusterExperiment := 100 * (nClusterExperiment / nCluster) ]
    dtlab <- unique(dt[,list(Clusters, nCluster)])
    g <- ggplot(dt, aes(x=percClusterExperiment, y=Clusters)) +
      geom_bar(stat='identity', aes(fill=get(variable)), width = .95) +
      scale_fill_manual(eval(variable),values = palette) +
      xlab("Percentage of cells") +  ylab(xlabel) +
      theme_ArchR() +
      coord_cartesian(xlim = c(0,100), expand = F, clip = 'off') +
      theme(panel.ontop = T,
            panel.grid = element_blank(),
            panel.grid.major.x = element_line(colour = 'black', linewidth = .2, linetype = 'longdash'),
            legend.direction = "horizontal",
            legend.position = "bottom",
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 12),
            text = element_text(size = 14),
            legend.text = element_text(size = 11),
            legend.background = element_blank(),
            plot.margin = unit(c(1,2,1,1), "cm"), #trbl
            plot.background = element_blank())
    g <- g + guides(fill=guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))
    if(addcounts)
      g <- g + geom_text(data = dtlab, aes(label=paste0("n=",nCluster)), x=100, color='black',
                         vjust=0.5, hjust=0, angle=0, size=3.5)
    if(filename !=""){
      fout <- file.path(getOutputDirectory(ArchRProj),'Plots',paste0(filename,'_vertical.pdf'))
      message("Plotting to:",fout)
      pdf(file = fout, width = 3*width_scale_perc, height = 5*height_scale_perc, bg = 'transparent')
      print(g)
      dev.off()
    }
  }
  g
}

#' Correlation of LSI dimensions with sequencing depth
PlotCorDimensionDepth <- function(ArchRProj){
  dt <- data.table(do.call(cbind, ArchRProj@reducedDims$IterativeLSI$corToDepth), keep.rownames = 'dimension')
  dt <- melt(dt, id.vars = 'dimension')
  dt[,dimension:= gsub("LSI","",dimension)]
  dt[,dimension := factor(x = dimension, levels = unique(dimension), ordered = T)]
  ggplot(dt, aes(x=dimension, y=value, color=variable)) + geom_point(size=3) + ggtitle(label = 'PCC Depth & dimension') +
    scale_color_brewer("LSI Dimension", palette = 'Set1') + xlab("LSI dimension") + ylab("PCC") + theme_bw()
}

#' Confusion heatmap of ATAC clusters versus predicted RNA clusters
PlotIntegrationHeatmap_ATAC_RNA <- function(ArchRProj, show.percentages=FALSE){
  ## Percent ATAC cells labeled
  pal <- as.vector(rna_cluster_palette)
  names(pal) <- names(rna_cluster_palette)

  dt <- data.table(table(data.frame(rna_group=as.character(ArchRProj@cellColData$gex_seurat_clusters), atac_cluster=ArchRProj@cellColData$Clusters)))

  #dt[,rna_group:=as.character(rna_group)]
  dt <- dt[rna_group!=-1]
  dt <- dt[order(atac_cluster)]
  #dt[,atac_percent:= 100*N/sum(N), by=atac_cluster]
  dt[,rna_percent:= 100*N/sum(N), by=rna_group]
  dtmat <- dcast.data.table(dt, atac_cluster ~ rna_group, value.var = 'rna_percent')
  mat <- as.matrix(dtmat[,-c('atac_cluster'), with=F])
  rownames(mat) <- dtmat$atac_cluster
  pal <- pal[!is.na(names(pal))]
  ha_column = HeatmapAnnotation(df = data.frame(rna_group=colnames(mat)),
                                show_annotation_name = F,
                                show_legend = F,
                                col = list(rna_group = pal),
                                annotation_legend_param = list(title="", labels_gp = gpar(fontsize = 14)))

  cm = cor(mat)

  od =  hclust(dist(mat))$order

  mat = mat[od, ]

  bS <- ArchR:::.binarySort(m = mat, clusterCols = TRUE)
  mat <- bS[[1]][, colnames(mat)]
  clusterCols <- 1:ncol(mat)
  clusterCols <- bS[[2]][['order']]
  mat <- mat[, clusterCols, drop = FALSE]

  default.cell.outline <- gpar(col = "gray80", lwd = 1)

  heatmap <- ComplexHeatmap::Heatmap(name = "RNA % shared in ATAC",
                                     matrix = mat,
                                     row_title = "snATAC clusters",
                                     row_title_side = "left",
                                     cluster_rows = F,
                                     cluster_columns = F,
                                     border = T,
                                     column_title = "snRNA clusters",
                                     column_title_side = "bottom",
                                     rect_gp = default.cell.outline, # creates "space" between cells
                                     #rect_gp = gpar(col = "black", lwd = .5),
                                     #col = unname(ArchR::ArchRPalettes$whitePurple),
                                     col = gist_heat_r,
                                     #bottom_annotation = ha_column,
                                     show_column_dend = F,
                                     show_row_dend = F,
                                     cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                       if(show.percentages==TRUE){
                                         ifelse(mat[i, j] > 40,
                                                shadowtext::grid.shadowtext(paste0(round(mat[i, j], digits = 1),"%"), x, y,
                                                                            bg.colour = 'grey60', gp = gpar(col='white'), bg.r = 0.05),
                                                ifelse(mat[i, j] <= 0,
                                                       grid.text("-", x, y,  gp=gpar(col='grey40')),
                                                       grid.text(paste0(round(mat[i, j]),"%"), x, y,  gp=gpar(col='grey10'))
                                                )
                                         )
                                       }
                                       grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "white", lwd = .5)) # outlines cells
                                     },
                                     heatmap_legend_param = list(at = c(0, 100),
                                                                 labels = c("0%", "100%"),
                                                                 title_position = 'leftcenter-rot',
                                                                 title_gp = gpar(fontsize = 12),
                                                                 grid_width = unit(.5, "cm"),
                                                                 border = 'black',
                                                                 legend_height = unit(4, "cm")),
                                     row_names_side = "left")
  pdf(file = file.path(getOutputDirectory(ArchRProj),'Plots','RNA_ATAC_Integration_Heatmap.pdf'),
      width = 8, height = 3.5, bg = 'transparent')
  draw(heatmap)
  dev.off()

  heatmap
}

#' Sequence logos of motifs
#'
#' @param motifList Named list of PWMs.
#' @param plotas "seqlogo" for TFBSTools::seqLogo, otherwise ggseqlogo.
PlotMotifSeqLogo <- function(motifList, motifnames=NA, plotas='ggseq'){
  library(TFBSTools)
  gl <- NULL
  if(plotas != 'seqlogo'){
    library(ggseqlogo)
    #list_fonts()
    if(is.na(motifnames)[1])
      motifnames <- names(motifList)
    gl <- lapply(seq_along(motifList), function(i){
      x <- motifList[[i]]
      title <- motifnames[i]
      if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
      mat <- exp(as(x,"matrix") - 0.008)*bg(x) / sum(bg(x))
      ggseqlogo(mat, font='roboto_medium') + ggtitle(title) +
        theme(plot.title = element_text(size=9),
              text = element_text(size=8),
              panel.background = element_rect(colour = 'transparent', fill = 'transparent'),
              plot.background = element_rect(colour = 'transparent', fill = 'transparent')) #
    })
    names(gl) <- motifnames
  }
  ## SeqLogo
  if(plotas=='seqlogo'){
    library(seqLogo)
    for(i in names(motifList)){
      x <- motifList[[i]]
      if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
      mat <- exp(as(x,"matrix") - 0.008)*bg(x) / sum(bg(x))
      icm <- ICMatrix(ID = i, name = i, profileMatrix = mat)
      seqLogo(x = icm)
    }
  }
  gl
}

#' Rank plot of motif enrichment (-log10 adjusted p-value by rank)
PlotMotifs <- function(se, topn=30){
  df <- data.frame(TF = rownames(se), mlog10Padj = assay(se)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))

  ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(topn)), ], aes(x = rank, y = mlog10Padj, label = TF),
      size = 4,
      nudge_x = 2,
      label.padding = 0.1,
      label.size = 0.1,
      color = "black"
    ) + theme_ArchR() +
    ylab("-log10(P-adj) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))

  ggUp
}

#' Tile grid of top motifs enriched in up- and down-regulated peaks
PlotMotifsTileGridSymmetric <- function(seListUp, seListDown, scale='group', title="", topn=15){

  dtup <- lapply(names(seListUp), function(i){
    se <- seListUp[[i]]
    dt <- data.table(do.call(cbind, se@assays@data)) %>% setnames(names(se@assays@data))
    #dt <- data.table(TF = rownames(se), mlog10Padj = assay(se)[,1])
    dt <- dt[order(mlog10Padj, decreasing = TRUE),]
    dt$rank <- seq_len(nrow(dt))
    dt$group <- i
    dt$direction <- 'up'
    dt
  }) %>% rbindlist

  dtdown <- lapply(names(seListDown), function(i){
    se <- seListDown[[i]]
    dt <- data.table(do.call(cbind, se@assays@data)) %>% setnames(names(se@assays@data))
    #dt <- data.table(TF = rownames(se), mlog10Padj = assay(se)[,1])
    dt <- dt[order(mlog10Padj, decreasing = TRUE),]
    dt[,mlog10Padj := -1*mlog10Padj]
    dt$rank <- seq_len(nrow(dt))
    dt$group <- i
    dt$direction <- 'down'
    dt
  }) %>% rbindlist

  # Top n by group
  #dtup <- dtup[ abs(mlog10Padj) >= -log10(0.05), head(.SD, topn), by=group]
  #dtdown <- dtdown[ abs(mlog10Padj) >= -log10(0.05), head(.SD, topn), by=group]

  # Let's filter later
  dtup <- dtup[abs(mlog10Padj) >= -log10(0.05)]
  dtdown <- dtdown[ abs(mlog10Padj) >= -log10(0.05)]

  dt <- rbindlist(list(dtup, dtdown))

  # It's possible that a motif is significantly enriched up and down in these comparisons
  # Makes me think the tests are a bit off
  # In any case, if it's found sig. in both directions then we remove it
  # group by 'feature', 'group' and if the nrow is equal to 1 (.N==1), we get the Subset of Data.table (.SD).
  message("Removing any potential duplicates")
  dt <- dt[, if(.N==1) .SD, .(feature, group)]

  # Top n features by group
  dt <- dt[order(rank, decreasing = FALSE)]
  topfeatures <- dt[, head(.SD, topn), by=list(group, direction)]$feature
  #dt <- dt[, head(.SD, topn), by=list(group, direction)]
  dt <- dt[feature %in% topfeatures]

  wide.filt <- dcast.data.table(dt, group ~ feature, value.var = 'mlog10Padj')
  m.wide.filt <- melt(wide.filt, id.vars = c("group"))
  m.wide.filt[,abs_max_val := max(abs(value), na.rm = T), by=variable]
  m.wide.filt <- m.wide.filt[order(abs_max_val, decreasing = T)]
  m.wide.filt[,variable := factor(variable, unique(variable), ordered = T)]

  legend.text <- 'Motif enrichment'
  if(!missing(scale)){
    if(scale=='group')
      m.wide.filt[,value := scale(value, center = T, scale = T), by=list(group)]
    else if(scale=='motif')
      m.wide.filt[,value := scale(value, center = T, scale = T), by=list(variable)]
    legend.text <- paste(legend.text,'\n',scale,'scaled')
  }

  scale_lim <- m.wide.filt[,max(abs(value), na.rm = T)]

  if(missing(title))
    title <- paste("Top",topn,"motifs per pairwise comparison")

  ggplot(m.wide.filt, aes(x=variable, y=group, fill=value)) +
    geom_point(shape=22, size=5, stroke=.1) +
    # scale_fill_gradientn(expression(-log[10](padj)),
    #                      colors = paletteContinuous('coolwarm'),
    #                      na.value = 'white',
    #                      limits=c(-scale_lim, scale_lim)) +
    scale_fill_gradientn(legend.text,
                         colors = rev(pals::brewer.rdbu(10)),
                         na.value = 'transparent',
                         breaks=c(-scale_lim, 0, scale_lim),
                         labels=c("Sig. down",0,"Sig. Up"),
                         limits=c(-scale_lim, scale_lim)) +
    coord_fixed() +
    ggtitle(label = title) +
    ylab("Pairwise comparison") +
    xlab(paste("Top",topn,"motifs by group")) +
    theme_bw() + theme(text = element_text(family = 'Arial', size = 14, colour = 'black'),
                       panel.grid = element_blank(),
                       axis.line = element_blank(),
                       axis.ticks = element_blank(),
                       panel.border = element_blank(),
                       plot.background = element_blank(),
                       legend.background = element_blank(),
                       # legend.position = 'right',
                       # legend.direction = 'vertical',
                       # legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                       legend.title = element_text(angle = 0),
                       legend.position = 'bottom',
                       legend.direction = 'horizontal',
                       plot.title = element_text(hjust = 0.5, size=16),
                       plot.margin = unit(x=c(.1,.1,.1,.1), units = 'cm'),
                       axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    coord_fixed() +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, title.vjust = 0.5, barwidth = unit(3.5, units = 'cm')))

}

#' Single-cell heatmap of motif deviation scores
#'
#' @param enrichMotifs Motif enrichment result.
#' @param group.by cellColData column used to group cells.
plotSingleCellMotifHeatmap <- function(enrichMotifs,
                                       cellset=NULL,
                                       cellset_ordered=F,
                                       group.by=NULL,
                                       sample.cell.ratio=1,
                                       mlog10Padj.cutOff=0,
                                       labelMotifs=T,
                                       labelMotifN=F,
                                       labelMotifFrac=F,
                                       labelSpecificMotifs,
                                       labelMotifClusters,
                                       remove.redundant=T,
                                       no.clustering=F,
                                       clustering.within=F,
                                       cluster.columns=T,
                                       cluster.rows=F,
                                       cluster.group.slices=F,
                                       heatmap.gap.mm = 0.5,
                                       add.score.to.name=T,
                                       col=unname(ArchRPalettes$whitePurple),
                                       na.col='gray90',
                                       show.group.dendrogram=T,
                                       show.cell.dendrogram=F,
                                       group.annotation.bar=T,
                                       group.bar.height=.5,
                                       group.colors = cluster_palette_list,
                                       pseudotimes = NULL,
                                       pseuodtime_cutoff = c(0,1),
                                       pseudotime_palette = pals::brewer.blues,
                                       reverse_pseudotime_palette = F,
                                       show.group.bar.legend=F,
                                       quantile.min=.1,
                                       quantile.max=.90,
                                       label.col=pals::polychrome(36),
                                       orientation_horizontal=T,
                                       legend.ncol=NA,
                                       title.legend ="Accessibility",
                                       title = "Enriched transcription factor motifs per family",
                                       row_title = "",
                                       row_name_fontsize=12,
                                       column_name_fontsize=7,
                                       outline=T,
                                       use.raster=T,
                                       scale_data=NA,
                                       label.outline=0.5){

  require(circlize)
  require(ComplexHeatmap)

  default.cell.outline <- gpar(col = NA)

  mat <- NULL
  if ('matrix' %in% class(enrichMotifs)){
    message("Using matrix")
    mat <- enrichMotifs
  }
  else if('SummarizedExperiment' %in% class(enrichMotifs))
    mat <- as.matrix(assays(enrichMotifs)[["mlog10Padj"]])
  else if ('data.table' %in% class(enrichMotifs)){
    mat <- as.matrix(enrichMotifs, rownames='motifID')
    mat <- mat[,grep("mlog10Padj", colnames(mat))]
    colnames(mat) <- gsub('mlog10Padj','',colnames(mat))
  } else if (class(enrichMotifs) %in% c('dgCMatrix','dgeMatrix')){
    message("Getting sparse matrix")
    mat <- as.matrix(enrichMotifs)
  }

  # remove empty
  ridx <- which(colSums(mat)==0)
  if(length(ridx)>0){
    message("Removing empty")
    mat <- mat[,-ridx]
  }

  mat[is.infinite(mat)] <- NA

  if(!is.null(cellset)){
    mat <- mat[,cellset]
  }

  group_anno<- factor(proj@cellColData[colnames(mat),group.by])
  ngroup <- length(unique(group_anno))
  if(sample.cell.ratio<1 & !is.null(cellset)){
    message("Sampling ",100*sample.cell.ratio,"% of cells in each group")
    d <- dim(mat)
    r <- sapply( levels(group_anno), function(i){
      sample(x = which(group_anno == i), size = length(which(group_anno == i))*sample.cell.ratio)
    })
    print(table(group_anno))
    print(lengths(r))
    idx <- unlist(r, use.names = F)
    # subset matrix to subsampled cells
    mat <- mat[, idx]
    # subset cluster annotation to subsampled cells
    group_anno <- group_anno[,idx]
    # re assign group size
    ngroup <- length(unique(group_anno))
  } else {
    message("Using all cells")
  }
  message("clusters:", paste(levels(group_anno), collapse = ", "))

  val.range <- c(0,1)

  if(is.null(pseudotimes)){
    q <- c(min(mat), max(mat))
    message("mat range:",paste(range(mat), collapse = ":"))
    message("lower-upper limit:",paste(round(q,2), collapse = ":"))

    message("capping color values at [",quantile.min," and ",quantile.max,"] quantiles")
    q <- quantile(mat, c(quantile.min, quantile.max))
    # what's the value range in the matrix
    message("Clipping and centering the range")
    val.range <- c(-max(abs(q)),max(abs(q)))
    cat("range quantile = ",val.range,"\n")
  } else {
    idx <- sort(which(pseudotimes >= pseuodtime_cutoff[1] & pseudotimes <= pseuodtime_cutoff[2]), decreasing = F)
    mat <- mat[,idx]
    message("Subsetting matrix using pseudotimes (n=",length(idx),")")
    cellset <- colnames(mat)
    pseudotimes <- pseudotimes[idx]

    pseudotime_palette <- pseudotime_palette(length(pseudotimes))
    if(reverse_pseudotime_palette)
      pseudotime_palette <- rev(pseudotime_palette)
    pseudotime_palette <- colorRamp2(seq(pseuodtime_cutoff[1], pseuodtime_cutoff[2], length = ncol(mat)), colors = pseudotime_palette)

    if(length(scale_data)>1){
      message("Scaling data between [",paste(scale_data, collapse = "-"),"]")
      mat <- t(apply(mat, 1, function(x) { scales::rescale(x, to=scale_data) }))
    }
    #val.range <- c(0,1)
  }

  group_anno<- factor(proj@cellColData[colnames(mat),group.by])
  bionames_anno<- factor(proj@cellColData[colnames(mat),'bioNames'])
  cluster_anno <- factor(proj@cellColData[colnames(mat),'Clusters'])
  condition_anno <- factor(proj@cellColData[colnames(mat),'Condition'])
  cellcluster_anno <- factor(proj@cellColData[colnames(mat),'gex_mammary_gland_cluster_classification'])
  celltype_anno <- factor(proj@cellColData[colnames(mat),'gex_mammary_gland_refine_clusters'])

  message("group:", paste(levels(group_anno), collapse = ", "))
  message("bionames:", paste(levels(bionames_anno), collapse=", "))
  message("clusters:", paste(levels(cluster_anno), collapse = ", "))
  message("condition:", paste(levels(condition_anno), collapse = ", "))
  message("cell clusters:", paste(levels(cellcluster_anno), collapse = ", "))
  message("cell types:", paste(levels(celltype_anno), collapse = ", "))

  heatmap.legend.ticks <- c(val.range[1],val.range[2])
  heatmap.legend.labels <- c("min", "max")

  if(min(val.range) < 0){
    val.range <- c(-max(abs(val.range)),max(abs(val.range)))
    heatmap.legend.ticks <- c(val.range[1],0,val.range[2])
    heatmap.legend.labels <- c("min",0,"max")
  }

  heatcolors <- colorRamp2(seq(val.range[1], val.range[2], length = length(col)), col)

  dendc <- NULL
  column_title <- NULL
  column_split <- NULL
  show_column_dend <- FALSE
  cluster_rows <- TRUE
  cluster_column_slices <- cluster.group.slices
  #show.group.bar.legend <- F
  ## We don't input single cell so clustering within is same as group clustering
  if(cellset_ordered){
    message("Using existing cell order")
    dendc <- F # prevents clustering of columns
    #cluster_rows <- F #as.dendrogram(hclust(dist(mat)))
    mat <- mat[,cellset]
  }
  else if( (ncol(mat) >= 65000)){
    message("Complexheatmap cannot handle dendrogram for large matrices, turning off dendrogram.")
    clustering.within <- F
  }
  else if(no.clustering==T){
    message("No clustering of groups")
    cluster_rows <- as.dendrogram(hclust(dist(mat)))
    show.group.bar.legend <- T # Show the legend
  }
  else{
    if(clustering.within){
      message("Performing within & between group clustering")
      #dendc <- cluster_within_group(mat, group_anno) # both groups and within groups are reordered, better for discerning clusters, but looks rougher
      dendc <- cluster_within_group(mat, group_anno)
      # For large heatmaps clustering rows with dendrogram works best
      #cluster_rows <- as.dendrogram(hclust(dist(mat, method = 'binary')))
      cluster_rows <- as.dendrogram(hclust(dist(mat)))
      column_split <- ngroup
      message("Finished clustering")
    }
    else{
      message("Performing between group clustering")
      #dendc <- cluster_between_groups(mat, group_anno) # only reordered between groups, looks more smooth
      dendc <- cluster_between_groups(mat, group_anno) # only reordered between groups, looks more smooth
      # For smaller heatmaps (with labels) clustering rows based on dendc works better
      cluster_rows <- as.dendrogram(hclust(dist(mat)))
      column_split <- ngroup
    }
  }

  if(show.group.dendrogram){
    message("Show group dendrogram")
    if(is.null(dendc)){
      dendc <- T
      column_title <- NULL
    }
    show_column_dend <- T
  } else if (show.cell.dendrogram){
    message("Show cell dendrogram")
    cluster_column_slices <- TRUE
    show_column_dend <- T
    column_split <- cluster_anno
  }

  nclup <- apply(mat, 1, function(x) { length(which(x>0)) } )
  nclempty <- apply(mat, 1, function(x) { length(which(is.na(x))) } )
  ncldown <- apply(mat, 1, function(x) { length(which(x<0)) } )
  fracup <- nclup / (nclup + ncldown + nclempty)
  fracdown <- ncldown / (nclup + ncldown + nclempty)
  clusterRows <- F

  # get families
  (cluster.names <- gsub(".*#","", rownames(mat)) %>% gsub("-[0-9]+-[01].*","",.))

  # clean up motif names
  (rownames(mat) <- gsub("#.*","",rownames(mat)))

  ht <- NULL

  # horizontal
  if(orientation_horizontal){
    message("\nCreating horizontal heatmap")
    # Barplot annotation
    ha <- NULL
    # these are number of events
    if(labelMotifN){
      ha = rowAnnotation(
        dist0 = anno_barplot(
          which = 'row',
          #x = nclup[clusterCols]+ncldown[clusterCols],
          x = nclup+ncldown,
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = 'black'),
          width = unit(2, "cm")
        )
        ,show_annotation_name = FALSE )
    } else if(labelMotifFrac){
      # This is fraction of events
      right_anno = rowAnnotation(
        dist3 = anno_barplot(
          which = 'row',
          #x = cbind(fracup[clusterCols], fracdown[clusterCols]),
          x = cbind(fracup, fracdown),
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = c(tail(col, n=1), head(col, n=1)) ), #up color, down color
          #gp = gpar(col = "white", fill = tail(col, n=1)),
          width = unit(2, "cm")
        ),
        show_annotation_name = FALSE)
    }

    # label clusters
    colAnn <- NULL
    if(!missing(labelMotifClusters)){
      # cluster annotation
      rloc <- grep(paste(labelMotifClusters, collapse = "-"), cluster.names, ignore.case = T)
      label.names <- unique(cluster.names[rloc])
      #labels <- label.col[seq_along(label.names)]
      labels <- label.col[label.names]
      names(labels) <- label.names
      alllabels <- cluster.names
      alllabels[setdiff(seq_along(alllabels), rloc)] <- NA
      if(is.na(legend.ncol))
        legend.ncol <- ifelse(length(labels) < 25, 2, ifelse(length(labels) < 18, 3, 4))
      colAnn <- HeatmapAnnotation(df=data.frame('Clustered'=alllabels),
                                  col = list('Clustered'=labels),
                                  which='col',
                                  show_annotation_name = FALSE,
                                  na_col = na.col,
                                  gp = gpar(col = "black", lwd=label.outline),
                                  # legend for cluster annotation
                                  annotation_legend_param = list(legend_height = unit(5, "in"),
                                                                 labels_gp = gpar(fontsize = 12),
                                                                 border = TRUE,
                                                                 ncol = legend.ncol))

    }
    if(outline)
      outline <- gpar(col = "grey10", lwd = 0.1)
    else outline <- gpar(col = NA)

    ha <- NULL
    if(group.annotation.bar){
      dfgroup <- NULL
      colgroup <- NULL
      if(group.by == 'Condition'){
        dfgroup <- data.frame(group = condition_anno,
                              bionames = bionames_anno,
                              cluster = cluster_anno,
                              celltype = celltype_anno)
        colgroup <- list(group = c(group.colors),
                         bionames = c(bionames_palette_list),
                         cluster = c(cluster_palette),
                         celltype = c(mammary_refine_clusters_palette))
      } else if(group.by == 'bioNames'){
        dfgroup <- data.frame(group = bionames_anno,
                              condition = condition_anno,
                              cluster = cluster_anno,
                              celltype = celltype_anno)
        colgroup <- list(group = c(group.colors),
                         condition = condition_palette,
                         cluster = cluster_palette,
                         celltype = mammary_refine_clusters_palette)
      } else if (group.by == 'Clusters'){
        dfgroup <- data.frame(group = cluster_anno,
                              condition = condition_anno,
                              bionames = bionames_anno,
                              celltype = celltype_anno)
        colgroup <- list(group = c(group.colors),
                         condition = condition_palette,
                         bionames = bionames_palette_list,
                         celltype = mammary_refine_clusters_palette)
      } else if (group.by == 'gex_mammary_gland_refine_clusters'){
        dfgroup <- data.frame(group = celltype_anno,
                              condition = condition_anno,
                              bionames = bionames_anno,
                              cluster = cluster_anno)
                              #pseudotime = pseudotimes)
        colgroup <- list(group = c(group.colors),
                         condition = condition_palette,
                         cluster = cluster_palette,
                         bionames = bionames_palette_list)
                         #pseudotime = pseudotime_palette)
      }
      ha = HeatmapAnnotation(df = dfgroup,
                             which = 'row',
                             col = colgroup,
                             width = unit(x = group.bar.height, units = 'cm'),
                             show_annotation_name = F,
                             show_legend = show.group.bar.legend,
                             border=T,
                             annotation_legend_param = list(legend_side='right',
                                                            labels_gp = gpar(fontsize = 12),
                                                            grid_border = 'black',
                                                            border='black'))
    }

    ht <- ComplexHeatmap::Heatmap(name=title.legend,
                                  matrix = t(mat),
                                  column_title = title,
                                  col = heatcolors,
                                  na_col = na.col,
                                  rect_gp = default.cell.outline, # creates "space" between cells
                                  border = T,
                                  border_gp = gpar(col='black', lwd=.2),
                                  column_gap = unit(heatmap.gap.mm, "mm"),
                                  show_column_dend = F,
                                  cluster_columns = cluster_rows,
                                  row_split = column_split,
                                  cluster_rows = dendc,
                                  use_raster = use.raster,
                                  show_column_names = labelMotifs,
                                  column_names_gp = gpar(fontsize=column_name_fontsize),
                                  column_names_side = 'top',
                                  show_row_names = F,
                                  row_names_side = 'left',
                                  heatmap_legend_param = list(at = heatmap.legend.ticks,
                                                              labels = heatmap.legend.labels,
                                                              title_position = 'topcenter',
                                                              title_gp = gpar(fontsize = 12),
                                                              grid_width = unit(2, "cm"),
                                                              grid_height = unit(.5, "cm"),
                                                              grid_border = 'black',
                                                              border='black',
                                                              gp = gpar(col = "black", lwd=label.outline),
                                                              legend_direction = 'horizontal'),
                                  #right_annotation = ha,
                                  top_annotation=colAnn,
                                  left_annotation = ha
    )

    if(labelMotifN){
      # just text
      cluster_anno = anno_mark(at = seq_len(ncol(mat)),
                               labels_gp = gpar(fontsize = row_name_fontsize),
                               labels = paste("n=",nclup[clusterCols]+ncldown[clusterCols]),
                               which = "row",
                               link_width = unit(3, "mm"))
      ht <- ht + rowAnnotation(mark = cluster_anno)
    }
    #ht <- draw(ht, heatmap_legend_side = "bottom", padding = unit(c(3, 2, 5, 8), "mm")) #bltr

    #### vertical orientation ####
  } else {
    message("\nCreating vertical heatmap of ",nrow(mat)," features and ",ncol(mat), " cells")
    # Barplot annotation
    right_anno <- NULL
    # these are number of events
    if(labelMotifN){
      right_anno = rowAnnotation(
        dist0 = anno_barplot(
          which = 'row',
          #x = nclup[clusterCols]+ncldown[clusterCols],
          x = nclup+ncldown,
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = 'black'),
          width = unit(2, "cm")
        )
        ,show_annotation_name = FALSE )
    } else if(labelMotifFrac){
      # This is fraction of events
      right_anno = rowAnnotation(
        dist3 = anno_barplot(
          which = 'row',
          x = cbind(fracup, fracdown),
          bar_width = 1,
          border = FALSE,
          axis_param = list("side"="top"),
          gp = gpar(col = "white", fill = c(tail(col, n=1), head(col, n=1)) ), #up color, down color
          numbers_rot = 0,
          width = unit(2, "cm")
        ),
        show_annotation_name = FALSE)
    }

    # label motif family
    colAnn <- NULL
    if(!missing(labelMotifClusters)){
      # cluster annotation
      rloc <- grep(paste(labelMotifClusters, collapse = "-"), cluster.names, ignore.case = T)
      label.names <- unique(cluster.names[rloc])
      #labels <- label.col[seq_along(label.names)]
      labels <- label.col[label.names]
      names(labels) <- label.names
      alllabels <- cluster.names
      alllabels[setdiff(seq_along(alllabels), rloc)] <- NA
      legend.nrow <- ifelse(length(labels) >= 40, 30, ifelse(length(labels) <= 25, 25, 20))
      #legend.nrow <- ifelse(length(labels) <= 30, 30, ifelse(length(labels) <= 25, 25, 20))
      colAnn <- HeatmapAnnotation(df=data.frame('Clustered'=alllabels),
                                  col = list('Clustered'=labels),
                                  which='row',
                                  show_annotation_name = FALSE,
                                  na_col = na.col,
                                  gp = gpar(col = "black", lwd=label.outline),
                                  # legend for cluster annotation
                                  annotation_legend_param = list(legend_height = unit(4, "in"),
                                                                 labels_gp = gpar(fontsize = 12),
                                                                 border = TRUE,
                                                                 nrow = legend.nrow))

    }
    if(outline)
      outline <- gpar(col = "grey10", lwd = 0.1)
    else outline <- gpar(col = NA)

    cat("group_anno: ", levels(group_anno))

    ha <- NULL
    if(group.annotation.bar){
      dfgroup <- NULL
      colgroup <- NULL
      if(group.by == 'Condition'){
        dfgroup <- data.frame(group = group_anno,
                              bionames = bionames_anno,
                              cluster = cluster_anno,
                              celltype = celltype_anno)
        colgroup <- list(group = c(group.colors),
                         bionames = c(bionames_palette_list),
                         cluster = c(cluster_palette),
                         celltype = c(mammary_refine_clusters_palette))
      } else if(group.by == 'bioNames'){
        dfgroup <- data.frame(group = group_anno,
                              condition = condition_anno,
                              cluster = cluster_anno,
                              celltype = celltype_anno)
        colgroup <- list(group = c(group.colors),
                         condition = condition_palette,
                         cluster = cluster_palette,
                         celltype = mammary_refine_clusters_palette)
      } else if (group.by == 'Clusters'){
        dfgroup <- data.frame(group = group_anno,
                              condition = condition_anno,
                              bionames = bionames_anno,
                              celltype = celltype_anno)
        colgroup <- list(group = group.colors,
                         condition = condition_palette,
                         bionames = bionames_palette_list,
                         celltype = mammary_refine_clusters_palette)
      } else if (group.by == 'gex_mammary_gland_refine_clusters'){
        dfgroup <- data.frame(group = group_anno)
                              #condition = condition_anno,
                              #bionames = bionames_anno)
                              #cluster = cluster_anno)
        colgroup <- list(group = c(group.colors))
                         #condition = condition_palette,
                         #bionames = bionames_palette_list)
                         #cluster = cluster_palette)
        if(!is.null(pseudotimes)){
          dfgroup <- cbind(dfgroup, data.frame(pseudotime=pseudotimes))
          colgroup[["pseudotime"]] <- pseudotime_palette
        }
      } else if (group.by == 'gex_mammary_gland_cluster_classification'){
        dfgroup <- data.frame(group = group_anno,
                              condition = condition_anno,
                              bionames = bionames_anno,
                              cluster = cluster_anno)
        colgroup <- list(group = c(group.colors),
                         condition = condition_palette,
                         cluster = cluster_palette,
                         bionames = bionames_palette_list)
      }

      ha = HeatmapAnnotation(df = dfgroup,
                             col = colgroup,
                             height = unit(x = group.bar.height, units = 'cm'),
                             show_annotation_name = F,
                             show_legend = show.group.bar.legend,
                             border=T, annotation_legend_param = list(
                               labels_gp = gpar(fontsize = 12),
                               grid_border = 'black',
                               border='black'
                             ))
    }

    ht <- ComplexHeatmap::Heatmap(name=title.legend,
                                  matrix = mat,
                                  column_title = title,
                                  col = heatcolors,
                                  na_col = na.col,
                                  border = T,
                                  border_gp = gpar(col='black', lwd=.2),
                                  rect_gp = default.cell.outline, # creates "space" between cells
                                  show_column_dend = show.group.dendrogram,
                                  cluster_columns = dendc,
                                  cluster_rows = cluster_rows,
                                  column_split = column_split,
                                  row_gap = unit(heatmap.gap.mm, "mm"),
                                  show_row_dend = F,
                                  use_raster = use.raster,
                                  show_row_names = labelMotifs,
                                  show_column_names = F,
                                  row_names_gp = gpar(fontsize = row_name_fontsize),
                                  column_names_side = 'top',
                                  column_names_centered = F,
                                  row_names_side = 'left',
                                  heatmap_legend_param = list(at = heatmap.legend.ticks,
                                                              labels = heatmap.legend.labels,
                                                              title_position = 'topleft',
                                                              #title_position = "leftcenter-rot",
                                                              title_gp = gpar(col = "black", lwd=label.outline),
                                                              #title_gp = gpar(fontsize = 12),
                                                              labels_gp = gpar(fontsize = 12),
                                                              grid_width = unit(.5, "cm"),
                                                              grid_height = unit(1.5, "cm"),
                                                              grid_border = 'black',
                                                              border='black'),
                                  top_annotation = ha,
                                  right_annotation = right_anno,
                                  left_annotation = colAnn
    )
    if(!missing(labelSpecificMotifs)){
      rloc <- grep(paste(labelSpecificMotifs, collapse = "-"), rownames(mat), ignore.case = T)
      labels <- rownames(mat)[rloc]
      motif_anno = anno_mark(at = rloc,
                             labels_gp = gpar(fontsize = 12),
                             labels = labels,
                             which = "row",
                             link_width = unit(8, "mm"))
      ht <- ht + rowAnnotation(mark = motif_anno)
    }
  }

  ht
}

#' Heatmap of motif enrichment across groups
#'
#' @param enrichMotifs Motif enrichment result (SummarizedExperiment).
#' @param n Number of top motifs per group.
#' @param mlog10Padj.cutOff Minimum -log10 adjusted p-value.
plotEnrichMotifHeatmap <- function(enrichMotifs,
                                   geneMotifSameMat=NULL,
                                   group.by=NULL,
                                   n=10,
                                   mlog10Padj.cutOff=0,
                                   labelMotifs=T,
                                   labelMotifN=T,
                                   labelMotifFrac=F,
                                   labelSpecificMotifs,
                                   labelMotifClusters,
                                   remove.redundant=T,
                                   cluster.columns=T,
                                   cluster.rows=F,
                                   scaling = 'None', #column
                                   add.score.to.name=T,
                                   col=unname(ArchRPalettes$whitePurple),
                                   na.col='gray90',
                                   cell.outline=T,
                                   label.col=pals::polychrome(36),
                                   orientation_horizontal=T,
                                   title.legend ="Average Accessibility",
                                   title = "Enriched transcription factor motifs per family",
                                   row_title = "",
                                   row_name_fontsize=12,
                                   column_name_fontsize=7,
                                   outline=T,
                                   quantile.min=.01,
                                   quantile.max=.99,
                                   label.outline=0.5){

  require(circlize)
  require(ComplexHeatmap)

  default.cell.outline <- gpar(col = "gray40", lwd = 2)

  if(!cell.outline)
    default.cell.outline <- gpar(col = NA)

  mat <- NULL
  if ('matrix' %in% class(enrichMotifs))
    mat <- enrichMotifs
  else if('SummarizedExperiment' %in% class(enrichMotifs))
    mat <- as.matrix(assays(enrichMotifs)[["mlog10Padj"]])
  else if ('data.table' %in% class(enrichMotifs)){
    mat <- as.matrix(enrichMotifs, rownames='motifID')
    mat <- mat[,grep("mlog10Padj", colnames(mat))]
    colnames(mat) <- gsub('mlog10Padj','',colnames(mat))
  }

  # remove empty
  ridx <- which(colSums(mat)==0)
  if(length(ridx)>0)
    mat <- mat[,-ridx]

  # remove.redundant
  # Keeps n items with padj value > mlog10Padj.cutOff
  keep <- lapply(seq_len(ncol(mat)), function(x) {
    sorted <- mat[order(mat[, x], decreasing = T),]
    mnames <- gsub("#.*","",rownames(sorted))
    if(remove.redundant){
      dups <- which(duplicated(mnames))
      if(length(dups) > 0){
        sorted <- sorted[-dups,]
      }
    }
    pass <- which(sorted[, abs(x)] >= abs(mlog10Padj.cutOff))
    head(names(pass), n)
  }) %>% unlist %>% unique
  if(length(keep)==0)
    stop("All entries filtered out, consider reducing cutoff")
  mat <- mat[keep, , drop = FALSE]

  mat[is.infinite(mat)] <- NA

  # subsetting
  if(cluster.columns){
    bS <- ArchR:::.binarySort(m = mat, clusterCols = cluster.columns)
    mat <- bS[[1]][, colnames(mat)]
    clusterCols <- 1:ncol(mat)
    clusterCols <- bS[[2]][['order']]
    mat <- mat[, clusterCols, drop = FALSE]
  } else{ # order by score
    cnames <- melt(data.table(mat, keep.rownames = 'motifid'), id.vars = 'motifid')[order(value, decreasing = T), motifid] %>% unique()
    mat <- mat[cnames,]
  }

  if(!is.null(geneMotifSameMat)){
    geneMotifSameMat <- geneMotifSameMat[rownames(mat), colnames(mat)]
  } else {
    geneMotifSameMat <- matrix(data = NA, nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
  }

  # get families

  #cluster.names <- gsub(".*#","", rownames(mat)) %>% gsub("-[0-9]+-[01].*","",.)
  cluster.names <- sapply(strsplit(x = rownames(mat), split='-'), function(x){ x[[ (length(x)-2) ]] })
  #cat(cluster.names)

  # clean up motif names
  rownames(mat) <- gsub("#.*","",rownames(mat))

  # clean up row names
  if(add.score.to.name){
    rownames(mat) <- round(apply(mat, 1, max)) %>% paste0(names(.)," (",.,")")
  }

  rowScale <- function(mat=NULL, min=NULL, max=NULL){
    if (!is.null(min)) {
      rMin <- min
    }
    else {
      rMin <- matrixStats::rowMins(mat, na.rm = T)
    }
    if (!is.null(max)) {
      rMax <- max
    }
    else {
      rMax <- matrixStats::rowMaxs(mat, na.rm = T)
    }
    cat("length rMin",length(rMin),"\n")
    rScale <- pmax(rMax, abs(rMin)) # centering point
    cat("length rScale",length(rScale),"\n")
    matDiff <- mat # - rMin
    matScale <- matDiff/rScale # centering
    out <- list(mat = matScale, min = rMin, max = rMax)
    return(out)
  }

  # scaling
  if(scaling=='sample'){
    #mat <- ArchR:::.rowScale(t(mat), min = 0)[[1]] %>% t(.) # scales across clusters
    mat <- rowScale(t(mat))[[1]] %>% t(.) # scales across clusters
    scaling <- "\nscaled by sample"
  }
  else if(scaling=='motif'){
    ##mat <- ArchR:::.rowScale(mat, min = 0)[[1]] # scales TF across motifs
    mat <- rowScale(mat)[[1]] # scales TF across motifs
    scaling <- "\nscaled by motif"
  }
  else{
    message("No scaling")
    scaling <- ""
  }
  val.range <- c(quantile(mat, quantile.min, na.rm=T), quantile(mat, quantile.max, na.rm=T))
  #val.range <- range(mat, na.rm = T)
  cat("scale.range = ",val.range,"\n")
  heatmap.legend.ticks <- c(val.range[1],val.range[2])
  heatmap.legend.labels <- c("min", "max")

  if(min(val.range, na.rm = T) < 0){
    val.range <- c(-max(abs(val.range), na.rm = T),max(abs(val.range), na.rm = T))
    heatmap.legend.ticks <- c(val.range[1],0,val.range[2])
    heatmap.legend.labels <- c("min",0,"max")
  }

  cat("new scale.range = ",val.range,"\n")

  heatcolors <- colorRamp2(seq(val.range[1], val.range[2], length = length(col)), col)

  #ncl <- apply(mat, 2, function(x) { length(which(!is.na(x))) } )
  nclup <- apply(mat, 2, function(x) { length(which(x>0)) } )
  nclempty <- apply(mat, 2, function(x) { length(which(is.na(x))) } )
  ncldown <- apply(mat, 2, function(x) { length(which(x<0)) } )
  fracup <- nclup / (nclup + ncldown + nclempty)
  fracdown <- ncldown / (nclup + ncldown + nclempty)
  cat(paste(colnames(mat),nclup+ncldown,collapse = ", "),"\n")
  clusterRows <- F
  if(cluster.rows){
    cmat <- mat
    cmat[is.na(cmat)] <- 0
    cmat[is.nan(cmat)] <- 0
    cmat[is.infinite(cmat)] <- 0
    clusterRows <- as.dendrogram(hclust(dist(cmat)))
  }

  clusterCols <- 1:ncol(mat)

  cluster_anno <- clusterCols

  if(!is.null(group.by))
    cluster_anno<- factor(proj@cellColData[colnames(mat),group.by])
  cat(levels(cluster_anno))

  # get families
  #cluster.names <- gsub(".*#","", rownames(mat)) %>% gsub("-[0-9]+-[01].*","",.)
  #cat(cluster.names)

  # clean up motif names
  #rownames(mat) <- gsub("#.*","",rownames(mat))

  ht <- NULL

  # horizontal
  if(orientation_horizontal){
    message("\nCreating horizontal heatmap")
    # Barplot annotation
    ha <- NULL
    # these are number of events
    if(labelMotifN){
      ha = rowAnnotation(
        dist0 = anno_barplot(
          which = 'row',
          x = nclup[clusterCols]+ncldown[clusterCols],
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = 'black'),
          width = unit(2, "cm")
        )
        # dist1 = anno_barplot(
        #   which = 'row',
        #   x = nclup[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = tail(col, n=1)),
        #   width = unit(2, "cm")
        # ),
        # dist2 = anno_barplot(
        #   which = 'row',
        #   x = ncldown[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = col[1]),
        #   width = unit(2, "cm")
        # )
        ,show_annotation_name = FALSE )
    } else if(labelMotifFrac){
      # This is fraction of events
      ha = rowAnnotation(
        dist3 = anno_barplot(
          which = 'row',
          x = cbind(fracup[clusterCols], fracdown[clusterCols]),
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = c(tail(col, n=1), head(col, n=1)) ), #up color, down color
          #gp = gpar(col = "white", fill = tail(col, n=1)),
          width = unit(2, "cm")
        ),
        show_annotation_name = FALSE)
    }

    # label clusters
    colAnn <- NULL
    if(!missing(labelMotifClusters)){
      # cluster annotation
      rloc <- grep(paste(labelMotifClusters, collapse = "-"), cluster.names, ignore.case = T)
      label.names <- unique(cluster.names[rloc])
      labels <- label.col[seq_along(label.names)]
      names(labels) <- label.names
      alllabels <- cluster.names
      alllabels[setdiff(seq_along(alllabels), rloc)] <- NA
      nc <- ifelse(length(labels) < 20, 2, 3)
      colAnn <- HeatmapAnnotation(df=data.frame('Clustered'=alllabels),
                                  col = list('Clustered'=labels),
                                  which='col',
                                  show_annotation_name = FALSE,
                                  na_col = na.col,
                                  gp = gpar(col = "black", lwd=label.outline),
                                  # legend for cluster annotation
                                  annotation_legend_param = list(legend_height = unit(5, "in"),
                                                                 labels_gp = gpar(fontsize = 12),
                                                                 border = TRUE,
                                                                 ncol = nc))

    }
    if(outline)
      outline <- gpar(col = "grey10", lwd = 0.1)
    else outline <- gpar(col = NA)

    ht <- ComplexHeatmap::Heatmap(name=paste0(title.legend, scaling),
                                  matrix = t(mat),
                                  column_title = title,
                                  #col = unname(ArchRPalettes$comet),
                                  col = heatcolors,
                                  na_col = na.col,
                                  rect_gp = default.cell.outline, # creates "space" between cells
                                  border = T,
                                  border_gp = gpar(col='black', lwd=.2),
                                  show_column_dend = F,
                                  cluster_columns = F,
                                  cluster_rows = F,
                                  use_raster = F,
                                  show_column_names = labelMotifs,
                                  column_names_gp = gpar(fontsize=column_name_fontsize),
                                  column_names_side = 'top',
                                  row_names_side = 'left',
                                  heatmap_legend_param = list(at = heatmap.legend.ticks,
                                                              labels = heatmap.legend.labels,
                                                              title_position = 'topcenter',
                                                              title_gp = gpar(fontsize = 12),
                                                              grid_width = unit(2, "cm"),
                                                              grid_height = unit(.5, "cm"),
                                                              grid_border = 'black',
                                                              border='black',
                                                              gp = gpar(col = "black", lwd=label.outline),
                                                              legend_position = 'right',
                                                              legend_direction = 'horizontal'),
                                  right_annotation = ha,
                                  top_annotation=colAnn,
                                  ## add text to each grid
                                  layer_fun = function(j, i, x, y, w, h, col) {
                                    #grid.points(x, y, pch=ifelse(geneMotifSameMat[j, i]=="","","*"), size = unit(4, "mm"), gp=gpar(col='white'))
                                    if(cell.outline)
                                      grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "white", lwd = .5)) # outlines cells
                                  }
    )

    if(labelMotifN){
      # just text
      cluster_anno = anno_mark(at = seq_len(ncol(mat)),
                               labels_gp = gpar(fontsize = row_name_fontsize),
                               labels = paste("n=",nclup[clusterCols]+ncldown[clusterCols]),
                               which = "row",
                               link_width = unit(3, "mm"))
      ht <- ht + rowAnnotation(mark = cluster_anno)
    }
    ht <- draw(ht, heatmap_legend_side = "bottom", padding = unit(c(3, 2, 5, 8), "mm")) #bltr

    # vertical orientation
  } else {
    message("\nCreating vertical heatmap")
    # Barplot indicating number of motifs top annotation
    ha <- NULL
    if(labelMotifN){
      ha = columnAnnotation(
        dist0 = anno_barplot(
          which = 'column',
          x = nclup[clusterCols]+ncldown[clusterCols],
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = 'black'),
          width = unit(2, "cm")
        )
        # dist1 = anno_barplot(
        #   which = 'row',
        #   x = nclup[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = tail(col, n=1)),
        #   width = unit(2, "cm")
        # ),
        # dist2 = anno_barplot(
        #   which = 'row',
        #   x = ncldown[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = col[1]),
        #   width = unit(2, "cm")
        # )
        ,show_annotation_name = FALSE )
    } else if(labelMotifFrac){
      # This is fraction of events
      ha = columnAnnotation(
        dist3 = anno_barplot(
          which = 'column',
          x = cbind(fracup[clusterCols], fracdown[clusterCols]), #up, down
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = c(tail(col, n=1), head(col, n=1)) ), #up color, down color
          width = unit(2, "cm")
        ),
        show_annotation_name = FALSE)
    }

    ht <- ComplexHeatmap::Heatmap(name=paste0(title.legend, scaling),
                                  matrix = mat,
                                  column_title = title,
                                  row_title = row_title,
                                  top_annotation = ha,
                                  col = heatcolors,
                                  na_col = na.col,
                                  border = T,
                                  border_gp = gpar(col='black', lwd=.2),
                                  rect_gp = default.cell.outline, # creates "space" between cells
                                  show_column_dend = T,
                                  cluster_columns = F,
                                  cluster_rows = clusterRows,
                                  show_row_dend = F,
                                  use_raster = F,
                                  show_row_names = labelMotifs,
                                  row_names_gp = gpar(fontsize = row_name_fontsize),
                                  column_names_side = 'top',
                                  column_names_centered = F,
                                  row_names_side = 'left',
                                  heatmap_legend_param = list(at = heatmap.legend.ticks,
                                                              labels = heatmap.legend.labels,
                                                              title_position = 'topleft',
                                                              #title_position = "leftcenter-rot",
                                                              title_gp = gpar(fontsize = 12),
                                                              labels_gp = gpar(fontsize = 12),
                                                              grid_width = unit(.5, "cm"),
                                                              grid_height = unit(1.5, "cm"),
                                                              grid_border = 'black',
                                                              border='black',
                                                              legend_side='bottom',
                                                              legend_direction = 'vertical'),
                                  # add text to each grid
                                  # add text to each grid
                                  layer_fun = function(j, i, x, y, w, h, col) {
                                    #grid.points(x, y, pch=ifelse(geneMotifSameMat[i, j]=="","","*"), size = unit(4, "mm"), gp=gpar(col='white'))
                                    if(cell.outline)
                                      grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "white", lwd = .5)) # outlines cells
                                  }
    )

    # Label
    if(!missing(labelSpecificMotifs)){
      rloc <- grep(paste(labelSpecificMotifs, collapse = "-"), rownames(mat), ignore.case = T)
      labels <- rownames(mat)[rloc]
      motif_anno = anno_mark(at = rloc,
                             labels_gp = gpar(fontsize = 12),
                             labels = labels,
                             which = "row",
                             link_width = unit(8, "mm"))
      ht <- ht + rowAnnotation(mark = motif_anno)
    }

    # Cluster colored tiles
    if(!missing(labelMotifClusters)){
      rloc <- grep(paste(labelMotifClusters, collapse = "-"), cluster.names, ignore.case = T)
      label.names <- unique(cluster.names[rloc])
      labels <- label.col[seq_along(label.names)]
      names(labels) <- label.names
      alllabels <- cluster.names
      alllabels[setdiff(seq_along(alllabels), rloc)] <- NA
      ht <- ht + HeatmapAnnotation(df=data.frame('Clustered'=alllabels),
                                   col = list('Clustered'=labels),
                                   which='row',
                                   show_annotation_name = FALSE,
                                   na_col = na.col,
                                   gp = gpar(col = "black", lwd=label.outline),
                                   annotation_legend_param = list(legend_height = unit(5, "in"),
                                                                  labels_gp = gpar(fontsize = 13)))
    }
    ht <- draw(ht, heatmap_legend_side = "right", padding = unit(c(3, 2, 5, 8), "mm")) #bltr
  }

  ht
}

#' Bar plot of peak genomic annotations (promoter, intron, ...)
PlotPeakAnnotations <- function(ArchRProj){
  lookup <- unique(data.table(data.frame(ArchRProj@cellColData[,c("Clusters","Clusters2")])))

  peakDT <- data.table(metadata(ArchRProj@peakSet)$PeakCallSummary)
  peakDT[,cluster := gsub("\\(.*","",Group)]
  peakDT$cluster_name <- lookup[match(x = peakDT$cluster, table = lookup$Clusters), Clusters2]
  peakDT[is.na(cluster_name), cluster_name := ""]
  peakDT[,cluster_name := paste(cluster_name, gsub(".*\\(","(",Group))]
  peakDT[,clust := as.integer(gsub("^C","",gsub("-.*","",cluster)))]
  peakDT <- peakDT[order(clust),]
  peakDT[,cluster_name := factor(x = cluster_name, levels = unique(cluster_name), ordered = T)]

  pal <- paletteDiscrete(values = peakDT$Var1)
  pal <- c(Distal = "#60BA64", Exonic = "#73C6FF", Intronic = "#620FA3", Promoter = "#FFC554", pal)

  p1 <- ggplot(peakDT, aes(x=cluster_name, y=Freq, fill=Var1)) +
    ylab(expression("Number of Peaks "~(x10^3))) + xlab("") +
    geom_col(col='black') +
    scale_fill_manual("Annotation",values = pal) +
    theme_ArchR(baseSize = 12) + theme(text = element_text(family = 'Arial'),
                                       panel.border = element_blank(),
                                       axis.line = element_line(linewidth = .4, colour = 'black'),
                                       plot.background = element_blank(),
                                       plot.margin = margin(.5, 0, .5, 2, "cm"),
                                       axis.text.x = element_text(angle = 45, hjust=1),
                                       legend.position = 'right',
                                       legend.text = element_text(size=11),
                                       legend.box.background = element_blank(),
                                       legend.background = element_blank())

  cairo_pdf(filename = file.path(getOutputDirectory(ArchRProj),'Plots','Peak_annotation_region.pdf'),
            width = 7, height = 5,
            bg = 'transparent')
  print(p1)
  dev.off()
}

#' Bar plot of up/down peaks in pairwise cluster comparisons
PlotBarDEPeaks <- function(ArchRProj, seDETests, absLog2FC=1, fdr=0.05, width=15, height=6, outfile.name="Peaks_in_pairwise_de_clusters_Barplot.pdf"){
  dtSummaryPairPlots <- lapply(names(seDETests), function(i){
    dtsummary <- data.table(FDR=seDETests[[i]]@assays@data$FDR$x,
                            Log2FC=seDETests[[i]]@assays@data$Log2FC$x)
    dtsummary <- dtsummary[FDR <= fdr & abs(Log2FC) >= absLog2FC, list(direction=sign(Log2FC))][, .N, by=list(direction)]
    dtsummary[,name:=i]
    dtsummary[,numerator:=gsub("_vs_.*","",name)]
    dtsummary[,denominator:=gsub(".*_vs_","",name)]
  }) %>% rbindlist()

  dtSummaryPairPlotsRev <- copy(dtSummaryPairPlots)
  dtSummaryPairPlotsRev[,direction := direction*-1]
  dtSummaryPairPlotsRev[,temp:=numerator]
  dtSummaryPairPlotsRev[,numerator:=denominator]
  dtSummaryPairPlotsRev[,denominator:=temp]
  dtSummaryPairPlotsRev[,temp:=NULL]

  dtSummaryPair <- rbindlist(list(dtSummaryPairPlots, dtSummaryPairPlotsRev))
  dtSummaryPair <- dtSummaryPair[order(numerator, denominator)]
  dtSummaryPair[,clust := as.integer(gsub("^C","",gsub("-.*","",denominator)))]
  dtSummaryPair <- dtSummaryPair[order(clust, direction)]
  dtSummaryPair[,denominator := factor(x = denominator, levels = unique(denominator), ordered = T)]

  dtSummaryPair[,clust := as.integer(gsub("^C","",gsub("-.*","",numerator)))]
  dtSummaryPair <- dtSummaryPair[order(clust, direction)]
  dtSummaryPair[,numerator := factor(x = numerator, levels = unique(numerator), ordered = T)]

  ## Barplot
  ylim <- 110000 #round(max(abs(dtSummaryPair$N)))
  g <- ggplot(dtSummaryPair, aes(x=denominator, y=direction)) +
    geom_col(data = dtSummaryPair[direction > 0], aes(y=direction*N), fill='indianred') +
    geom_col(data = dtSummaryPair[direction < 0], aes(y=direction*N), fill='dodgerblue') +
    ggtitle("Significantly Differential Peak Accessibility", subtitle = bquote(FDR<= ~ .(fdr)~","~Log2FC>= ~ .(absLog2FC))) +
    theme_ArchR() +
    geom_hline(yintercept = 0, linetype='solid') +
    scale_y_continuous(breaks = seq(from=-ylim, to=ylim, by=10000), labels = abs(seq(from=-ylim/1000, to=ylim/1000, by=10))) +
    ylab(expression("Number of Peaks (x"~10^3~")")) + xlab("") +
    #facet_wrap( ~group, nrow = 1, scales = "free_x") +
    facet_grid(~numerator, scales = "free_x", space = 'free_x') +
    theme(strip.placement = "inside",
          plot.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, angle = 90),
          axis.text.x = element_text(angle=90, hjust = 1, vjust=.5),
          panel.spacing.y = unit(.1, "lines"),
          panel.spacing.x = unit(.1, "lines"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
  #g <- g + ggforce::facet_row(vars(numerator), scales = 'free_x', space = 'free', strip.position = 'right')
  pdf(file = file.path(getOutputDirectory(ArchRProj),'Plots',outfile.name), width = width, height = height, bg = 'transparent')
  print(g)
  dev.off()
  dtSummaryPair
}

#' Density of chromVAR deviation z-scores for selected motifs
PlotDeviationZScoreDensity <- function(motifmat, motifs){
  z <- assays(motifmat)
  cl <- colData(motifmat)$Clusters2
  ucl <- names(rna_group_palette) #colData(motifmat)$Clusters2 %>% unique %>% sort
  p1 <- lapply(motifs, function(i){
    ## get motif row
    motif <- grep(i, rownames(z$z), value = T)[1]
    zm <- z$z[motif,]
    dens <- lapply(ucl, function(j){
      zidx <- zm[which(cl==j)]
      zd <- density(zidx)
      data.table(x=zd$x, y=zd$y, cluster=j)
    }) %>% rbindlist()
    p <- ggplot(dens, aes(x=x, y=y, color=cluster)) +
      geom_line(size=.5) +
      xlab("deviation z-score") + ylab("Density") + ggtitle(label = motif) +
      scale_color_manual("",values = rna_group_palette) +
      theme_ArchR(baseSize = 10,
                  legendTextSize = 10,
                  legendPosition = 'right') +
      theme(plot.margin = unit(x = c(0,0,0,0.1), units = 'cm'),
            plot.background = element_blank(),
            legend.background = element_blank())
    p + lemon::coord_capped_cart(expand = F)
  }) %>% set_names(motifs)
  p1
}

#' Bar plot of the top motif deviations
PlotDeviations <- function(dt, n=30, text.adjust=T, title){

  dt <- setDT(dt)
  dt <- dt[order(rank, decreasing = F)]
  dt[,label := paste(gsub("#.*","", motifID), gsub("-.*","",gsub('.*#','',motifID)), sep = ':')]

  if(missing(title))
    title <- paste("Top",n,"deviations")

  g <- ggplot(dt, aes(x = rank, y = combinedVars, color = combinedVars)) +
    geom_point(size = 1) + scale_color_gradientn(colors = paletteContinuous(set = "comet"))

  if(text.adjust)
    g <- g + ggrepel::geom_label_repel(data = dt[seq_len(n) ],
                               aes(x = rank, y = combinedVars, label = label), max.overlaps = 100,
                               size = 3, box.padding = 0.1, label.padding = 0.1, color = "black",
                               segment.size = .1,
                               nudge_x = 1000, direction = "y", hjust = "left")
  else
    g <- g + ggrepel::geom_label_repel(data = dt[seq_len(n)],
                                       aes(x = rank, y = combinedVars, label = label), max.overlaps = 100,
                                       nudge_x = 1000, min.segment.length = 2, direction = "both",
                                       size = 3, box.padding = 0.1, label.padding = 0.1, color = "black", segment.size = .1)

  g <- g + ylab("Variability") + xlab("Rank Sorted Annotations") +
    guides(color = guide_colourbar(title = "Scaled score",direction='horizontal', barwidth = 10, barheight = 1)) +
    ggtitle(title) +
    theme_ArchR(plotMarginCm = .1)
}

#' Rank plot of the top correlations
PlotRankCorrelations <- function(dt, n=30, title="", text.adjust=T){

  dt <- dt[order(rank, decreasing = F)]
  dt[,label:=paste(gsub("#.*","", motifID), family, sep = ':')]

  # we only want to label the top n so clear labels from others
  # need geom_repel to label subset but be aware of all other points
  dt[-seq_len(n), label:=""]

  g <- ggplot(dt, aes(x = rank, y = mlog10Padj, color = rho)) +
    #geom_point(aes(size=-log10(padj)/10), show.legend = F) +
    geom_point(size=3, show.legend = T) +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
    scale_y_continuous(expand = expansion(mult = c(0.1,0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.1,3.5)))

  if(text.adjust)
    g <- g + ggrepel::geom_label_repel(aes(label = label),
                                       max.overlaps = 100,
                                       size = 3, box.padding = 0.1, label.padding = 0.1, color = "black",
                                       segment.size = .1,
                                       nudge_x = 300, direction = "y", hjust = "right")
  else
    g <- g + ggrepel::geom_label_repel(aes(label = label),
                                       max.overlaps = 100,
                                       nudge_x = 200, min.segment.length = 2, direction = "both",
                                       size = 3, box.padding = 0.1, label.padding = 0.1, color = "black", segment.size = .1)

  if(title=="")
    title <- paste("Top",n,"correlated deviations")

  g <- g + ylab("-log10(Padjusted) accessibililty") + xlab("Rank Sorted Annotations") +
    guides(color = guide_colourbar(title = "Correlation",direction='horizontal', barwidth = 10, barheight = 1)) +
    ggtitle(title) +
    theme_ArchR(plotMarginCm = .1)
}

#' Arrange a list of ggplots with one shared legend
#'
#' @param plots List of ggplot objects (legend taken from the first one).
#' @param ncol,nrow Layout of the grid.
#' @param legend.position One of "bottom", "right" or "topright".
grid_arrange_list_shared_legend <- function(plots, ncol = length(plots),
                                            nrow = 1,
                                            legend.position = c("bottom", "right","topright")) {

  # plots <- list(...)
  position <- match.arg(legend.position)
  position.justification <- "center"
  if(position=="topright"){
    position <- "right"
    position.justification <- c(1,.85)
  }
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position, legend.justification = position.justification))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))

  if(ncol < length(plots))
    nrow <- ceiling(length(plots) / ncol)
  else if (nrow > 1)
    ncol <- ceiling(length(plots) / nrow)

  gl <- c(gl, ncol = ncol, nrow = nrow)

  library(grid)
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight),
                                            padding = unit.c(unit(0,"lines"))),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(.95, "npc") - lwidth, lwidth),
                                           #widths = unit.c(unit(.95, "npc"), unit(.05, "npc"), lwidth),
                                           padding = unit.c(unit(0,"lines"))))

  #grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}

#' Patched copy of ArchR::addDeviationsMatrix
#'
#' Local modification of the ArchR function (MIT license, GreenleafLab/ArchR).
addDeviationsMatrix <- function (ArchRProj = NULL, peakAnnotation = NULL, matches = NULL,
                                 bgdPeaks = getBgdPeaks(ArchRProj, method = "chromVAR"), matrixName = NULL,
                                 out = c("z", "deviations"), binarize = FALSE, threads = getArchRThreads(),
                                 verbose = TRUE, parallelParam = NULL, force = FALSE, logFile = createLogFile("addDeviationsMatrix"))
{
  message("Using my own implementation to fix following error:\n")
  message('as(<lgCMatrix>, "dgCMatrix") is deprecated since Matrix 1.5-0; do as(., "dMatrix") instead')

  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = peakAnnotation, name = "peakAnnotation",
                      valid = c("character", "null"))
  ArchR:::.validInput(input = matches, name = "matches", valid = c("SummarizedExperiment",
                                                                   "null"))
  ArchR:::.validInput(input = bgdPeaks, name = "bgdPeaks", valid = c("SummarizedExperiment"))
  ArchR:::.validInput(input = matrixName, name = "matrixName", valid = c("character",
                                                                         "null"))
  ArchR:::.validInput(input = out, name = "out", valid = c("character"))
  ArchR:::.validInput(input = binarize, name = "binarize", valid = c("boolean"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = parallelParam, name = "parallelParam",
                      valid = c("parallelparam", "null"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  if (!inherits(ArchRProj, "ArchRProject")) {
    stop("Error Needs to be ArchR Project for Input!")
  }
  ArrowFiles <- getSampleColData(ArchRProj)$ArrowFiles
  threads <- min(length(ArrowFiles), threads)
  allCells <- rownames(getCellColData(ArchRProj))
  outDir <- getOutputDirectory(ArchRProj)
  if (!all(file.exists(ArrowFiles))) {
    stop("Error Input Arrow Files do not all exist!")
  }
  print(matches)
  if (is.null(matches)) {
    anno <- getPeakAnnotation(ArchRProj, peakAnnotation)
    matches <- readRDS(anno$Matches)
    if (is.null(matrixName)) {
      matrixName <- paste0(anno$Name, "Matrix")
    }
  }
  else {
    if (is.null(matrixName)) {
      matrixName <- paste0("MotifMatrix")
    }
  }
  annotationsMatrix <- SummarizedExperiment::assay(matches)
  rownames(annotationsMatrix) <- paste0(seqnames(matches),
                                        "_", start(matches), "_", end(matches))
  #annotationsMatrix <- as(annotationsMatrix, "dgCMatrix")
  annotationsMatrix <- as(annotationsMatrix, "dMatrix")
  gc()
  ArrowFiles <- getArrowFiles(ArchRProj)
  useMatrix <- "PeakMatrix"
  availableChr <- ArchR:::.availableSeqnames(ArrowFiles, useMatrix)
  rS <- suppressMessages(ArchR:::.getRowSums(ArrowFiles = ArrowFiles,
                                             seqnames = availableChr, useMatrix = useMatrix, filter0 = FALSE))
  rownames(rS) <- paste0(rS$seqnames, "_", rS$idx)
  rS <- rS[paste0(seqnames(ArchRProj@peakSet), "_", mcols(ArchRProj@peakSet)$idx),
  ]
  rS$start <- start(ArchRProj@peakSet)
  rS$end <- end(ArchRProj@peakSet)
  rS$GC <- ArchRProj@peakSet$GC
  rownames(rS) <- paste0(rS$seqnames, "_", rS$start, "_", rS$end)
  annotationsMatrix <- annotationsMatrix[rownames(rS), ]
  args <- mget(names(formals()), sys.frame(sys.nframe()))
  args$peakAnnotation <- NULL
  rm(peakAnnotation)
  args$annotationsMatrix <- annotationsMatrix
  args$featureDF <- rS
  args$useMatrix <- useMatrix
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$matrixName <- matrixName
  args$X <- seq_along(ArrowFiles)
  args$FUN <- ArchR:::.addDeviationsMatrix
  args$logFile <- logFile
  args$registryDir <- file.path(getOutputDirectory(ArchRProj),
                                paste0(matrixName, "DeviationsRegistry"))
  args$ArchRProj <- NULL
  args$matches <- NULL
  ArchR:::.logThis(args, "Deviations-Args", logFile = logFile)
  outList <- ArchR:::.batchlapply(args)
  ArchR:::.logDiffTime("Completed Computing Deviations!", tstart, addHeader = TRUE,
                       logFile = logFile)
  ArchR:::.endLogging(logFile = logFile)
  gc()
  return(ArchRProj)
}

