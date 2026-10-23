################################################################################
# ARID1A multiome (snRNA-seq) analysis with Seurat
#
# Processes the RNA modality of 10x Genomics Multiome (RNA + ATAC) samples from
# mouse mammary gland: Arid1a WT vs KO, with and without estradiol (E2).
# The ATAC modality is processed separately with ArchR; ArchR is used here only
# to restrict the RNA data to cells that also pass ATAC QC.
#
# Contents
#   1. Configuration (paths, sample labels, palettes, marker gene lists)
#   2. mergeMultiomeSeurat()  - build, QC-filter, normalize and cluster the
#                               merged Seurat object from Cell Ranger output
#   3. main()                 - downstream analysis: UMAPs, markers, DE,
#                               heatmaps, signature scoring, cell classification
#   4. Analysis functions     - ExportToSingleCell(), ScanpyImportReductions(),
#                               refineLuminalBasalClassification(),
#                               SaekiMammaryMarkers(), miloAnalysis(),
#                               hotspot_go(), ...
#   5. Plotting helpers       - UmapPlotEnhance(), DoHeatmapEnhance(),
#                               ViolinPlotEnhance(), MyDotPlot(), ...
#
# main() and mergeMultiomeSeurat() are not meant to be run top to bottom:
# they are exploratory notebooks-in-a-function. Run them section by section
# (each "#### Section ####" header) in an interactive session.
#
# Requires:  R >= 4.2, Seurat v4/v5 (the code uses the @assays slot syntax),
#            see the library() calls below and the "Dependencies" note in the
#            README for the complete package list.
################################################################################

options(bitmapType = 'cairo')

#### Libraries ####
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(data.table)
library(ggrastr)
library(magrittr)

#### Configuration ####
# All locations are configurable. Set the environment variables below (e.g. in
# ~/.Renviron) or edit the defaults. Relative defaults are resolved against the
# directory the script is run from.
#   ARID1A_BASE_DIR   project root; robj/ and plot/ are created beneath it
#   ARID1A_DATA_DIR   input data (marker spreadsheets, scanpy exports, ...)
#   ARID1A_SHEETS_DIR output spreadsheets
#   ARID1A_ARCHR_DIR  ArchR project used to keep only cells that pass ATAC QC
#   ARID1A_RAW_DIR    Cell Ranger Arc output; one sub-directory per sample
#   ARID1A_ANNOT_DIR  annotation files (ribosomal genes, cell-cycle genes, GMTs)
#   ARID1A_LIB_DIR    directory holding venn_helpers.R (default: working directory)
getDir <- function(var, default) {
  path <- Sys.getenv(var, unset = default)
  normalizePath(path, mustWork = FALSE)
}
base.dir         <- getDir('ARID1A_BASE_DIR', '.')
robj.dir         <- file.path(base.dir, 'robj')
plot.dir         <- file.path(base.dir, 'plot')
data.dir         <- getDir('ARID1A_DATA_DIR', file.path(base.dir, 'data'))
spreadsheets.dir <- getDir('ARID1A_SHEETS_DIR', file.path(base.dir, 'spreadsheets'))
archr.dir        <- getDir('ARID1A_ARCHR_DIR', file.path(base.dir, 'ArchR'))
raw.dir          <- getDir('ARID1A_RAW_DIR', file.path(base.dir, 'cellranger'))
annot.dir        <- getDir('ARID1A_ANNOT_DIR', file.path(base.dir, 'annotation'))
lib.dir          <- getDir('ARID1A_LIB_DIR', '.')
for (d in c(robj.dir, plot.dir, spreadsheets.dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Merged, normalized and clustered Seurat object produced by mergeMultiomeSeurat()
arid1a.object <- file.path(robj.dir, 'arid1a_sctransform_multiome_regress.rds')

# gProfiler custom-organism IDs for the Hallmark and GSEA G5 gene sets used in
# hotspot_go(). Create your own with gprofiler2::upload_GMT_file() (see the
# commented calls in hotspot_go()) and export the returned IDs:
#   GPROFILER_HALLMARK_ID, GPROFILER_G5_ID

# all murine specific
bach_marker_list <- list(
  'Bach HS.progenitor'=c('Esr1','Prlr','Pgr','S100a6','Cited1','Aldh1a3','Cd14','Kit'),
  'Bach HS.differentiated'=c('Esr1','Prlr','Pgr','S100a6','Cited1'),
  'Bach Lum.progenitor'=c('Aldh1a3','Cd14','Kit'),
  'Bach Alv.differentiated'=c('Wap','Lalba','Csn2','Glycam1'),
  'Bach Alv.progenitor'=c('Wap','Lalba','Csn2','Glycam1','Aldh1a3','Cd14','Kit'),
  'Bach Basal'=c('Krt4','Krt14','Pdpn','Etv5','Acta2'),
  'Bach Myo'=c('Oxtr','Acta2','Krt4','Krt14'),
  'Bach Procr'=c('Procr','Igfbp4','Gng11','Zeb2')
)
markers_list <- list(
  #https://www.nature.com/articles/s41467-020-15548-7#Sec2 (Supplemental figure 7)
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4331521/
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5696379/ - Visavader
  'Mammary Luminal progenitor'=c('Anpep','Brca1','Cd14','Cd47','Elf5','Egfr','Ltf','Mfge8','Prom1','Tnfrsf11a'),
  'Mammary Luminal mature'=c('Cd24a','Epcam',"Erbb2","Esr1","Foxa1","Gata3",'Krt8','Krt18',"Krt19","Muc1","Pgr","Trp53"),
  'Mammary basal'=c('Anxa8','Itga6','Krt4','Krt5','Krt14','Krt15','Krt17','Mki67','Trp63','Trp73'),
  'mammary epithelium 2-week'=c('Sparc','Tgfbr3','Lgals1','Sorbs1','Timp3','Scn7a','Vim','Cadm1','Prdm16'), #Rxrg
  'mammary epithelium mature luminal'=c('Areg','Cited1','Cldn4','Epcam','Krt8','Ly6d','Prlr','Ptn','Stc2','Krt18'),
  'mammary epithelium progenitor'=c('Aldh1a3','Cd14','Cd81','Csn1s1','Csn3','Elf5','Jund','Lalba','Muc15','Tspan8','Wnt5b'),
  'mammary epithelium basal'=c('Fstl1','Palm','Cd55','Adamts2','Lrp1','Col5a1','Tpm2','Myl9','Myh11','Cnn1','Krt14','Acta2',
                               'Krt5','Cxcl14','Krt17','Snai2','Col18a1','Vcan','Nrg1','Lama1','Wif1','Bmp7','Tgfb1i1',
                               'Cacna2d1','Tenm2','Synpo2'), #'Serpinb11', 'Mrvi1'
  'AR dependent'=c('Ar','Fkbp5','Hpgd','Fgfr2','Lgr4','Egf','Egfr'),
  'Cell Cycle'=c("Pcna","Mki67","Top2a"),
  ### Hallmark gene sets for breast
  'HALLMARK_ESTROGEN_RESPONSE_EARLY'=c('Slc19a2','Adcy1','Cish','Dhcr7','Slc26a2','Abca3','Fkbp5','Flnb','Med13l','Hspb8','B4galt1','Gja1','Gla','Tubb2b','Ar','Jak2','Mapt',
                                       'Myc','Hr','Ovol2','Rara','Rbbp8','Ccnd1','Bcl2','Ret','Scnn1a','Slc1a1','Reep1','Slc22a5','Wfs1','Mlph','Calcr','Papss2','Dlc1','Cant1',
                                       'Abat','Foxc1','Krt13','Msmb','Retreg1','Elovl5','Cxcl12','Tmprss3','Slc2a1','Slc16a1','Sult2b1','Slc9a3r1','Opn3','Gab2','Med24','Il17rb',
                                       'Krt8','Xbp1','Kcnk5','Ptges','Pdzk1','Cldn7','Aqp3','Rasgrp1','Tob1','Mybbp1a','Slc27a2','Adcy9','Add3','Celsr2','Egr3','Elf1','Elf3',
                                       'Fasn','Fdft1','Fhl2','Fkbp4','Fos','Frk','Amfr','Gfra1','Sfn','Hes1','Igf1r','Igfbp4','Il6st','Inhbb','Areg','Krt15','Krt18','Krt19',
                                       'Ablim1','Arl3','Aff1','Muc1','Myb','Mybl1','Npy1r','P2ry2','Pgr','Pmaip1','Podxl','Bag1','Rps6ka2','Siah2','Slc7a2','Sox3',
                                       'Tfap2c','Tff3','Tgm2','Klf10','Tiam1','Tpbg','Ugcg','Zfp185','Sema3b','Calb2','Slc7a5','Akap1','Nrip1','Bhlhe40','Stc2','Cbfa2t3',
                                       'Pex11a','Ccn5','Unc119','Syngr1','Sh3bp5','Scarb1','Cd44','Ncor2','Celsr1','Greb1','Bcl11b','Slc1a4','Inpp5f','Nxt1','Lrig1',
                                       'Slc24a3','Dhrs3','Klf4','Olfm1','Clic3','Pdlim3','Mpped2','Syt12','Cyp26b1','Anxa9','Plaat3','Sec14l2','Abhd2','Rhobtb3',
                                       'Dynlt3','Wwc1','Chpt1','Farp1','Svil','Blvrb','Mindy1','Slc37a1','Asb13','Tbc1d30','Kazn','Myof','Nadsyn1','Tpd52l1','Snx24',
                                       'Deptor','Nav2','Thsd4','Itpk1','Slc39a6','Klk10','Car12','Kdm4b','Rab17','Tjp3','Rhod','Sybu','Isg20l2','Olfml3','Lad1','Tgif2','Tsku',
                                       'Ppif','Fcmr','Esrp2','Aldh3b1','Rrp12','Tmem164','Elovl2','Tiparp','Ttc39a','Dhrs2','Endod1','Kcnk15','Rab31','Prss23','Mast4',
                                       'Fam102a','Rapgefl1','Mreg'),
  'HALLMARK_ESTROGEN_RESPONSE_LATE'=c('Rnaseh2a','Cish','Dhcr7','Slc26a2','Abca3','Etfb','Fgfr3','Fkbp5','Flnb','Gale','Hspb8','Gjb3','Gla','Hprt','Idh2','Jak2','Lamc2','Mapt',
                                      'Ascl1','Mocs2','Ass1','Hr','Ovol2','Rbbp8','Ccnd1','Bcl2','Ret','Scnn1a','Stil','Slc22a5','Wfs1','Xrcc3','Car2','Calcr','Acox2',
                                      'Cav1','Papss2','Large1','Cdh1','Plk4','Emp2','Aldh3a2','Foxc1','Snx10','Igsf1','Krt13','Gal','Dcxr','Elovl5','Cxcl12',
                                      'Tmprss3','Slc16a1','St14','Sult2b1','Tnnc1','Slc9a3r1Homer2','Cdc6','Opn3','Agr2','Slc2a8','Il17rb','Xbp1','Plac1',
                                      'Mettl3','Kcnk5','Ptges','Isg20','Pdzk1','Tob1','Nmu','Btg3','Slc27a2','Klk11','Pkp3','Ckb','Add3','Cpe','Dusp2',
                                      'Celsr2','Egr3','Slc29a1','Fdft1','Fkbp4','Fos','Frk','Amfr','Sfn','Hmgcs2','Id2','Igfbp4','Il6st','Jak1',
                                      'Areg','Krt19','Arl3','Ltf','Mdk','Mest','Aff1','Myb','Nab2','Npy1r','Serpina5Pcp4','Pgr','Prkar2b','Mapk13',
                                      'Prlr','Bag1','Ptger3','Ptpn6','Rps6ka2','S100a9','Sgk1','Siah2','Sord','Sox3','Tfap2c','Tff3','Th','Tiam1',
                                      'Top2a','Tpbg','Tst','Ugdh','Zfp36','Sema3b','Tfpi2','Slc7a5','Nrip1','Ccn5','Cacna2d2','Cd9','Scarb1','Cd44',
                                      'Ncor2','Unc13b','Cxcl14','Slc1a4','Chst8','Scube2','Dnajc12','Lsr','Nxt1','Perp','Trim29','Slc24a3','Impa2',
                                      'Cdc20','Pdcd4','Klf4','Kif20a','Olfm1','Ccna1','Clic3','Dlg5','Pdlim3','Plxnb1','Fabp5','Cyp26b1','Anxa9',
                                      'Plaat3','Abhd2','Dynlt3','Batf','Chpt1','Farp1','Blvrb','Gfus','Dnajc1','Gper1','Myof','Tspan13','Tpd52l1',
                                      'Rabep1','Llgl2','Itpk1','Klk10','Car12','Cox6c','St6galnac2','Tjp3','Ppif','Hspa4l','Aldh3b1','Dhrs2','Gins2',
                                      'Rab31','Prss23','Fam102a','Atp2b4','Rapgefl1','Tpsb2'),
  'HALLMARK_ANDROGEN_RESPONSE' = c('Mertk','Slc26a2','Akt1','Fkbp5','Gpd1l','B4galt1','Gucy1a1','Hpgd','Maf','Ngly1','Ank','Ccnd1','Rps6ka3',
                                   'Bmpr1b','Sms-ps','Actn1','Cdk6','Ndrg1','Aldh1a3','Mak','B2m','Elovl5','Homer2','Sec24d','Map7','Krt8',
                                   'Adamts1','Ube2j1','Ncoa4','Adrm1','Pdlim5','Dbi','Elk4','Ptk2b','Xrcc6','Gnai3','Gsr','H1f0','Hmgcr',
                                   'Hmgcs1','Itgav','Krt19','Lifr','Nkx3-1','Pgm3','Azgp1','Rab4a','Sat1','Selenop','Sgk1','Sord','Srf',
                                   'Tpd52','Ube2i','Xrcc5','Pias1','Tsc22d1','Ccnd3','Vapa','Akap12','Tmprss2','Spdef','Acsl3','Dhcr24',
                                   'Stk39','Dnajb9','Ptpn21','Lman1','Herc3','Camkk2','Inpp4b','Abhd2','Slc38a2','Tmem50a','Insig1','Idi1',
                                   'Abcc4','Tnfaip8','Arxes2','Hsd17b14','Srp19','Uap1','Appbp2','Pmepa1','Pa2g4','Cdc14b','Arid5b','Ell2',
                                   'Fads1','Zbtb10','Cenpn','Zmiz1','Rrp12','Steap4','Iqgap2','Plpp1')
)

markers_list[['Mammary Luminal mature']] <- sort(unique(c(markers_list[['Mammary Luminal mature']],
                                                          markers_list[['mammary epithelium mature luminal']])))
markers_list[['Mammary Luminal progenitor']] <- sort(unique(c(markers_list[['Mammary Luminal progenitor']],
                                                              markers_list[['mammary epithelium progenitor']])))
markers_list[['Mammary basal']] <- sort(unique(c(markers_list[['Mammary basal']])))

# AR targets https://www.ncbi.nlm.nih.gov/pmc/articles/PMC129364/
ar_targets_list = list(
  'Metabolism' = c('Dhcr24','Idi1','Scd2','Scd4','Sms','Phyh','Hpgd','Aldh1a3','Hmgcs1','Sord','Ugdh','B4galt1','Dbi'), #Ddc
  'Proliferation/differentiation/apoptosis' =  c('Ell2','Maf','Ndrg1'), #,'Cdk8','Hes6','Id2','Ctbp1'
  'Transport/trafficking' = c('Slc26a2','Abcc4','Srp19','Ank','Lman1','Vapa','Sec24d','Fkbp5','Snap25','Azgp1','Appbp2'), #'Scnn1a'
  'Transcription regulation' = c('Nkx3-1','Klf4','Arid5b'), #Myc MRF2=Arid5b
  'Protease/protease Inhibitor' = c('Klk4','Cpd','Tmprss2'
                                    #'Klk1b8','Klk1b9',
                                    #"Klk1b11","Klk1b27","Klk1b21","Klk1b24","Klk1b3",
                                    #'Klk1b3','Klk1b16','Klk1b26','Klk1b4','Klk1b5', #Serpin1
  ),
  'Signal transduction' = c('Lifr','Inpp4b','Camkk2','Ptpn21','Tpd52','Gm49337','Akap12','Iqgap2','Sdc4','Mertk'), #'Bche', 'Peg3', 'Pik3r3',
  'Structure/motility/adhesion' = c('Adrm1','Krt19','Actn1','Krt8','Gm9774','Itgav','Dctn3'), #'Mmp16','Fn1'
  'Stress response' = c('Orm1','Dnajb9','Gsr'), #'Orm3','Orm2'
  'Other functions' = c('Rbm10','Uap1','Herc3','B2m') #'St7'
)
nelson_mouse_andr_up <- c("Slc26a2","Orm1","Orm3","Orm2","Gsr","Aldh1a3","Gucy1a1","Hpgd","Actn1","Azgp1","Cpd","B4galt1","Klk1b8","Klk1b1","Klk1b9","Klk1b11","Klk1b26","Klk1b27","Klk1b21","Klk1b22","Klk1b16","Klk1b24","Klk1b3","Klk1b4","Klk1b5","Klk1","Hmgcs1","Id2","Itgav","Krt8","Krt19","Lifr","Sat1","Snap25","Sord","Uap1","Srp19","Ugdh","Vapa","Cdc14b","Plpp1","Inpp4b","B2m","Fkbp5","Klf4","Acsl3","Idi1","Rab4a","Sms","Homer2","Klk4","Scd2","Scd4","Tpd52","Akap12","H1f0","Maf","Selenop","Lman1","Sgk1","Tmprss2","Rbm10","Abcc4","Tsc22d1","Ndrg1","Nkx3-1","Phyh","Unc13b","Appbp2","Pdlim5","Camkk2","Iqgap2","Adamts1","Adrm1","Gm9774","Abhd2","Ptpn21","Pakap","Dctn3","Ell2","Dnajb9","Spdef","Stk39","Tnfaip8","Herc3","Dhcr24","Sec24d","Pgm3","Pias1","Hsd17b14","Hes6","Slc38a2","Pmepa1","Dbi","Elovl5","Zbtb10","Steap4","1700025G04Rik")

#### Color palettes (continuous) ####
rocket_palette <- c("linen", viridis::viridis(direction = -1, option = "F", n=20))

solarExtraWhite_palette <- c("#3361A5","#248AF3","#14B3FF","#88CEEF","white","#EAD397","#FDB31A","#E42A2A","#A31D1D")

color_blind_parula_palette <- colorRampPalette(c("deepskyblue4","deepskyblue2","deepskyblue","grey90","goldenrod1","goldenrod2","goldenrod4"), space="rgb")(120)

blue_gray_purple_diverging_palette <- c(rev(pals::brewer.pubu(10)[-c(1:2)]),'grey99',pals::brewer.purd(10)[-c(1:2)]) # pals::coolwarm(25),

godsnot_102 = c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#6A3A4C","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900",
                "#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75",
                "#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1",
                "#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72")

#### Utilities ####
# Write an R list of gene sets to a .gmt file (function by Levi Waldron, from the
# GSEA helper collection).
writeGMT <- function #Create a gmt (gene matrix transposed) file
### Createss a gmt (gene matrix transposed) file such as those
### provided by mSigDB or geneSigDB, from an R list object.
### Function by Levi Waldron.
(object,
 ### R list object that will be converted to GMT file.  Each element
 ### should contain a vector of gene names, and the names of the
 ### elements will used for the gene set names
 fname
 ### Output file name for .gmt file
){
  if (class(object) != "list") stop("object should be of class 'list'")
  if(file.exists(fname)) unlink(fname)

  #d <- rbindlist(lapply(names(object), function(i){ data.table(i,i, list(object[[i]])) }))
  d <- rbindlist(lapply(names(object), function(i){
    data.table(i,i, c(object[[i]]))
  }))

  fwrite(d, file=fname, sep=" ", sep2 = c("","\t",""),
         quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)
  d
}

INPUT_LABELS = list(
  "MMO1_CRE_Control_multiome" = "Arid1a(-/-)",
  "MMO1_Control_multiome" = "Arid1a(+/+)",
  "MMO1_CRE_E2_Estrogen_multiome" = "Arid1a(-/-) E2",
  "MMO1_E2_Estrogen_multiome" = "Arid1a(+/+) E2"
  # "MMO2_CRE_E2_Estrogen_multiome" = "Arid1a(-/-) E2 mm02"
)

CONDITION_LABELS <- list(
  "MMO1_CRE_Control_multiome" = "KO",
  "MMO1_Control_multiome" = "WT",
  "MMO1_CRE_E2_Estrogen_multiome" = "KO",
  "MMO1_E2_Estrogen_multiome" = "WT"
  # "MMO2_CRE_E2_Estrogen_multiome" = "KO"
)

TREATMENT_LABELS <- list(
  "MMO1_CRE_Control_multiome" = "None",
  "MMO1_Control_multiome" = "None",
  "MMO1_CRE_E2_Estrogen_multiome" = "E2",
  "MMO1_E2_Estrogen_multiome" = "E2"
  #"MMO2_CRE_E2_Estrogen_multiome" = "E2"
)

#sample_palette <- c('#B27BB4', '#00B050', '#34BDEF', '#F3AB19') #,'#E60012') # removed red
sample_palette <- c('#EE1C24','#03A45D','#203F9A','#F9A619')
names(sample_palette) <- names(INPUT_LABELS)
bionames_palette <- sample_palette
names(bionames_palette) <- unlist(INPUT_LABELS[names(bionames_palette)], use.names = F)
bionames_order <- c('Arid1a(+/+)',
                    'Arid1a(+/+) E2',
                    'Arid1a(-/-)',
                    'Arid1a(-/-) E2')

seurat_cluster_palette <- c(as.vector(pals::alphabet(18)))
names(seurat_cluster_palette) <- c(1:18)
condition_palette <- c('#CBA00D','#272E6A','#D51F26')
names(condition_palette) <- c("KO","WT","Mutant")

saeki_mammary_scores <- list("masc_mammary_scores" = "MaSC",
                             "basal_mammary_scores" = "Basal",
                             "la_pro_mammary_scores" = "LA-Pro",
                             "l_alv_mammary_scores" = "L-Alv",
                             "lh_pro_mammary_scores" = "LH-Pro",
                             "l_hor_mammary_scores" = "L-Hor")

mammary_cluster_class_palette <- c("#F0A0FF","#0075DC","#F3C300","#33FFCC","#875692") #, "#CCCCCC")
names(mammary_cluster_class_palette) <- c("Basal/MaSC","L-Alv","L-Hor","LA-Pro","LH-Pro") #, NA)
mammary_cluster_class_palette_list <- as.list(mammary_cluster_class_palette)

mammary_refine_clusters_palette <- c("#F0A0FF","#0075DC","#F3C300","#33FFCC","#875692","#bd0000") #, "#CCCCCC")
names(mammary_refine_clusters_palette) <- c("Basal/MaSC","L-Alv","L-Hor","LA-Pro","LH-Pro","L-Hor-Arid1a-") #, NA)
mammary_refine_clusters_palette_list <- as.list(mammary_refine_clusters_palette)

treatment_palette <- c("#B0A0FC", "#A3C300")
names(treatment_palette) <- c('E2','None')

luminal_basal_palette <- c("#33FFCC","#F3C300","#875692") # is pals::kelly()-1

names(luminal_basal_palette) <- c("Mammary basal","Mammary Luminal mature","Mammary Luminal progenitor")

heatmap_palette <- colorRampPalette(c("royalblue4","royalblue2","deepskyblue",
                                     "grey90","orangered3","red3","darkred"), space="rgb")(1024);
#' Build the merged, QC-filtered and clustered Seurat object
#'
#' Reads the RNA matrix (filtered_feature_bc_matrix.h5) of every multiome sample
#' found under `raw.dir`, merges the samples, annotates them (genotype, treatment,
#' QC percentages), filters low-quality cells, restricts to cells that also pass
#' ATAC QC in the ArchR project (`archr.dir`), then runs SCTransform (with
#' mitochondrial, library size and cell-cycle regression), PCA, UMAP, clustering
#' at multiple resolutions (with silhouette scores) and renumbers/merges clusters.
#' Intermediate QC plots are written to `plot.dir`; the final object is saved to
#' `arid1a.object`.
#'
#' Run section by section (see the "####" headers); the sections depend on one
#' another through objects held in the function environment.
#'
#' @return Nothing; the Seurat object is written to `arid1a.object`.
mergeMultiomeSeurat <- function(){
  min_cells <- 10 # "min_cells" parameter is the minimum number of cells where the gene must be detected in order to be included in the data scaling
  multiome.data <- list.files(path = raw.dir, 'multiome', full.names = T)
  (multiome.sample <- basename(multiome.data))
  (rna.counts <- list.files(path = multiome.data, pattern = '^filtered_feature_bc_matrix.h5$', recursive = T, full.names = T))
  # (atac.frags <- list.files(path = multiome.data, pattern = '^atac_fragments.tsv.gz$', recursive = T, full.names = T))
  sample.list <- lapply(multiome.sample, function(x){
    list('rna'=grep(x, rna.counts, value = T))
    #'atac'=grep(x, atac.frags, value = T))
  }) %>% set_names(multiome.sample)

  print("Found the following samples")
  print(unlist(sample.list))
  print("Excluding sample MM02")
  sample.list <- sample.list[-which(names(sample.list)=='MMO2_CRE_E2_Estrogen_multiome')]
  print("Processing the following samples")
  print(unlist(sample.list))

  #### 1. Create seurat object ####
  for(i in names(sample.list)){
    x <- sample.list[[i]]
    mat <- Read10X_h5(x$rna)
    sample.list[[i]]$seurat_rna <- CreateSeuratObject(counts =  mat$`Gene Expression`, project = i, min.cells = min_cells)
  }
  names(sample.list)

  #### 2. Merge unnormalized seurat objects ####
  # https://satijalab.org/seurat/articles/merge_vignette.html
  arid1a.multiome <- merge(x = sample.list$MMO1_Control_multiome$seurat_rna,
                          y = c(sample.list$MMO1_CRE_Control_multiome$seurat_rna,
                                sample.list$MMO1_CRE_E2_Estrogen_multiome$seurat_rna,
                                sample.list$MMO1_E2_Estrogen_multiome$seurat_rna),
                          add.cell.ids = c("MMO1_Control_multiome",
                                           "MMO1_CRE_Control_multiome",
                                           "MMO1_CRE_E2_Estrogen_multiome",
                                           "MMO1_E2_Estrogen_multiome"),
                          project = 'Arid1a')
  dim(arid1a.multiome@assays$RNA@counts)

  #### 3. Attributes ####
  ## Better for Integration in ArchR
  arid1a.multiome <- RenameCells(object = arid1a.multiome, new.names = gsub("multiome_","multiome#",Cells(arid1a.multiome)) )

  ## Check for cells
  # "MMO1_CRE_E2_Estrogen_multiome#GCAGGACCAAGCTTAT-1" "MMO1_CRE_E2_Estrogen_multiome#GCCTCGACACCCTCAC-1" "MMO1_CRE_E2_Estrogen_multiome#CTTCTAAGTTAACGAT-1"
  # "MMO1_CRE_E2_Estrogen_multiome#GTCTAATCACTCAACA-1" "MMO1_CRE_E2_Estrogen_multiome#ACCAGGGAGTGTTGTA-1" "MMO1_CRE_E2_Estrogen_multiome#ACACTTAAGTGCACGC-1"
  # "MMO1_CRE_E2_Estrogen_multiome#CGTACGGGTATTGAGT-1" "MMO1_CRE_E2_Estrogen_multiome#TGTGATCAGCATTAGC-1" "MMO1_CRE_E2_Estrogen_multiome#CGGTTCCGTCTAACCT-1"
  # "MMO1_CRE_E2_Estrogen_multiome#ACTCGCGCACAAGCCT-1" "MMO1_CRE_E2_Estrogen_multiome#CTCCTGAGTTCATTTG-1" "MMO1_CRE_E2_Estrogen_multiome#GTTAATGTCAGCAAGA-1"
  # "MMO1_CRE_E2_Estrogen_multiome#GAGTGATCATAGACCC-1" "MMO1_CRE_E2_Estrogen_multiome#GTGTCCAAGTTGCGTG-1" "MMO1_CRE_E2_Estrogen_multiome#ATTTAGGTCGTTTCGC-1"
  # arid1a.multiome@meta.data['MMO1_CRE_E2_Estrogen_multiome#GTTAATGTCAGCAAGA-1',]

  ## Bionames
  arid1a.multiome@meta.data$bioNames <- sapply(arid1a.multiome@meta.data$orig.ident, function(i) INPUT_LABELS[[i]])
  table(arid1a.multiome@meta.data$bioNames)
  ## Condition
  arid1a.multiome@meta.data$condition <- sapply(arid1a.multiome@meta.data$orig.ident, function(i) CONDITION_LABELS[[i]])
  table(arid1a.multiome@meta.data$condition)
  ## Treatment
  arid1a.multiome@meta.data$treatment <- sapply(arid1a.multiome@meta.data$orig.ident, function(i) TREATMENT_LABELS[[i]])
  table(arid1a.multiome@meta.data$treatment)
  ## Group
  arid1a.multiome@meta.data$group <- paste(arid1a.multiome@meta.data$condition, arid1a.multiome@meta.data$treatment,sep='_')
  table(arid1a.multiome@meta.data$group)

  identical(arid1a.multiome@assays$RNA@counts, arid1a.multiome@assays$RNA@data)

  head(colnames(arid1a.multiome))
  tail(colnames(arid1a.multiome))
  table(gsub('_multiome.*',"",colnames(arid1a.multiome)))
  table(arid1a.multiome$orig.ident)

  #### 3. Percentages ####
  # https://satijalab.org/seurat/articles/sctransform_vignette.html
  # store mitochondrial percentage in object meta data
  ribosomal_genes = fread(file.path(annot.dir, 'mouse_ribosomal_genes.txt'), header = F)
  ribosomal_genes <- paste(paste0("^",unname(unlist(ribosomal_genes))), collapse = "|")
  grep("^mt-",rownames(arid1a.multiome), ignore.case = T, value = T)
  grep(ribosomal_genes,rownames(arid1a.multiome), ignore.case = T, value = T)
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^mt-", col.name = "percent_mito")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = ribosomal_genes, col.name = "percent_ribo")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Malat1", col.name = "percent_malat1")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Arid1a", col.name = "percent_arid1a")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Foxa1", col.name = "percent_foxa1")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Krt4$", col.name = "percent_krt4")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Krt5$", col.name = "percent_krt5")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Krt8$", col.name = "percent_krt8")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Krt14$", col.name = "percent_krt14")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Krt18$", col.name = "percent_krt18")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Krt19$", col.name = "percent_krt19")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Trp63$", col.name = "percent_trp63")
  arid1a.multiome <- PercentageFeatureSet(arid1a.multiome, pattern = "^Epcam", col.name = "percent_epcam")

  #### 4. Basic QC filtering ####
  # Visualize QC metrics as a violin plot
  MAX.MITO.PERC <- 25
  MIN.FEATURE <- 500 #2000
  MIN.RNA.COUNT <- 300 #2000
  MAX.RNA.COUNT <- 100000
  dt_cutoff <- data.table(variable=c('percent_mito',
                                     'nFeature_RNA',
                                     'nCount_RNA_MIN',
                                     'nCount_RNA'), y=c(MAX.MITO.PERC,
                                                        MIN.FEATURE,
                                                        MIN.RNA.COUNT,
                                                        MAX.RNA.COUNT))

  dtplot <- data.table(arid1a.multiome@meta.data)[,list(bioNames, nFeature_RNA, nCount_RNA, percent_malat1,
                                                       percent_mito, percent_ribo, percent_foxa1, percent_arid1a,
                                                       percent_krt8, percent_krt5, percent_trp63, percent_krt14)]
  dtplot$nFeature_RNA <- as.numeric(dtplot$nFeature_RNA)
  dtplot <- melt(dtplot, id.vars = "bioNames")
  png(filename = file.path(plot.dir, 'pre_processing_violin.png'), width = 15, height = 7, bg = 'transparent', res = 150, units = 'in')
  ggplot(dtplot, aes(x=bioNames, y=value, fill=bioNames, label=variable)) +
    geom_violin(trim = F, scale = 'width') +
    scale_fill_manual("", values = unname(sample_palette)) +
    xlab("") + ylab("") +
    geom_boxplot(width=.4, fill='white', alpha=0.6, outlier.size = 0.0) +
    geom_hline(data = dt_cutoff, aes(yintercept=y), color='red', linetype='dashed') +
    facet_wrap(~ variable, nrow = 2, scales = 'free_y') + theme_bw(14) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          panel.grid = element_blank(), panel.background = element_blank(), plot.background = element_blank(), legend.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = 'bold'), )
  dev.off()

  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                          group.by = 'bioNames', pt.size = 0.5, cols = bionames_palette) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='nFeature_RNA',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())

  plot2 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "percent_arid1a",
                          group.by = 'bioNames', cols = bionames_palette, pt.size = 0.5) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='percent_foxa1',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())

  plot3 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "percent_malat1",
                          group.by = 'bioNames', cols = bionames_palette, pt.size = 0.5) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='percent_malat1',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())

  plot4 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "percent_mito",
                          group.by = 'bioNames', cols = bionames_palette, pt.size = 0.5) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='percent_mito',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())

  png(filename = file.path(plot.dir, 'pre_processing_scatter.png'), width = 18, height = 8, bg = 'transparent', res = 150, units = 'in')
  grid.arrange(plot1, plot2, plot3, plot4, nrow=2)
  dev.off()

  plot5 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                          group.by = 'bioNames', pt.size = 0.5, cols = bionames_palette) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='nFeature_RNA',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank()) + facet_wrap( ~ colors)

  png(filename = file.path(plot.dir, 'pre_processing_scatter2.png'), width = 9, height = 6, bg = 'transparent', res = 150, units = 'in')
  plot5
  dev.off()

  g1 <- ggplot(dtplot[variable %in% 'nCount_RNA' & between(value, MIN.RNA.COUNT, MAX.RNA.COUNT)],
         aes(x=value, color=bioNames, label=variable)) +
    geom_density() +
    scale_color_manual("", values = unname(sample_palette)) +
    xlab("") + ylab("") +
    facet_wrap(~ variable, nrow = 2, scales = 'free') + theme_bw(14) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          strip.background = element_blank(), strip.text = element_text(face = 'bold'), )

  g2 <- ggplot(dtplot[variable=='nFeature_RNA'], aes(x=value, color=bioNames)) +
    geom_density() + xlab("nFeature_RNA") +
    geom_vline(xintercept = MIN.FEATURE) +
    scale_color_manual(values = unname(sample_palette)) + theme_bw()

  png(filename = file.path(plot.dir, 'pre_processing_scatter3.png'), width = 9, height = 9, bg = 'transparent', res = 150, units = 'in')
  grid.arrange(g1, g2, nrow=2)
  dev.off()

  #### Save object ####
  arid1a.object
  saveRDS(arid1a.multiome, file = arid1a.object)

  #### filter counts ####
  lowRNA <- WhichCells(arid1a.multiome, expression = nFeature_RNA < dt_cutoff[variable=='nFeature_RNA',y])
  lowCounts <- WhichCells(arid1a.multiome, expression = nCount_RNA < dt_cutoff[variable=='nCount_RNA_MIN',y])
  exaggeratedCounts <- WhichCells(arid1a.multiome, expression = nCount_RNA > dt_cutoff[variable=='nCount_RNA',y])
  exaggeratedMito <- WhichCells(arid1a.multiome, expression = percent_mito > dt_cutoff[variable=='percent_mito',y])
  # remove cells
  message("Filtering cells with low RNA transcripts and low counts and exaggerated counts")
  arid1a.multiome <- subset(arid1a.multiome, cells=setdiff(Cells(arid1a.multiome), lowRNA))
  arid1a.multiome <- subset(arid1a.multiome, cells=setdiff(Cells(arid1a.multiome), lowCounts))
  arid1a.multiome <- subset(arid1a.multiome, cells=setdiff(Cells(arid1a.multiome), exaggeratedCounts))
  arid1a.multiome <- subset(arid1a.multiome, cells=setdiff(Cells(arid1a.multiome), exaggeratedMito))

  identical(arid1a.multiome@assays$RNA@counts, arid1a.multiome@assays$RNA@data)

  #### Remove cells not in ATAC ####
  ## Bring in ArchR cells
  proj <- loadArchRProject(path = archr.dir)
  archRcells <- rownames(proj@cellColData)

  table(gsub("#.*","",proj$Sample))
  table(gsub("#.*","",Cells(arid1a.multiome)))
  message("Filtering cells filtered by ATAC")
  arid1a.multiome <- subset(arid1a.multiome, cells=intersect(Cells(arid1a.multiome), archRcells))

  table(gsub("#.*","",Cells(arid1a.multiome)))
  table(gsub("#.*","",proj$Sample))
  ## done filtering

  dtplot <- data.table(arid1a.multiome@meta.data)[,list(bioNames, nFeature_RNA, nCount_RNA, percent_malat1,
                                                        percent_mito, percent_ribo, percent_foxa1, percent_arid1a,
                                                        percent_krt8, percent_krt5, percent_trp63, percent_krt14)]
  dtplot$nFeature_RNA <- as.numeric(dtplot$nFeature_RNA)
  dtplot <- melt(dtplot, id.vars = "bioNames")
  png(filename = file.path(plot.dir, 'post_processing_violin.png'), width = 15, height = 7, bg = 'transparent', res = 150, units = 'in')
  ggplot(dtplot, aes(x=bioNames, y=value, fill=bioNames, label=variable)) +
    geom_violin(trim = F, scale = 'width') +
    scale_fill_manual("", values = unname(sample_palette)) +
    xlab("") + ylab("") +
    geom_boxplot(width=.4, fill='white', alpha=0.6, outlier.size = 0.0) +
    geom_hline(data = dt_cutoff, aes(yintercept=y), color='red', linetype='dashed') +
    facet_wrap(~ variable, nrow = 2, scales = 'free_y') + theme_bw(14) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          panel.grid = element_blank(), panel.background = element_blank(), plot.background = element_blank(), legend.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = 'bold'), )
  dev.off()

  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                          group.by = 'bioNames', pt.size = 0.5, cols = bionames_palette) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='nFeature_RNA',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())

  plot2 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "percent_arid1a",
                          group.by = 'bioNames', cols = bionames_palette, pt.size = 0.5) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='percent_foxa1',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())

  plot3 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "percent_malat1",
                          group.by = 'bioNames', cols = bionames_palette, pt.size = 0.5) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='percent_malat1',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())

  plot4 <- FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "percent_mito",
                          group.by = 'bioNames', cols = bionames_palette, pt.size = 0.5) +
    scale_x_log10() +
    geom_hline(yintercept = dt_cutoff[variable=='percent_mito',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA_MIN',y], color='red', linetype='dashed') +
    geom_vline(xintercept = dt_cutoff[variable=='nCount_RNA',y], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())

  png(filename = file.path(plot.dir, 'post_processing_scatter.png'), width = 18, height = 8, bg = 'transparent', res = 150, units = 'in')
  grid.arrange(plot1, plot2, plot3, plot4, nrow=2)
  dev.off()

  ## Cell cycle
  cell_cycle_markers_df = fread(file.path(annot.dir, 'Human_Cell_Cycle_Macosko_et_al.txt'), sep='\t')
  s.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$S[cell_cycle_markers_df$S!=""]) )
  g2m.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$G2.M[cell_cycle_markers_df$G2.M!=""]) )
  m.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$M[cell_cycle_markers_df$M!=""]) )

  # run sctransform, default is 3000 variable features
  arid1a.multiome
  arid1a.multiome <- SCTransform(arid1a.multiome,
                                 vst.flavor="v2",
                                 assay = 'RNA',
                                 new.assay.name = 'SCT' ,
                                 method = "glmGamPoi",
                                 vars.to.regress = c("percent_mito","nFeature_RNA","nCount_RNA"),
                                 return.only.var.genes = FALSE,
                                 verbose = T)

  ## Regress out signal from SCT assay
  arid1a.multiome <- CellCycleScoring(object = arid1a.multiome, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT')
  ## now that we have SCT normalized cell-cycle scores we run normalization again and regress these out
  arid1a.multiome <- SCTransform(arid1a.multiome,
                                 vst.flavor="v2",
                                 assay = 'RNA',
                                 new.assay.name = 'SCT' ,
                                 method = "glmGamPoi",
                                 return.only.var.genes = FALSE,
                                 vars.to.regress = c('percent_mito', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

  # inspect model
  library(sctransform)
  #### 5. RNA assay ####
  ## Normalize RNA assay
  arid1a.multiome <- NormalizeData(object = arid1a.multiome, assay='RNA', normalization.method = "LogNormalize")
  arid1a.multiome <- CellCycleScoring(object = arid1a.multiome,
                                     s.features = s.genes,
                                     g2m.features = g2m.genes,
                                     assay = 'RNA')
  arid1a.multiome <- FindVariableFeatures(arid1a.multiome, selection.method = "vst", assay = 'RNA')
  ## Regress out signal from RNA assay
  arid1a.multiome <- ScaleData(arid1a.multiome,
                              assay = 'RNA',
                              vars.to.regress = c('percent_mito','nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

  #### Save object ####
  saveRDS(arid1a.multiome, file = arid1a.object)

  # Ensure dimension are all same
  # Must match for import to scanpy
  dim(arid1a.multiome@assays$RNA@counts)
  dim(arid1a.multiome@assays$SCT@counts)

  dim(arid1a.multiome@assays$SCT@data)
  dim(arid1a.multiome@assays$RNA@data)

  dim(arid1a.multiome@assays$SCT@scale.data)
  # okay this isn't same as we don't import it to python
  dim(arid1a.multiome@assays$RNA@scale.data)

  #### 4. Attributes ####
  #### START HERE ####
  arid1a.multiome <- readRDS(arid1a.object)

  ## Better for Integration in ArchR
  Cells(arid1a.multiome)
  ## Bionames
  table(arid1a.multiome@meta.data$bioNames)
  ## Condition
  table(arid1a.multiome@meta.data$condition)

  #### 5. PCA and UMAP ####
  # arid1a.multiome[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA
  # arid1a.multiome[["SCT"]]@counts are corrected UMI counts
  # arid1a.multiome[["SCT"]]@data are log-normalized corrected UMI counts

  # library(future)
  # # max ram allowed here is 1000 * 1024^2 = 1GB
  # options('future.globals.maxSize' = 1000 * 1024^2)
  # # check the current active plan
  # plan()
  # # change the current plan to access parallelization
  # #plan("multiprocess", workers = 2)
  # plan("sequential")
  # plan()
  # ### RNA
  # the final UMAP 'looks' better when using assay = 'RNA' for PCA (only), however it is recommended to use SCT:
  arid1a.multiome.rna <- RunPCA(arid1a.multiome, verbose = FALSE, assay = 'RNA')
  # determine dimensionality
  pdf(file = file.path(plot.dir, 'pca_dimension_elbow_regress_RNA.pdf'), width = 6, height = 6)
  print(ElbowPlot(arid1a.multiome.rna, ndims = 50))
  dev.off()
  png(filename = file.path(plot.dir, 'pca_dimension_heatmap_regress_RNA.png'), width = 8, height = 16, units = 'in', bg = 'transparent', res = 150)
  DimHeatmap(arid1a.multiome.rna, dims=1:15, cells = 2000, assays = 'RNA', balanced = T)
  dev.off()

  # evaluate umap
  DefaultAssay(arid1a.multiome.rna) <- 'RNA'
  arid1a.multiome.rna <- RunUMAP(arid1a.multiome.rna, dims= 1:30, n.neighbors = 30L,
                                 min.dist = 0.3, spread = 1, verbose = FALSE,
                                 assay = 'RNA', metric = 'cosine')

  g1 <- DimPlot(arid1a.multiome.rna, label = F, group.by = 'nCount_RNA') + NoLegend()
  g1$data$nCount_RNA <- as.numeric(g1$data$nCount_RNA)
  UmapPlotEnhance(obj = arid1a.multiome.rna, p1 = g1, addLabels = F, label.size = 10, break.unit = 2000,
                  axis.label = NA, legend.title = 'RNA count', scale.color = F, pal=rocket_palette,
                  filename.prefix = 'UMAP_count_RNA')

  g1 <- DimPlot(arid1a.multiome.rna, label = F, group.by = 'bioNames') + NoLegend()
  UmapPlotEnhance(obj = arid1a.multiome.rna, p1 = g1, addLabels = F, label.size = 10,
                  axis.label = NA, legend.title = '', pt.size = 0.5, pal=bionames_palette,
                  filename.prefix = 'UMAP_bionames_RNA')

  ### SCT
  arid1a.multiome.sct <- RunPCA(arid1a.multiome, verbose = FALSE, assay = 'SCT')
  # determine dimensionality
  pdf(file = file.path(plot.dir, 'pca_dimension_elbow_regress_SCT.pdf'), width = 6, height = 6)
  print(ElbowPlot(arid1a.multiome.sct, ndims = 50))
  dev.off()
  png(filename = file.path(plot.dir, 'pca_dimension_heatmap_regress_SCT.png'), width = 8, height = 16, units = 'in',
      bg = 'transparent', res = 150)
  DimHeatmap(arid1a.multiome.sct, dims=1:15, cells = 2000, balanced = T)
  dev.off()

  # evaluate umap
  DefaultAssay(arid1a.multiome.sct)
  arid1a.multiome.sct <- RunUMAP(arid1a.multiome.sct, dims= 1:30, n.neighbors = 30L, seed.use = 21,
                                 min.dist = 0.3, spread = 1, verbose = FALSE, assay = 'SCT', metric = 'cosine')

  g1 <- DimPlot(arid1a.multiome.sct, label = F, group.by = 'bioNames') + NoLegend()
  UmapPlotEnhance(obj = arid1a.multiome.sct, p1 = g1, addLabels = F, label.size = 10, axis.label = NA,
                  legend.title = '',  pt.size = 0.5,
                  pal=bionames_palette, filename.prefix = 'UMAP_bionames_SCT', as.pdf = F)

  g1 <- DimPlot(arid1a.multiome.sct, label = F, group.by = 'nCount_RNA') + NoLegend()
  g1$data$nCount_RNA <- as.numeric(g1$data$nCount_RNA)
  UmapPlotEnhance(obj = arid1a.multiome.sct, p1 = g1, addLabels = F, label.size = 10, axis.label = NA,
                  legend.title = 'RNA count', scale.color = F, filename.prefix = 'UMAP_count_SCT')

  g1 <- DimPlot(arid1a.multiome.sct, label = F, group.by = 'condition') + NoLegend()
  UmapPlotEnhance(obj = arid1a.multiome.sct, p1 = g1, addLabels = F, label.size = 10, legend.title = '', axis.label = NA,
                  pal=condition_palette, filename.prefix = 'UMAP_day_condition', as.pdf = F)

  g1 <- FeaturePlot(object = arid1a.multiome.sct, features = c('Arid1a'), slot = 'data')
  UmapPlotEnhance(obj = arid1a.multiome.sct, p1 = g1, addLabels = F, label.size = 10, axis.label = '',
                  legend.title = '', scale.color = F,
                  pal=rocket_palette, filename.prefix = 'UMAP_Arid1a', as.pdf = F)

  g1 <- FeaturePlot(object = arid1a.multiome.sct, features = c('Krt8'), slot = 'data')
  UmapPlotEnhance(obj = arid1a.multiome.sct, p1 = g1, addLabels = F, label.size = 10, axis.label = '',
                  legend.title = '', scale.color = F,
                  pal=rocket_palette, filename.prefix = 'UMAP_Krt8', as.pdf = F)

  g1 <- FeaturePlot(object = arid1a.multiome.sct, features = 'Foxa1', slot = 'data')
  UmapPlotEnhance(obj = arid1a.multiome.sct, p1 = g1, addLabels = F, label.size = 10, axis.label = '',
                  legend.title = '', scale.color = F,
                  pal=rocket_palette, filename.prefix = 'UMAP_Foxa1', as.pdf = F)

  g1 <- FeaturePlot(object = arid1a.multiome.sct, features = c('Krt5'), slot = 'data')
  UmapPlotEnhance(obj = arid1a.multiome.sct, p1 = g1, addLabels = F, label.size = 10, axis.label = '',
                  legend.title = '', scale.color = F,
                  pal=rocket_palette, filename.prefix = 'UMAP_Krt5', as.pdf = F)

  #### Choose Assay ####
  arid1a.multiome <- arid1a.multiome.sct

  #### Save object ####
  saveRDS(arid1a.multiome, file = arid1a.object)

  #
  arid1a.multiome <- FindNeighbors(arid1a.multiome,
                                  dims = 1:30,
                                  annoy.metric = 'cosine',
                                  assay = 'SCT',
                                  verbose = FALSE)

  ## Evaluate clustering
  library(future)
  plan()
  res <- seq(from = 0.1, to = 1.5, by = 0.1)
  arid1a.multiome <- FindClusters(arid1a.multiome, verbose = FALSE, algorithm = 4, resolution = res)

  res <- grep("SCT_snn_res", colnames(arid1a.multiome@meta.data), value = T)
  (res <- res[order(gsub("SCT_snn_res.","",res))])

  sils <- lapply(res, function(r){
    sil <- silhouetteScore(obj = arid1a.multiome, reduction = 'umap', clusters =  r, dims = 1:2)
    #summary(sil)$avg.width
  }) %>% set_names(res)

  um <- lapply(res, function(r){
    score <- summary(sils[[r]])$avg.width
    g1 <- DimPlot(arid1a.multiome, label = TRUE, group.by = r) + NoLegend()
    g1 <- UmapPlotEnhance(obj = arid1a.multiome,
                          p1 = g1,
                          label.size = 10,
                          title.size = 12,
                          pal=seurat_cluster_palette,
                          legend.position = "right",
                          title = paste("Cluster res.:", gsub("SCT_snn_res.","",r),
                                        "\nmean silhouette ", formatC(score, digits = 3)))
  }) %>% set_names(res)

  png(filename = file.path(plot.dir, 'silhouette_scores_regress.png'), width = 25, height = 16, units = 'in', res=150, bg = 'transparent')
  do.call(grid.arrange, c(um, nrow=4, top='resolutions'))
  dev.off()

  pdf(file = file.path(plot.dir, 'silhouette_scores_regress.pdf'), width = 10, height = 8, onefile=T)
  for(r in res){
    g <- um[[r]]
    p1 <- ~plot(sils[[r]], main=paste("resolution:",r), cex.names = par("cex.axis")) # RStudio sometimes does not display silhouette plots correctly
    grid.arrange(cowplot::as_grob(p1), g, nrow=1, widths=c(.35,.65))
  }
  dev.off()

  #### Pick clustering resolution ####
  cresolution <- 0.6 # 1.0
  ##arid1a.multiome <- FindClusters(arid1a.multiome, verbose = FALSE, algorithm = 4, resolution = cresolution, method = 'igraph')

  cresolution <- paste0('SCT_snn_res.', cresolution)

  library(clustree)
  # https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
  clustree(arid1a.multiome@meta.data, prefix = "SCT_snn_res.")

  arid1a.multiome@meta.data$seurat_clusters <- arid1a.multiome@meta.data[,cresolution]
  identical(arid1a.multiome@meta.data$seurat_clusters, arid1a.multiome@meta.data[[cresolution]])

  Idents(object = arid1a.multiome) <- cresolution

  # Plot UMAPS
  g1 <- DimPlot(arid1a.multiome, label = TRUE) + NoLegend()
  g1 <- UmapPlotEnhance(obj = arid1a.multiome, p1 = g1, addLabels = T,
                        label.size = 4, axis.label = NA,
                        pal=seurat_cluster_palette, title = "clusters", legend.position = "right",
                        filename.prefix = 'umap_by_cluster_sctassay_final', legend.title = '')

  #### 8. Merge Clusters ####
  ## Merge clusters 3, 6 and merge 7,8 and merge 4,3,9 and merge 10,11
  ## Then renumber clusters based on size
  table(arid1a.multiome@meta.data$seurat_clusters)
  arid1a.multiome@meta.data$seurat_clusters[arid1a.multiome@meta.data$seurat_clusters %in% c(3,6)] <- 3
  arid1a.multiome@meta.data$seurat_clusters[arid1a.multiome@meta.data$seurat_clusters %in% c(7,8)] <- 7
  is.factor(arid1a.multiome@meta.data$seurat_clusters)
  arid1a.multiome@meta.data$seurat_clusters <- as.integer(arid1a.multiome@meta.data$seurat_clusters)
  is.factor(arid1a.multiome@meta.data$seurat_clusters)
  table(arid1a.multiome@meta.data$seurat_clusters)
  cl.order <- sort(table(arid1a.multiome@meta.data$seurat_clusters), decreasing = T)
  # Re number clusters
  for(i in seq_len(length(cl.order))){
    cur.clust <- names(cl.order[i])
    arid1a.multiome@meta.data$seurat_clusters[arid1a.multiome@meta.data$seurat_clusters==cur.clust] <- paste0("C",i)
  }
  table(arid1a.multiome@meta.data$seurat_clusters)
  arid1a.multiome@meta.data$seurat_clusters <- as.integer(gsub("C","",arid1a.multiome@meta.data$seurat_clusters))
  table(arid1a.multiome@meta.data$seurat_clusters)
  arid1a.multiome@meta.data$seurat_clusters <- factor(x = arid1a.multiome@meta.data$seurat_clusters,
                                                      levels = sort(unique(arid1a.multiome@meta.data$seurat_clusters), decreasing=F),
                                                      ordered = T)
  table(arid1a.multiome@meta.data$seurat_clusters)
  unique(arid1a.multiome@meta.data$seurat_clusters)
  unique(arid1a.multiome@meta.data$orig.ident)
  is.factor(arid1a.multiome@meta.data$seurat_clusters)

  Idents(object = arid1a.multiome) <- arid1a.multiome$seurat_clusters

  ## After the fact
  ## Merge 3,4,8
  table(arid1a.multiome@meta.data$seurat_clusters)
  arid1a.multiome@meta.data$seurat_clusters[arid1a.multiome@meta.data$seurat_clusters %in% c(3,4,8)] <- 3
  is.factor(arid1a.multiome@meta.data$seurat_clusters)
  arid1a.multiome@meta.data$seurat_clusters <- as.integer(arid1a.multiome@meta.data$seurat_clusters)
  (cl.order <- sort(table(arid1a.multiome@meta.data$seurat_clusters), decreasing = T))
  # Re number clusters
  for(i in seq_len(length(cl.order))){
    cur.clust <- names(cl.order[i])
    arid1a.multiome@meta.data$seurat_clusters[arid1a.multiome@meta.data$seurat_clusters==cur.clust] <- paste0("C",i)
  }
  table(arid1a.multiome@meta.data$seurat_clusters)
  arid1a.multiome@meta.data$seurat_clusters <- as.integer(gsub("C","",arid1a.multiome@meta.data$seurat_clusters))
  table(arid1a.multiome@meta.data$seurat_clusters)
  arid1a.multiome@meta.data$seurat_clusters <- factor(x = arid1a.multiome@meta.data$seurat_clusters,
                                                      levels = sort(unique(arid1a.multiome@meta.data$seurat_clusters), decreasing=F),
                                                      ordered = T)
  table(arid1a.multiome@meta.data$seurat_clusters)
  unique(arid1a.multiome@meta.data$seurat_clusters)
  unique(arid1a.multiome@meta.data$orig.ident)
  is.factor(arid1a.multiome@meta.data$seurat_clusters)

  Idents(object = arid1a.multiome) <- arid1a.multiome$seurat_clusters

  #### 9. save object ####
  message("Saving seurat object to: ",arid1a.object)
  saveRDS(arid1a.multiome, arid1a.object)

}

# Description of Assay and Slots
# https://github.com/satijalab/seurat/wiki/Assay
# RNA: for gene level : D.E., violin plots
# counts = raw counts
# data = log1p(counts) after normalization divide counts by 10K (for instance)
# scale.data = normalized and scaled for dim reduction
#
# SCT: cell level, UMAP clustering
# counts = (corrected UMI counts if all cells were sequenced with same depth) counts
# data = log1p(counts)
# scale.data = pearson residuals

## Which assay for what?
# https://github.com/satijalab/seurat/issues/4082
# RNA@data : D.E., Feature/Violin plots
# SCT@scale.data : cell level, UMAP clustering
#
# https://github.com/satijalab/seurat/issues/2899
# You can use the corrected Pearson residuals (stored in the scale.data slot of the integrated assay) for joint analysis, like PCA and clustering
# You should never use integrated/corrected values when performing DE.
# We typically use the @data slot (log-normalized corrected UMI counts) instead of scale.data (pearson residuals).
# Not suggested but if you do run DE on the scale.data slot (for any assay), the logFC values are not meaningful.
# Typically use:
# scale.data slot for heat maps
# Data slot for Feature/Violin plots.

#' Downstream analysis of the clustered multiome RNA object
#'
#' Exploratory analysis notebook, organised in "####" sections: UMAP figures,
#' marker expression, differential expression (cluster markers, KO vs WT,
#' E2 vs no E2), heatmaps and dot plots, volcano plots, comparison with sorted
#' bulk populations, AR / ER target scoring, luminal/basal classification and
#' export of the object. Requires the object produced by mergeMultiomeSeurat().
#' Figures are written to `plot.dir`, tables to `spreadsheets.dir`.
#'
#' Run section by section in an interactive session.
main <- function(){
  library(future)
  # check the current active plan
  plan()
  # max ram allowed here is 1000 * 1024^2 = 1GB
  options('future.globals.maxSize' = 1000 * 2048^2)
  # change the current plan to access parallelization
  plan("multisession", workers = 4)
  #plan(strategy = 'sequential')

  #### Load object ####
  arid1a.multiome <- readRDS(file=arid1a.object)

  #### Export to Scanpy ####
  ExportToSingleCell(obj = arid1a.multiome)

  FeatureScatter(arid1a.multiome, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

  g <- ggplot(arid1a.multiome@meta.data) +
    geom_bar(aes(x=bioNames, y= stat(count), fill=bioNames)) +
    scale_fill_manual("", values=bionames_palette) +
    theme_classic(14) + xlab("") + ylab("number of cells") +
    theme(axis.text.x = element_text(angle=90), legend.position = 'none')
  gp <- ViolinPlotQC(arid1a.multiome)
  cairo_pdf(filename = file.path(plot.dir, 'Arid1a_multiome_qc_emptydrop_compare.pdf'), width = 14, height = 4, bg = 'transparent')
  do.call(grid.arrange, c( c(list(g), gp), nrow=1, top='Foxa1 multiome'))
  dev.off()

  arid1a.multiome$atac_clusters <- proj@cellColData[colnames(arid1a.multiome),"Clusters"]

  #### Basic UMAP Plot ####
  glabel <- DimPlot(arid1a.multiome, label = TRUE, group.by = 'seurat_clusters')
  glabel <- UmapPlotEnhance(obj = arid1a.multiome, p1 = glabel, addLabels = T,
                        label.size = 4, axis.label = NA,
                        pal=NA, title = "clusters", legend.title = '')

  g1 <- DimPlot(arid1a.multiome, label = F, group.by = 'bioNames') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, pt.size = 0.9,
                    axis.label = NA, legend.title = '', title='snRNA samples',
                    pal=bionames_palette, panel.border.size = 0)
  pdf(file = file.path(plot.dir,'bionames_umap_basic.pdf'), width = 5, height = 2.6, bg = 'transparent')
  g1 + glabel$layers
  dev.off()

  g1 <- DimPlot(arid1a.multiome, label = F, group.by = 'bioNames') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, panel.border.size = 0,
                    pt.size = 0.9, axis.label = NA, legend.title = '',
                    pal=bionames_palette, show.legend = T, title = "snRNA samples", as.pdf = T,
                    filename.prefix = "umap_bionames")

  g2 <- DimPlot(arid1a.multiome, label = F, group.by = 'seurat_clusters') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., addLabels = T, panel.border.size = 0,
                    label.size = 5, pt.size = 0.9, axis.label = NA, legend.title = '',
                    scale.color = F, pal=seurat_cluster_palette, show.legend = F,
                    title = "snRNA clusters", filename.prefix = "umap_clusters")

  g3 <- DimPlot(arid1a.multiome, label = F, group.by = 'condition') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10,
                    panel.border.size = 0, order.random=F, order.reverse = T,
                    pt.size = 0.9, axis.label = NA, legend.title = '',
                    as.pdf = T,
                    pal=condition_palette, title = "Genotype",
                    filename.prefix = "umap_condition")

  g4 <- DimPlot(arid1a.multiome, label = F, group.by = 'condition') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    label.size = 10,
                    panel.border.size = 0,
                    order.random = T,
                    order.reverse = T,
                    pt.size = 0.9,
                    axis.label = "condition",
                    legend.title = '',
                    as.pdf = T,
                    nrow=1,
                    pal=condition_palette,
                    title = "Genotype",
                    filename.prefix = "umap_condition_facet")

  cluster_palette <- c(ArchRPalettes$stallion, ArchRPalettes$stallion2)
  names(cluster_palette) <- paste0("C",seq_along(cluster_palette))
  levels(cluster_palette) <- names(cluster_palette)
  cluster_palette_list <- as.list(cluster_palette)

  g4 <- DimPlot(arid1a.multiome, label = F, group.by = 'atac_clusters')
  g4$data <- g4$data[g4$data$atac_clusters!="C5",]
  UmapPlotEnhance(obj = arid1a.multiome, p1 = g4, addLabels = F, panel.border.size = 0,
                  label.size = 5, pt.size = 0.9, axis.label = NA, legend.title = '',
                  scale.color = F, pal=cluster_palette_list, show.legend = F, title = "snRNA atac clusters",
                  as.pdf = T,
                  filename.prefix = "umap_atac_clusters")

  ##  once we have this
  g4 <- DimPlot(arid1a.multiome, label = F, group.by = 'mammary_gland_cluster_classification') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, panel.border.size = 0,
                    pt.size = 0.9, axis.label = NA, legend.title = '',
                    pal=mammary_cluster_class_palette_list, title = "snRNA mammary cell types", filename.prefix = "umap_mammary_gland_cluster_classification")
  png(filename = file.path(plot.dir,'umap_mammary_gland_cluster_classification_umap_basic.png'), width = 5, height = 2.6, units = 'in', res = 200, bg = 'transparent')
  g4 + glabel$layers
  dev.off()

  g4 <- DimPlot(arid1a.multiome, label = F, group.by = 'mammary_gland_refine_clusters') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    label.size = 10,
                    panel.border.size = 0,
                    pt.size = 0.9,
                    axis.label = NA,
                    legend.title = '',
                    pal=mammary_refine_clusters_palette_list, title = "Mammary cell types",
                    as.pdf = T,
                    filename.prefix = "umap_mammary_gland_refine_clusters")
  png(filename = file.path(plot.dir,'umap_mammary_gland_refine_clusters_basic.png'), width = 5, height = 2.6, units = 'in', res = 200, bg = 'transparent')
  g4 + glabel$layers
  dev.off()

  pdf(file = file.path(plot.dir,'umap_mammary_gland_refine_clusters_basic.pdf'), width = 5, height = 2.6, bg = 'transparent')
  g4 + glabel$layers
  dev.off()

  png(filename = file.path(plot.dir,'umap_basic.png'), width = 15, height = 4, units = 'in', res = 200, bg = 'transparent')
  grid.arrange(g1+glabel$layers, g2, g3, nrow=1)
  dev.off()

  #### UMAP density ####

  g <- UmapPlotDensity(obj = arid1a.multiome, feature = 'Phase',
                       filename.prefix = 'CCPhaseDensityPlots', as.pdf = T,
                       nrow = 1, pt.size = 0.8, pal = rev(pals::ocean.dense(10)))

  g <- UmapPlotDensity(obj = arid1a.multiome,
                       feature = 'condition',
                       filename.prefix = 'ConditionDensityPlots',
                       as.pdf = T, fixed_coord = T, nrow = 1, pt.size = 0.8) #, pal = pals::brewer.purples(10))

  dense.pal <- list(pals::ocean.amp(10),
                    rev(pals::ocean.dense(10)),
                    rev(pals::ocean.algae(10)),
                    rev(pals::ocean.solar(10)))
  names(dense.pal) <- c('Arid1a(-/-)','Arid1a(-/-) E2','Arid1a(+/+)','Arid1a(+/+) E2')

  g <- UmapPlotDensity(obj = arid1a.multiome, feature = 'bioNames',
                       pal = dense.pal, legend.height.factor = .6,
                       filename.prefix = 'BioNamesDensityPlots',
                       as.pdf = T,
                       fixed_coord = T, nrow = 2, pt.size = 0.5)
  pdf(file = file.path(plot.dir,"BioNamesDensityPlots2.pdf"), width=9, height=7)
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  g <- UmapPlotDensity(obj = arid1a.multiome,
                       feature = 'mammary_gland_cluster_classification',
                       #pal = rev(pals::ocean.dense(10)),
                       as.pdf = T,
                       fixed_coord = T, nrow = 2, pt.size = 0.5)

  pdf(file = file.path(plot.dir,"MammaryGlandDensityPlots.pdf"), width=9, height=7)
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  #### Figure 2C-alt UMAP Facet ####
  dense.pal <- list(pals::ocean.amp(10),
                    rev(pals::ocean.dense(10)),
                    rev(pals::ocean.algae(10)),
                    rev(pals::ocean.solar(10)))
  names(dense.pal) <- c('Arid1a(-/-)','Arid1a(-/-) E2','Arid1a(+/+)','Arid1a(+/+) E2')

  g <- DimPlot(arid1a.multiome, group.by = 'bioNames') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    label.size = 10,
                    panel.border.size = 0,
                    pt.size = 0.9,
                    axis.label = "bioNames",
                    legend.title = '',
                    pal=bionames_palette,
                    panel_width_in = 5,
                    #pal = dense.pal,
                    title = "Mammary cell types",
                    as.pdf = T,
                    filename.prefix = "MammaryGlandBioNamesPlots")

  g <- DimPlot(arid1a.multiome, group.by = 'condition') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    label.size = 10,
                    panel.border.size = 0,
                    pt.size = 0.9,
                    axis.label = "condition",
                    legend.title = '',
                    pal=condition_palette,
                    panel_width_in = 5,
                    nrow = 2,
                    #pal = dense.pal,
                    title = "Mammary conditions",
                    as.pdf = T,
                    filename.prefix = "MammaryGlandConditionPlots")

  ## Cell cycle check
  cell_cycle_markers_df = fread(file.path(annot.dir, 'Human_Cell_Cycle_Macosko_et_al.txt'), sep='\t')
  s.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$S[cell_cycle_markers_df$S!=""]) )
  g2m.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$G2.M[cell_cycle_markers_df$G2.M!=""]) )
  m.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$M[cell_cycle_markers_df$M!=""]) )

  ## PCA (same for both assays) ##
  cc <- RunPCA(arid1a.multiome, features = c(s.genes, g2m.genes))
  g1 <- DimPlot(cc, reduction = 'pca', group.by = 'Phase', split.by = 'Phase')
  g2 <- DimPlot(cc, reduction = 'pca', group.by = 'Phase')
  png(file = file.path(plot.dir,'sct_cell_cycle.png'), width = 10, height = 4, bg = 'transparent', units = 'in', res = 150)
  grid.arrange(g1, g2, nrow=1, widths=c(.65,.35))
  dev.off()

  DimHeatmap(arid1a.multiome, dims=1:2, cells = 1000, balanced = T, combine = T, fast = T)
  DimHeatmap(arid1a.multiome, dims=1:2, cells = 1000, balanced = T, combine = T, fast = T)

  ### Cell cycle
  items <- list('S.Score'='S Phase score',
                'G2M.Score'='G2M phase score',
                'Phase'='Cell Cycle phase score')
  palettes <- list( pals::brewer.rdpu(10), pals::brewer.pubu(10), c(pals::brewer.set1(3))) # cols25(25)[c(1,5,3)]))
  g <- lapply(1:3, function(i){
    x <- names(items)[i]
    FeaturePlot(object = arid1a.multiome, features = x, slot = 'data') %>%
      UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, vmin=0, na.color = "gray95",
                      pt.size = 1, axis.label = NA, title = items[[x]], break.unit = 0.1,
                      show.legend = T, scale.color = F, pal = palettes[[i]])
  })
  png(filename = file.path(plot.dir,"UMAP_cell_cycle.png"), width = 13, height = 4.5, units = 'in', bg = 'transparent', res = 200)
  do.call(grid.arrange, c(g, ncol=3))
  dev.off()

  PlotBarVariableInClusters(arid1a.multiome, variable='bioNames', palette = bionames_palette,
                            width_multiplier = 0.7, height_multiplier = 0.7,
                            filename = 'Seurat_cluster_regress_barplot', x.label = "Clusters")
  PlotBarVariableInClusters(arid1a.multiome, variable='condition', palette = condition_palette,
                            width_multiplier = 0.7, height_multiplier = 0.7,
                            filename = 'Seurat_condition_regress_barplot', x.label = "Clusters")
  PlotBarVariableInClusters(arid1a.multiome,
                            by.group = "bioNames",
                            variable='mammary_gland_cluster_classification',
                            palette = mammary_cluster_class_palette,
                            orientation = "vertical",
                            width_multiplier = 1.35, height_multiplier = 0.6,
                            filename = 'Seurat_mammary_gland_cluster_classification_regress_barplot',
                            y.label = "bioNames",
                            x.label = "Percentage of cells")

  gene_list <- rownames(arid1a.multiome@assays$SCT@data)
  test.list <- list("Arid"= grep('Arid', gene_list, value = T),
                    "Fox"=grep('Fox', gene_list, value = T))

  gsct <- lapply(names(test.list), function(x){
    g <- VlnPlot(arid1a.multiome, features = test.list[[x]], assay = 'RNA',pt.size = 0, flip = T,
                 group.by = 'bioNames', stack = T)
    g$theme$legend.position='none'
    g$theme$plot.title=element_text(hjust = 0.5)
    g
  })

  pdf(file = file.path(plot.dir,'Arid_Foxa1_expression_bionames_violin.pdf'), width = 5, height = 12, bg = 'transparent', onefile = T)
  gsct
  dev.off()

  #### Arid1a violin/beeswarm ####
  library(ggsignif)
  g0 <- VlnPlot(arid1a.multiome, features = 'Arid1a', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'bioNames', stack = F, cols = bionames_palette) + xlab("")
  df <- g0$data
  g0 <- ggplot(df, aes(x=ident, y=Arid1a, fill=ident)) +
    geom_violin(trim = T, scale = 'width') +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    scale_fill_manual("", values = bionames_palette) +
    xlab("") + ylab("Expression level") + ggtitle("Arid1a") +
    scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  g0 <- g0 + geom_signif(comparisons = list(c("Arid1a(-/-)", "Arid1a(+/+)"), c("Arid1a(-/-) E2", "Arid1a(+/+) E2")),
                         y_position = c(3.0, 3.5), tip_length = 0, map_signif_level=F)
  pdf(file = file.path(plot.dir,'arid1a_violin.pdf'), width = 5, height = 4, bg = 'transparent')
  g0
  dev.off()

  library(ggbeeswarm)
  g0 <- VlnPlot(arid1a.multiome, features = 'Arid1a', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'bioNames', stack = F, cols = bionames_palette) + xlab("")
  df <- g0$data
  df$ident <- factor(x = df$ident, levels = c("Arid1a(-/-)", "Arid1a(+/+)", "Arid1a(-/-) E2", "Arid1a(+/+) E2"), ordered = T)
  df$ident2 <- df$ident
  #df$ident <- as.character(df$ident)
  g0 <- ggplot(df, aes(x=ident, y=Arid1a, fill=ident, color=ident)) +
    geom_quasirandom(aes(x=ident2), cex = .01, varwidth = .4, width = 1, pch = 1) +
    geom_boxplot(width=0.1, color="black", alpha=0.2) +
    scale_x_discrete(breaks = c("Arid1a(-/-)", "Arid1a(+/+)", "Arid1a(-/-) E2", "Arid1a(+/+) E2"), expand = c(0,0)) +
    scale_fill_manual("", values = bionames_palette) +
    scale_color_manual("", values = bionames_palette) +
    xlab("") + ylab("Expression level") + ggtitle("Arid1a") +
    scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
    theme_classic(14) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16), legend.position = 'none')
  g0 <- g0 + geom_signif(comparisons = list(c("Arid1a(-/-)", "Arid1a(+/+)"), c("Arid1a(-/-) E2", "Arid1a(+/+) E2")),
                         y_position = c(3.0, 3.5), tip_length = 0, map_signif_level=F)
  cairo_pdf(filename = file.path(plot.dir,'arid1a_violin.pdf'), width = 4.5, height = 4.5, bg = 'transparent')
  g0
  dev.off()

  g1 <- VlnPlot(arid1a.multiome, features = 'Foxa1', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'orig.ident', stack = F)
  g2 <- VlnPlot(arid1a.multiome, features = 'Foxa1', assay = 'SCT', slot = 'data',
                pt.size = 0, flip = T, group.by = 'orig.ident', stack = F)
  grid.arrange(g1, g2, ncol=2)

  g0 <- DotPlot(object = arid1a.multiome, assay = 'SCT', features = 'dsRED', group.by = 'orig.ident')
  g1 <- DotPlot(object = arid1a.multiome, assay = 'RNA', features = 'dsRED', group.by = 'orig.ident')
  # by cluster
  g2 <- DotPlot(object = arid1a.multiome, assay = 'RNA', features = 'dsRED')
  grid.arrange(g0, g1, g2, nrow=1)

  #### Markers Violin ####
  #### Violin feature plot ####
  g1 <- lapply(names(markers_list)[c(1:3)], function(x){
    g <- ViolinPlotEnhance(arid1a.multiome,
                           features = markers_list[[x]],
                           assay = 'RNA',
                           pt.size = 0,
                           flip = T,
                           group.by = 'bioNames',
                           stack = T,
                           fill.by = 'ident',
                           title = gsub("_"," ",x),
                           x_axis_text_angle = 90,
                           pal = pals::brewer.purd(20)
    ) + xlab("")
    g # theme(legend.position = 'none')
  })
  g1[[1]] <- g1[[1]] + theme(legend.position = 'none')
  g1[[2]] <- g1[[2]] + theme(legend.position = 'none')
  cairo_pdf(filename = file.path(plot.dir,"Violin_markers_expression_by_sample.pdf"),
            width = 10.5, height = 9, bg = 'transparent', onefile = T)
  do.call(grid.arrange, c(g1, nrow=1))
  dev.off()

  g2 <- lapply(names(markers_list)[c(1:3)], function(x){
    g <- ViolinPlotEnhance(arid1a.multiome,
                           features = markers_list[[x]],
                           assay = 'RNA',
                           pt.size = 0,
                           flip = T,
                           group.by = 'seurat_clusters',
                           stack = T,
                           fill.by = 'ident',
                           title = gsub("_"," ",x),
                           pal = pals::brewer.purd(20)
    ) + xlab("")
    g #+ theme(legend.position = 'none')
  })
  g2[[1]] <- g2[[1]] + theme(legend.position = 'none')
  g2[[2]] <- g2[[2]] + theme(legend.position = 'none')
  cairo_pdf(filename = file.path(plot.dir,"Violin_markers_expression_by_cluster.pdf"),
            width = 17, height = 9, bg = 'transparent', onefile = T)
  do.call(grid.arrange, c(g2, nrow=1))
  dev.off()

  ##### Bach markers #####
  g3 <- lapply(names(bach_marker_list), function(x){
    g <- ViolinPlotEnhance(arid1a.multiome,
                           features =   bach_marker_list[[x]],
                           assay = 'RNA',
                           pt.size = 0,
                           flip = T,
                           group.by = 'seurat_clusters',
                           stack = T,
                           fill.by = 'ident',
                           title = gsub("_"," ",x),
                           pal = pals::brewer.purd(20),
    ) + xlab("")
    g
  })
  for(x in 1:(length(g3)-1)){
    g3[[x]] <- g3[[x]] + theme(legend.position = 'none')
  }
  cairo_pdf(filename = file.path(plot.dir,"Violin_bach_markers_expression_by_cluster.pdf"),
            width = 17, height = 9, bg = 'transparent', onefile = T)
  do.call(grid.arrange, c(g3, nrow=1))
  dev.off()

  maxn <- max(lengths(markers_list)[1:3])
  g1 <- lapply(names(markers_list)[c(1:3)], function(x){
    g <- VlnPlot(arid1a.multiome, features = markers_list[[x]], assay = 'RNA', slot = 'data', pt.size = 0,
                 flip = T, group.by = 'bioNames', stack = T, cols = as.vector(pals::brewer.blues(20)))
    g$labels$title <- gsub("_"," ",x)
    g$theme$legend.position='none'
    g$theme$plot.title=element_text(hjust = 0.5)
    g
  })

  cairo_pdf(filename = file.path(plot.dir,"Violin_markers_expression_RNA_regress.pdf"),
            width = 14, height = 8, bg = 'transparent', onefile = T)
  do.call(grid.arrange, c(g1, nrow=1))
  dev.off()

  g2 <- lapply(names(markers_list)[c(1:3)], function(x){
    g <- VlnPlot(arid1a.multiome, features = markers_list[[x]], assay = 'RNA', pt.size = 0, flip = T,
                 group.by = 'seurat_clusters', stack = T, cols = as.vector(pals::kelly(maxn)))
    g$labels$title <- gsub("_"," ",x)
    g
    #g + scale_fill_viridis(discrete = TRUE)
  })
  #cairo_pdf(filename = file.path(plot.dir,"Violin_markers_expression_cluster.pdf"), width = 26, height = 8, bg = 'transparent', onefile = T)
  png(filename = file.path(plot.dir,"Violin_markers_expression_cluster.png"), width = 22, height = 8, units = 'in', bg = 'transparent', res = 200)
  do.call(grid.arrange, c(g2, nrow=1))
  dev.off()

  maxn <- length(markers_list$`Bach markers`)
  g <- VlnPlot(object = arid1a.multiome, features = markers_list$`Bach markers`,assay = 'RNA', pt.size = 0, flip = T,
          group.by = 'seurat_clusters', stack = T, cols = as.vector(pals::brewer.bugn(n = maxn)))
  g$labels$title <- gsub("_"," ",x)

  DotPlot(object = arid1a.multiome, assay = 'RNA', features = markers_list$`Mammary Luminal progenitor`)

  #### Violin Cell Cycle ####
  f <- c("S.Score","G2M.Score")
  g1 <- ViolinPlotEnhance(arid1a.multiome,
                          features = f,
                          assay = 'RNA',
                          slot = 'data',
                          pt.size = 0,
                          flip = T,
                          group.by = 'seurat_clusters',
                          stack = T,
                          fill.by = 'ident',
                          y_axis_label = "Normalized scores",
                          title = gsub("_"," ","cell cycle scores"),
                          pal = pals::coolwarm(25),
                          x_axis_text_angle = 90,
                          strip.text.size = 14,
                          #pal = rev(pals::kovesi.linear_blue_5_95_c73(20)),
  ) + xlab("")

  pdf(file = file.path(plot.dir,"Cellcycle_violin.pdf"), width = 6, height = 4, bg = "transparent")
  g1
  dev.off()

  #### Markers on UMAP ####
  f <- 'Arid1a'
  (g1 <- FeaturePlot(object = arid1a.multiome, features = f, slot = 'data') %>%
      UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, axis.label = '',
                      legend.height.factor = 0.4, pt.size = .8, outline.size = 0.01,
                      title = paste0(f),
                      scale.color = F,
                      order.random = F,
                      legend.title = "snRNA expression",
                      pal = rocket_palette,
                      filename.prefix = paste0('UMAP_',f), as.pdf = T))

  f <- 'Arid1b'
  (g1 <- FeaturePlot(object = arid1a.multiome, features = f, slot = 'data') %>%
      UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, axis.label = '',
                      legend.height.factor = 0.4, pt.size = .8, outline.size = 0.01,
                      title = paste0(f),
                      scale.color = F,
                      order.random = F,
                      legend.title = "snRNA expression",
                      pal = rocket_palette,
                      filename.prefix = paste0('UMAP_',f), as.pdf = T))

  f <- 'Stc2'
  (g1 <- FeaturePlot(object = arid1a.multiome, features = f, slot = 'data') %>%
      UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, axis.label = '',
                      legend.height.factor = 0.4, pt.size = .8, outline.size = 0.01,
                      title = paste0(f),
                      scale.color = F,
                      order.random = F,
                      legend.title = "snRNA expression",
                      pal = rocket_palette,
                      filename.prefix = paste0('UMAP_',f), as.pdf = T))

  f <- 'Smarca4' # Smarca4 = Brg1 gene
  #DefaultAssay(arid1a.multiome) <- "SCT"
  g1 <- FeaturePlot(object = arid1a.multiome, features = f, slot = 'data') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    label.size = 10,
                    axis.label = '',
                    outline.size = 0.05,
                    legend.height.factor = 0.4,
                    title = paste0(f), scale.color = F, legend.title = "snRNA expression",
                    pal=rocket_palette, filename.prefix = paste0('UMAP_',f), as.pdf = F)

  g <- lapply(markers_list$`Cell Cycle`, function(x){
    FeaturePlot(object = arid1a.multiome, features = x, slot = 'data') %>%
      UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10,
                      pt.size = 1, axis.label = NA, title = x,
                      show.legend = T, scale.color = F, pal = rev(pals::ocean.dense(10)))
  })
  png(filename = file.path(plot.dir,"UMAP_cell_cycle.png"), width = 12, height = 4, units = 'in', bg = 'transparent', res = 200)
  do.call(grid.arrange, c(g, ncol=3))
  dev.off()

  g1 <- DimPlot(arid1a.multiome, label = F, group.by = 'nCount_RNA') + NoLegend()
  g1$data$nCount_RNA <- as.numeric(g1$data$nCount_RNA)
  UmapPlotEnhance(obj = arid1a.multiome, p1 = g1, label.size = 10, axis.label = NA, legend.height.factor = 0.4,
                  break.unit = 2000, pt.size = 1, pal = pals::brewer.purples(10) , #rev(pals::ocean.dense(10)),
                  title = 'No. molecules per cell', scale.color = F, filename.prefix = 'UMAP_count')

  g1 <- DimPlot(arid1a.multiome, label = F, group.by = 'nFeature_RNA') + NoLegend()
  g1$data$nFeature_RNA <- as.numeric(g1$data$nFeature_RNA)
  UmapPlotEnhance(obj = arid1a.multiome, p1 = g1, label.size = 10, axis.label = NA, break.unit = 2000, legend.height.factor = 0.4,
                  pt.size = 1, pal = pals::brewer.purples(10), #rev(pals::ocean.dense(10)),
                  title = 'No. genes per cell', scale.color = F, filename.prefix = 'UMAP_RNAcount')

  lapply(c('Pgr','Csn3',"Ly6d",'Krt7','Xrcc4','Stc2'), function(f){
    FeaturePlot(arid1a.multiome, features = f,
                       pt.size = 0.2, ncol = 1, slot = 'data') %>%
       UmapPlotEnhance(obj = arid1a.multiome, p1 = ., addLabels = F, axis.label = "",
                       pal = rocket_palette, scale.color = F, filename.prefix = paste0("umap_",f), as.pdf = T)
  })

  # Procr set from Bach et al.
  (g2 <- FeaturePlot(arid1a.multiome, features = c('Tns1','Gng11','Procr','Igfbp7','Nrip2','Notch3'),
                     pt.size = 0.5, ncol = 2, slot = 'data'))

  x <- c('Krt8','Foxa1','Krt5','Krt14')
  RidgePlot(arid1a.multiome, features = x, ncol = 3)
  g1 <- FeaturePlot(arid1a.multiome, features = x[c(1,2)], pt.size = 0.2, ncol = 1, slot = 'data', blend=T)

  g1 <- lapply(x, function(i){
    g1 <- FeaturePlot(object = arid1a.multiome, features = i, slot = 'data') %>%
      UmapPlotEnhance(obj = arid1a.multiome,
                      p1 = .,
                      label.size = 10,
                      axis.label = '',
                      title = i, scale.color = F,
                      pal=rocket_palette,
                      as.pdf = T,
                      filename.prefix = paste0("UMAP_",i))
  })

  features <- c("Mat2b","Kdm1a","Prdm2","Prdm16")
  g <- lapply(features, function(f){
    FeaturePlot(object = arid1a.multiome, features = f, slot = 'data') %>%
      UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, axis.label = '',
                    legend.height.factor = 0.4,
                    title = paste0(f), scale.color = F, legend.title = "snRNA expression",
                    pal=pals::brewer.purd(20))
  })

  png(filename = file.path(plot.dir,"Amaia_genes.png"), width = 10, height = 7, units = 'in', res = 200)
  do.call(grid.arrange, c(g, nrow=2))
  dev.off()

  #### DE Markers (1-vs-all) ####
  # https://github.com/immunogenomics/presto
  library(presto)

  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description"))
  gene_annot <- gene_annot[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]

  # DE bionames vs others
  DEsample_v_others_wilcox <- FindAllMarkers(arid1a.multiome, group.by = 'bioNames', test.use = "wilcox", minmin.pct = 0.01, logfc.threshold = 0)
  DEsample_v_others_roc <- FindAllMarkers(arid1a.multiome, group.by = 'bioNames', test.use = "roc", minmin.pct = 0.01, logfc.threshold = 0)
  ##DEsample_v_others <- wilcoxauc(arid1a.multiome, group_by = 'bioNames', seurat_assay = 'RNA', assay = 'data')
  DEsample_v_others <- setDT(DEsample_v_others)
  setnames(DEsample_v_others,
           old=c("cluster","pct.1","pct.2"),
           new = c("cluster1","percent_cluster1","percent_cluster2"))
  DEsample_v_others[,cluster2 := 'others']
  DEsample_v_others <- DEsample_v_others[order(p_val_adj, decreasing = T)]
  DEsample_v_others[, head(.SD, 10), by=cluster1]
  ## Add gene annotation and sort outputs
  # DE sample
  DEsample_v_others <- merge(DEsample_v_others, gene_annot, by='gene', all.x=T)
  DEsample_v_others[,rank:= avg_logFC*log(statistic)]
  DEsample_v_others <- DEsample_v_others[order(rank, decreasing = T)]

  plot(DEsample_v_others[cluster1=='Arid1a(-/-) E2' ,rank],
       DEsample_v_others[cluster1=='Arid1a(-/-) E2' ,auc])

  # DE cluster vs others
  # DE bionames vs others
  DEcluster_v_others_wilcox <- FindAllMarkers(arid1a.multiome,
                                              group.by = 'seurat_clusters',
                                              test.use = "wilcox",
                                              minmin.pct = 0.01, logfc.threshold = 0)

  DEcluster_v_others_wilcox <- setDT(DEcluster_v_others_wilcox)
  DEcluster_v_others_wilcox[,cluster2 := 'others']
  DEcluster_v_others_wilcox <- DEcluster_v_others_wilcox[order(p_val_adj, decreasing = T)]
  ## Add gene annotation and sort outputs
  # DE sample
  DEcluster_v_others <- merge(DEcluster_v_others, gene_annot, by='gene', all.x=T)
  DEcluster_v_others_wilcox <- DEcluster_v_others_wilcox[order(p_val_adj, decreasing = F)]
  DEcluster_v_others_wilcox[, head(.SD, 10), by=cluster1]

  #### DE Markers KO vs WT (all) ####
  table(arid1a.multiome@meta.data$condition)
  #bionames
  DEsample_KO_vs_WT <- wilcoxauc(arid1a.multiome, group_by = 'condition', seurat_assay = 'RNA', assay = 'data')
  setDT(DEsample_KO_vs_WT)
  setnames(DEsample_KO_vs_WT,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEsample_KO_vs_WT[cluster1=='KO', cluster2 := 'WT']
  DEsample_KO_vs_WT[cluster1=='WT', cluster2 := 'KO']
  DEsample_KO_vs_WT <- merge(DEsample_KO_vs_WT, gene_annot, by='gene', all.x=T)
  DEsample_KO_vs_WT[,rank:= avg_logFC*log(statistic)]
  DEsample_KO_vs_WT <- DEsample_KO_vs_WT[order(rank, decreasing = T)]
  DEsample_KO_vs_WT[gene=="Ly6d"]
  DEsample_KO_vs_WT[gene=="Ctcf"]

  #### DE compare mammary_gland_cluster_classification
  require(presto)
  table(arid1a.multiome$mammary_gland_cluster_classification, arid1a.multiome$seurat_clusters)
  DE_mammary_gland_cluster_classification_vs_others <- wilcoxauc(arid1a.multiome, group_by = 'mammary_gland_cluster_classification', seurat_assay = 'RNA', assay = 'data')
  setDT(DE_mammary_gland_cluster_classification_vs_others)
  setnames(DE_mammary_gland_cluster_classification_vs_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DE_mammary_gland_cluster_classification_vs_others[, cluster2 := 'others']
  ## Add gene annotation and sort outputs
  DE_mammary_gland_cluster_classification_vs_others <- merge(DE_mammary_gland_cluster_classification_vs_others, gene_annot, by='gene', all.x=T)
  DE_mammary_gland_cluster_classification_vs_others[,rank:= avg_logFC*log(statistic)]
  setorderv(x = DE_mammary_gland_cluster_classification_vs_others, cols = c('cluster1', 'rank'), order = c(1, -1))
  DE_mammary_gland_cluster_classification_vs_others

  #### DE compare mammary_gland_refine_clusters
  require(presto)
  table(arid1a.multiome$mammary_gland_refine_clusters, arid1a.multiome$seurat_clusters)
  DE_mammary_gland_refine_cluster_vs_others <- wilcoxauc(arid1a.multiome, group_by = 'mammary_gland_refine_clusters', seurat_assay = 'RNA', assay = 'data')
  setDT(DE_mammary_gland_refine_cluster_vs_others)
  setnames(DE_mammary_gland_refine_cluster_vs_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DE_mammary_gland_refine_cluster_vs_others[, cluster2 := 'others']
  ## Add gene annotation and sort outputs
  DE_mammary_gland_refine_cluster_vs_others <- merge(DE_mammary_gland_refine_cluster_vs_others, gene_annot, by='gene', all.x=T)
  DE_mammary_gland_refine_cluster_vs_others[,rank:= avg_logFC*log(statistic)]
  setorderv(x = DE_mammary_gland_refine_cluster_vs_others, cols = c('cluster1', 'rank'), order = c(1, -1))
  DE_mammary_gland_refine_cluster_vs_others

  #### DE compare WT L-alv L-hormonal clusters ####
  require(presto)
  table(arid1a.multiome$mammary_gland_cluster_classification, arid1a.multiome$seurat_clusters)
  # 1,6 vs 4,5
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(1,6,4,5) )
  table(ss_arid$mammary_gland_cluster_classification, ss_arid$seurat_clusters)
  DE_lalv_vs_lhor <- wilcoxauc(ss_arid, group_by = 'mammary_gland_cluster_classification', seurat_assay = 'RNA', assay = 'data')
  setDT(DE_lalv_vs_lhor)
  setnames(DE_lalv_vs_lhor,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DE_lalv_vs_lhor[cluster1=='L-Alv', cluster2 := 'L-Hor']
  DE_lalv_vs_lhor[cluster1=='L-Hor', cluster2 := 'L-Alv']
  ## Add gene annotation and sort outputs
  DE_lalv_vs_lhor <- merge(DE_lalv_vs_lhor, gene_annot, by='gene', all.x=T)
  DE_lalv_vs_lhor[,rank:= avg_logFC*log(statistic)]
  setorderv(x = DE_lalv_vs_lhor, cols = c('cluster1', 'rank'), order = c(1, -1))
  DE_lalv_vs_lhor

  #### DE compare WT E2 vs not-E2, 4 vs 5 clusters in L-Hor ####
  require(presto)
  table(arid1a.multiome$treatment, arid1a.multiome$seurat_clusters)
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(4,5) & bioNames %in% c("Arid1a(+/+)", "Arid1a(+/+) E2") )
  table(ss_arid$treatment, ss_arid$seurat_clusters)
  table(ss_arid$bioNames)
  DE_WT_E2_vs_WT <- wilcoxauc(ss_arid, group_by = 'treatment', seurat_assay = 'RNA', assay = 'data')
  setDT(DE_WT_E2_vs_WT)
  setnames(DE_WT_E2_vs_WT,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DE_WT_E2_vs_WT[cluster1=="E2", cluster1:="WT_E2"]
  DE_WT_E2_vs_WT[cluster1=="None", cluster1:="WT"]

  DE_WT_E2_vs_WT[cluster1=='WT_E2', cluster2 := 'WT']
  DE_WT_E2_vs_WT[cluster1=='WT', cluster2 := 'WT_E2']
  ## Add gene annotation and sort outputs
  DE_WT_E2_vs_WT <- merge(DE_WT_E2_vs_WT, gene_annot, by='gene', all.x=T)
  DE_WT_E2_vs_WT[,rank:= avg_logFC*log(statistic)]
  setorderv(x = DE_WT_E2_vs_WT, cols = c('cluster1', 'rank'), order = c(1, -1))
  DE_WT_E2_vs_WT

  #### DE compare KO E2 vs not-E2, 2 vs 3 clusters ####
  require(presto)
  table(arid1a.multiome$treatment, arid1a.multiome$seurat_clusters)
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(2,3) & bioNames %in% c("Arid1a(-/-)", "Arid1a(-/-) E2") )
  table(ss_arid$treatment, ss_arid$seurat_clusters)
  table(ss_arid$bioNames)
  DE_KO_E2_vs_KO <- wilcoxauc(ss_arid, group_by = 'treatment', seurat_assay = 'RNA', assay = 'data')
  setDT(DE_KO_E2_vs_KO)
  setnames(DE_KO_E2_vs_KO,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DE_KO_E2_vs_KO[cluster1=="E2", cluster1:="KO_E2"]
  DE_KO_E2_vs_KO[cluster1=="None", cluster1:="KO"]

  DE_KO_E2_vs_KO[cluster1=='KO_E2', cluster2 := 'KO']
  DE_KO_E2_vs_KO[cluster1=='KO', cluster2 := 'KO_E2']
  ## Add gene annotation and sort outputs
  DE_KO_E2_vs_KO <- merge(DE_KO_E2_vs_KO, gene_annot, by='gene', all.x=T)
  DE_KO_E2_vs_KO[,rank:= avg_logFC*log(statistic)]
  setorderv(x = DE_KO_E2_vs_KO, cols = c('cluster1', 'rank'), order = c(1, -1))
  DE_KO_E2_vs_KO

  #### DE compare L-Hor KO vs WT 2,3 vs 4,5 clusters ####
  require(presto)
  table(arid1a.multiome$treatment, arid1a.multiome$seurat_clusters)
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(2,3,4,5) )
  table(ss_arid$condition, ss_arid$seurat_clusters)
  DE_KO_Lhor_vs_WT_Lhor <- wilcoxauc(ss_arid, group_by = 'condition', seurat_assay = 'RNA', assay = 'data')
  setDT(DE_KO_Lhor_vs_WT_Lhor)
  setnames(DE_KO_Lhor_vs_WT_Lhor,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DE_KO_Lhor_vs_WT_Lhor[cluster1=='KO', cluster2 := 'WT']
  DE_KO_Lhor_vs_WT_Lhor[cluster1=='WT', cluster2 := 'KO']
  ## Add gene annotation and sort outputs
  DE_KO_Lhor_vs_WT_Lhor <- merge(DE_KO_Lhor_vs_WT_Lhor, gene_annot, by='gene', all.x=T)
  DE_KO_Lhor_vs_WT_Lhor[,rank:= avg_logFC*log(statistic)]
  setorderv(x = DE_KO_Lhor_vs_WT_Lhor, cols = c('cluster1', 'rank'), order = c(1, -1))
  DE_KO_Lhor_vs_WT_Lhor

  groups_use <- c(2,4)
  DE_2_4_cluster_compare <- wilcoxauc(arid1a.multiome, groups_use = groups_use, group_by = 'seurat_clusters', seurat_assay = 'RNA', assay = 'data')
  setDT(DE_2_4_cluster_compare)
  setnames(DE_2_4_cluster_compare,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DE_2_4_cluster_compare[cluster1==groups_use[1], cluster2 := groups_use[2]]
  DE_2_4_cluster_compare[cluster1==groups_use[2], cluster2 := groups_use[1]]
  ## Add gene annotation and sort outputs
  DE_2_4_cluster_compare <- merge(DE_2_4_cluster_compare, gene_annot, by='gene', all.x=T)
  DE_2_4_cluster_compare[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  setorderv(x = DE_2_4_cluster_compare, cols = c('cluster1', 'rank'), order = c(1, -1))
  DE_2_4_cluster_compare

  #### DE spreadsheet ####
  require(writexl)
  require(readxl)
  # write
  writexl::write_xlsx(x = DE_WT_E2_vs_WT,
                      path = file.path(spreadsheets.dir,'scRNA_WTE2_vs_WT_4_5_wilcox_de.xlsx'))

  writexl::write_xlsx(x = DE_KO_E2_vs_KO,
                      path = file.path(spreadsheets.dir,'scRNA_KOE2_vs_KO_3_2_wilcox_de.xlsx'))

  writexl::write_xlsx(x = DEsample_v_others,
                      path = file.path(spreadsheets.dir,'scRNA_sample_vs_others_wilcox_de.xlsx'))

  writexl::write_xlsx(x = DEcluster_v_others_wilcox,
                      path = file.path(spreadsheets.dir,'scRNA_cluster_vs_others_wilcox_de.xlsx'))

  writexl::write_xlsx(x = DEsample_KO_vs_WT,
                      path = file.path(spreadsheets.dir,'scRNA_sample_KO_vs_WT_wilcox_de.xlsx'))

  writexl::write_xlsx(x = DE_lalv_vs_lhor,
                      path = file.path(spreadsheets.dir,'scRNA_lalv_vs_lhor_wilcox_de.xlsx'))

  writexl::write_xlsx(x = DE_mammary_gland_cluster_classification_vs_others,
                      path = file.path(spreadsheets.dir,'scRNA_mammary_gland_cluster_classification_vs_others_wilcox_de.xlsx'))

  writexl::write_xlsx(x = DE_mammary_gland_refine_cluster_vs_others, path = file.path(spreadsheets.dir,'scRNA_mammary_gland_refine_cluster_vs_others_wilcox_de.xlsx'))

  writexl::write_xlsx(x = DE_KO_Lhor_vs_WT_Lhor, path = file.path(spreadsheets.dir,'scRNA_KO_Lhor_vs_WT_Lhor_wilcox_de.xlsx'))

  # read
  DE_WT_E2_vs_WT <- setDT(read_xlsx(file.path(spreadsheets.dir,'scRNA_WTE2_vs_WT_4_5_wilcox_de.xlsx')))
  DE_KO_E2_vs_KO <- setDT(read_xlsx(file.path(spreadsheets.dir,'scRNA_KOE2_vs_KO_3_2_wilcox_de.xlsx')))

  DEsample_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_sample_vs_others_wilcox_de.xlsx')))
  DEcluster_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_cluster_vs_others_wilcox_de.xlsx')))
  DEsample_KO_vs_WT <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_sample_KO_vs_WT_wilcox_de.xlsx')))
  DEsample_KO_None_vs_WT_None <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_sample_KO_None_vs_WT_None_wilcox_de.xlsx')))
  DEsample_KO_E2_vs_WT_E2 <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_sample_KO_E2_vs_WT_E2_wilcox_de.xlsx')))
  DE_lalv_vs_lhor <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_lalv_vs_lhor_wilcox_de.xlsx')))
  DE_mammary_gland_cluster_classification_vs_others <- setDT(read_xlsx(file.path(spreadsheets.dir,'scRNA_mammary_gland_cluster_classification_vs_others_wilcox_de.xlsx')))
  DE_mammary_gland_refine_cluster_vs_others <- setDT(read_xlsx(file.path(spreadsheets.dir,'scRNA_mammary_gland_refine_cluster_vs_others_wilcox_de.xlsx')))

  DE_KO_Lhor_vs_WT_Lhor <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'scRNA_KO_Lhor_vs_WT_Lhor_wilcox_de.xlsx')))

  saveRDS(arid1a.multiome, file = file.path(robj.dir, 'arid1a_sctransform_multiome_regress.rds'))

  arid1a.multiome <- readRDS(file = file.path(robj.dir, 'arid1a_sctransform_multiome_regress.rds'))

  #### Heatmaps ####
  x <- markers_list[c(1:3)]
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = unlist(x),
                         group.by = 'seurat_clusters',
                         assay = 'RNA',
                         slot = 'data',
                         show.group.dendrogram = T,
                         show.feature.names = T,
                         group.fontsize = 14, feature.fontsize = 12,
                         group.colors = seurat_cluster_palette)
  png(filename = file.path(plot.dir,'Heatmap_Luminal_basal_cluster.png'),
      width = 12, height = 12, units = 'in', res = 200, bg = 'transparent')
  draw(ht,
       column_title = 'Luminal/Basal markers',
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()

  #### Saeki markers Heatmap ####
  s <- readxl::excel_sheets(file.path(data.dir,'Mammary_Significant_GeneMarkers_by_Cluster_Saeki_et_al.xlsx'))
  geneLists <- lapply(s, function(i){
    wb <- readxl::read_xlsx(file.path(data.dir,'Mammary_Significant_GeneMarkers_by_Cluster_Saeki_et_al.xlsx'),
                            sheet = i, col_names = T)
    setDT(wb)
  }) %>% set_names(s)

  x <- sapply(s, function(i){
    head(geneLists[[i]], 10)
  })
  unlist(x)

  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = unlist(x),
                         group.by = 'seurat_clusters',
                         assay = 'RNA',
                         layer = 'data',
                         clustering.within = F,
                         cluster_group_slices = T,
                         cluster.features = F,
                         show.cell.dendrogram = F,
                         feature.fontsize = 9,
                         group.colors =  seurat_cluster_palette)

  png(filename = file.path(plot.dir,'Heatmap_Mammary_Markers_cluster.png'),
      width = 12, height = 10, units = 'in', res = 200, bg = 'transparent')
  draw(ht,
       column_title = 'Mammary Markers (top 10 per category)',
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()

  #### Samples heatmap ####
  rem <- grep("^Gm[0-9]+", rownames(arid1a.multiome), value=T)
  rem <- c(rem, grep("^[0-9]+",rownames(arid1a.multiome), value=T))
  rem <- c(rem, grep("^mt-",rownames(arid1a.multiome), value=T, ignore.case = T))
  ribosomal_genes = fread(file.path(annot.dir, 'mouse_ribosomal_genes.txt'), header = F)
  ribosomal_genes <- paste(paste0("^",unname(unlist(ribosomal_genes))), collapse = "|")
  rem <- c(rem, grep(ribosomal_genes, rownames(arid1a.multiome), value=T))

  DEsample_v_others[p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, .N, by=cluster1]
  x <- DEsample_v_others[p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, gene] %>% unique()
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x,
                         show.feature.names = F,
                         group.by = 'bioNames',
                         assay = 'RNA',
                         slot = 'data',
                         clustering.within = F,
                         heatmap.orientation = "horizontal",
                         group.colors = bionames_palette)

  f <- file.path(plot.dir,'DEsample_v_others_heatmap_pval_0.05_avglog_0.5.png')

  require(magick)
  message("Plotting to:",f)
  res <- 250
  fig <- image_graph(width = 8*res, height = 3.5*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("sample vs others\n(FDR ≤ 0.05, | Log2FC ≥0.05)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  topn <-  30
  x <- DEsample_v_others[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x,
                         heatmap.orientation = 'horizontal',
                         group.by = 'bioNames',
                         assay = 'RNA',
                         slot = 'data',
                         feature.fontsize = 6,
                         group.colors = bionames_palette)

  f <- file.path(plot.dir,paste0('DEsample_v_others_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 10*res, height = 3*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"sample vs others"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')
  dev.off()

  # KO vs WT all
  topn <- 50
  x <- DEsample_KO_vs_WT[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x,
                         group.by = 'condition',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         show.group.dendrogram = T,
                         #show.cell.dendrogram = T,
                         group.colors = condition_palette)
  f <- file.path(plot.dir,paste0('DE_KO_v_WT_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 12*res, height = 2.5*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"KO vs WT"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')
  dev.off()

  #### Cluster heatmap / dotplot ####
  DEcluster_v_others[auc >= 0.7 | auc <= 0.3, .N, by=cluster1][order(N, decreasing = T)]
  x <- DEcluster_v_others[auc >= 0.7 | auc <= 0.3, unique(gene)]
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x,
                         group.by = 'seurat_clusters',
                         assay = 'RNA',
                         slot = 'data',
                         show.dendrogram = T,
                         show.features.names = F,
                         group.font.size = 12,
                         disp.heatmap.vertically = F,
                         group.colors = seurat_cluster_palette)
  png(filename = file.path(plot.dir,'Heatmap_Cluster_vs_others_DE.png'),
      width = 8, height = 7, units = 'in', res = 200, bg = 'transparent')
  draw(ht,
       column_title = paste0("Cluster vs Others n=",length(x),"\n(auc ≤ 0.3 | auc ≥ 0.7)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1),
       heatmap_legend_side = "right")
  dev.off()

  topn <- 15
  #x <- DEcluster_v_others[!gene %in% rem][order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, .(gene=head(gene, topn)), by=cluster1]
  x <- DEcluster_v_others[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, .(gene=head(gene, topn)), by=cluster1]
  x <- unique(x$gene)
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x,
                         heatmap.orientation = 'vertical',
                         group.by = 'seurat_clusters',
                         assay = 'RNA',
                         slot = 'data',
                         feature.fontsize = 6,
                         clustering.within = F,
                         group.colors = seurat_cluster_palette)

  f <- file.path(plot.dir,paste0('DEcluster_v_others_heatmap_features_vertical_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 4.5*res, height = 11*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"cluster vs others"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  ##### Dotplot #####
  f <- DEcluster_v_others_wilcox[avg_log2FC > 0, head(.SD, 10), by=cluster1]
  uclust <- unique(f$cluster1)
  demarkersL <- lapply(uclust, function(i){
    f[cluster1==i, gene]
  }) %>% set_names(uclust)

  dp <- MyDotPlot(object = arid1a.multiome,
                  features = demarkersL,
                  plot.orientation = "h",
                  legend.box = "vertical",
                  dot.scale = 5,
                  scale = T, scale.min = 0, scale.max = 1,
                  cols = pals::brewer.reds(20),
                  #group.by = "saeki_mammary_classification",
                  group.by = "mammary_gland_refine_clusters",
                  cluster.idents = T)
  pdf(file = file.path(plot.dir,"DEdotplot_by_cluster.pdf"), width = 16, height = 3)
  dp
  dev.off()

  #### Heatmap WT L-Alv vs L-hormone in WT ####
  topn <- 50
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(1,6,4,5) )
  x <- DE_lalv_vs_lhor[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = ss_arid,
                         features = x,
                         group.by = 'mammary_gland_cluster_classification',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         show.group.dendrogram = T,
                         #show.cell.dendrogram = T,
                         group.colors = mammary_cluster_class_palette)
  f <- file.path(plot.dir,paste0('DE_lalv_vs_lhor_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 12*res, height = 2.1*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"L-Alv vs L-Hor in WT"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')
  dev.off()

  #### Upregulated in E2 WT vs WT and also not upregulated in E2 K vs KO ####
  # I imagine 2 scenarios:
  # 1. Gene is significantly up or down regulated in E2 WT vs WT, but it is not significantly different in E2 KO vs KO.
  # 2. Opposite regulation:
  #    2.1 Gene is significantly upregulated in E2 WT vs WT, but it is downregulated in E2 KO vs KO.
  #    2.2 Gene is significantly downregulated in E2 WT vs WT, but it is upregulated in E2 KO vs KO.
  library(readxl)
  DE_WT_E2_vs_WT <- setDT(read_xlsx(file.path(spreadsheets.dir,'scRNA_WTE2_vs_WT_4_5_wilcox_de.xlsx')))
  DE_KO_E2_vs_KO <- setDT(read_xlsx(file.path(spreadsheets.dir,'scRNA_KOE2_vs_KO_3_2_wilcox_de.xlsx')))

  # 1. scenario
  # Gene is significantly up or down regulated in E2 WT vs WT, but it is not significantly different in E2 KO vs KO.
  #X1 <- DE_WT_E2_vs_WT[cluster1=="WT_E2" & p_val_adj < 0.05 & abs(avg_logFC) >= 0.1, gene]
  X1 <- DE_WT_E2_vs_WT[cluster1=="WT_E2" & p_val_adj < 0.05, gene]
  X2 <- DE_KO_E2_vs_KO[cluster1=="KO_E2" & p_val_adj < 0.05, gene]
  # Get KO not in WT
  #X0 <- setdiff(X2, X1)
  # Get WT not in KO
  X0 <- setdiff(X1, X2)
  # Get KO not in WT
  X3 <- setdiff(X2, X1)

  DE_WT_E2_vs_WT[]

  length(X2)
  length(X1)
  length(X0)
  length(X3)

  source(file.path(lib.dir, 'venn_helpers.R'))
  plotVenn3(list.items = list("WTE2vWT_no_KOE2vKO"=X0,
                              "WTE2vWT"=X1,
                              "KOE2vKO_no_WTE2vWT"=X3))

  pdf(file = file.path(plot.dir,"Venn_diagram_KOE2vKO_WTE2vWT.pdf"), width = 5, height = 5, bg = "transparent")
  plotVenn2(list.items = list("WTE2vWT"=X1,
                              "KOE2vKO"=X2),
            fill.colors = c("dodgerblue2","red"),
            cat.colors = c("dodgerblue2","red"))
  dev.off()

  DE_WT_E2_vs_WT[gene %in% X0][order(gene, rank)]
  DE_KO_E2_vs_KO[gene %in% X0][order(gene, rank)]

  DE_WT_E2_vs_WT[gene %in% "Foxa1"]
  DE_KO_E2_vs_KO[gene %in% "Foxa1"]

  # Genes possible for the heatmap
  writexl::write_xlsx(x = DE_WT_E2_vs_WT[gene %in% X0], path = file.path(spreadsheets.dir,"DE_WTE2_VS_WT_NOT_KOE2_VS_KO xlsx"))
  DEset <- read_xlsx(file.path(spreadsheets.dir,"DE_WTE2_VS_WT_NOT_KOE2_VS_KO xlsx"))

  topn <- 100
  ss_arid <- subset(arid1a.multiome, (seurat_clusters %in% c(4,5) & bioNames %in% c("Arid1a(+/+)","Arid1a(+/+) E2")) |
                      (seurat_clusters %in% c(2,3) & bioNames %in% c("Arid1a(-/-)","Arid1a(-/-) E2")))
  x <- head(X0, topn)

  #### Heatmap ####
  ht <- DoHeatmapEnhance(object = ss_arid,
                         features = x,
                         group.by = 'bioNames',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 11,
                         show.feature.names = T,
                         clustering.within = T,
                         show.group.dendrogram = T,
                         title = paste0("Top ",length(x)," genes differential Arid1a WT+E2 vs Arid1a WT but not in KO+E2 vs KO"),
                         group.colors = bionames_palette)
  f <- file.path(plot.dir,paste0('DE_WTE2_not_in_KOE2_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 16*res, height = 4*res, res = res, bg = 'transparent')
  draw(ht)
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  pdf(file = file.path(plot.dir,paste0('DE_WTE2_not_in_KOE2_heatmap_features_top',topn,'.pdf')),
      width = 17, height = 4, bg = "transparent")
  draw(ht)
  dev.off()

  #### Dotplot ####
  ##### Dotplot #####
  #f <- DE_WT_E2_vs_WT[gene %in% head(X0, topn), list(gene, cluster1)]
  f <- DE_WT_E2_vs_WT[gene %in% head(X0, topn)]
  uclust <- unique(f$cluster1)
  demarkersL <- lapply(uclust, function(i){
    f[cluster1==i & avg_logFC>0, gene]
  }) %>% set_names(uclust)

  dp <- MyDotPlot(object = ss_arid,
                  features = demarkersL,
                  plot.orientation = "h",
                  legend.box = "vertical",
                  dot.scale = 5,
                  font.base.size = 16,
                  scale = T, scale.min = 0, scale.max = 1,
                  cols = pals::brewer.reds(20),
                  #group.by = "saeki_mammary_classification",
                  group.by = "bioNames",
                  cluster.idents = T)
  pdf(file = file.path(plot.dir,"DEdotplot_WTE2_vs_WT_not_KO_by_biotype.pdf"), width = 19, height = 2.6)
  dp
  dev.off()

  # 2. scenario (confusing)
  # Opposite regulation:
  # 2.1 Gene is significantly upregulated in E2 WT vs WT, but it is downregulated in E2 KO vs KO.
  # 2.2 Gene is significantly downregulated in E2 WT vs WT, but it is upregulated in E2 KO vs KO.
  X1 <- DE_WT_E2_vs_WT[cluster1=="WT_E2" & p_val_adj < 0.05 & avg_logFC > 0, gene]
  X2 <- DE_KO_E2_vs_KO[cluster1=="KO_E2" & p_val_adj < 0.05 & avg_logFC < 0, gene]
  X0 <- intersect(X1, X2)

  length(unique(X0))
  DE_WT_E2_vs_WT[gene %in% X0][order(gene, rank)]
  DE_KO_E2_vs_KO[gene %in% X0][order(gene, rank)]

  DE_WT_E2_vs_WT[gene %in% "Kif5b"]
  DE_KO_E2_vs_KO[gene %in% "Kif5b"]

  topn <- 50
  ss_arid <- subset(arid1a.multiome, (seurat_clusters %in% c(4,5) & bioNames %in% c("Arid1a(+/+)","Arid1a(+/+) E2")) |
                      (seurat_clusters %in% c(2,3) & bioNames %in% c("Arid1a(-/-)","Arid1a(-/-) E2")))
  x <- head(X0, topn)
  ht <- DoHeatmapEnhance(object = ss_arid,
                         features = x,
                         group.by = 'bioNames',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 11,
                         show.feature.names = T,
                         clustering.within = F,
                         show.group.dendrogram = T,
                         title = paste0("Top ",length(x)," genes UP E2WT vs WT and DOWN DE E2KO vs KO"),
                         group.colors = bionames_palette)
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 16*res, height = 4*res, res = res, bg = 'transparent')
  draw(ht)
  #column_title = paste("TOP",topn,"E2 vs None in WT L-Hor (clusters 4, 5)"),
  #column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  #image_write(image = fig, path = f, format = 'png')

  #### Heatmap DE_KO_E2_vs_None in L-Hor ####
  topn <- 50
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(2,3) )
  x <- DE_KO_E2_vs_None[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = ss_arid,
                         features = x,
                         group.by = 'treatment',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         show.group.dendrogram = T,
                         #show.cell.dendrogram = T,
                         group.colors = treatment_palette)
  f <- file.path(plot.dir,paste0('DE_KO_E2_vs_None_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 8*res, height = 2.1*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"E2 vs None in KO L-Hor (clusters 2, 3)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  #### Heatmap DE_WT_E2_vs_WT in L-Hor ####
  DE_WT_E2_vs_WT <- setDT(read_xlsx(file.path(spreadsheets.dir,'scRNA_WTE2_vs_WT_4_5_wilcox_de.xlsx')))
  DE_KO_E2_vs_KO <- setDT(read_xlsx(file.path(spreadsheets.dir,'scRNA_KOE2_vs_KO_3_2_wilcox_de.xlsx')))

  topn <- 50
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(4,5) )
  x <- DE_WT_E2_vs_None[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = ss_arid,
                         features = x,
                         group.by = 'treatment',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         show.group.dendrogram = T,
                         #show.cell.dendrogram = T,
                         group.colors = treatment_palette)
  f <- file.path(plot.dir,paste0('DE_WT_E2_vs_None_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 12*res, height = 2.1*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"E2 vs None in WT L-Hor (clusters 4, 5)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  #### Heatmap DE_KO_E2_vs_None in L-Hor ####
  topn <- 50
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(2,3) )
  x <- DE_KO_E2_vs_None[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = ss_arid,
                         features = x,
                         group.by = 'treatment',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         show.group.dendrogram = T,
                         #show.cell.dendrogram = T,
                         group.colors = treatment_palette)
  f <- file.path(plot.dir,paste0('DE_KO_E2_vs_None_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 8*res, height = 2.1*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"E2 vs None in KO L-Hor (clusters 2, 3)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  #### Heatmap mammary_gland_cluster_classification vs others ####
  #DE_mammary_gland_cluster_classification_vs_others
  topn <- 15
  x <- DE_mammary_gland_cluster_classification_vs_others[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x,
                         group.by = 'mammary_gland_cluster_classification',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         group.fontsize = 12,
                         show.group.dendrogram = T,
                         #show.cell.dendrogram = T,
                         group.colors = mammary_cluster_class_palette_list)
  f <- file.path(plot.dir,paste0('mammary_gland_cluster_vs_others_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 10*res, height = 3.1*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"mammary gland clusters"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  pdf(file = file.path(plot.dir,paste0('mammary_gland_cluster_vs_others_heatmap_features_top',topn,'.pdf')), width = 10, height = 3.1)
  draw(ht,
       column_title = paste("TOP",topn,"mammary gland clusters"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()

  #### Heatmap mammary_gland_refine_cluster vs others ####
  topn <- 15
  x <- DE_mammary_gland_refine_cluster_vs_others[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x,
                         group.by = 'mammary_gland_refine_clusters',
                         assay = 'RNA',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         group.fontsize = 12,
                         show.group.dendrogram = T,
                         #show.cell.dendrogram = T,
                         group.colors = mammary_refine_clusters_palette)
  f <- file.path(plot.dir,paste0('mammary_gland_refine_cluster_vs_others_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 11*res, height = 3.1*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"mammary gland cell types"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  pdf(file = file.path(plot.dir,paste0('mammary_gland_refine_cluster_vs_others_heatmap_features_top',topn,'.pdf')), width = 11, height = 3.1)
  draw(ht,
       column_title = paste("TOP",topn,"mammary gland cell types"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()

  ##### DE Dotplot mammary_gland_refine_cluster #####
  topn <- 15
  x <- DE_mammary_gland_refine_cluster_vs_others[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, .(gene=head(gene, topn)), by=cluster1]
  x <- split(x, by="cluster1", keep.by=F)
  feature.list <- lapply(x, function(y) { unlist(y, use.names = F)})
  gd <- MyDotPlot(object = arid1a.multiome,
                  features = feature.list,
                  plot.orientation = "h",
                  legend.box = "vertical",
                  group.by = "mammary_gland_refine_clusters",
                  dot.scale = 5,
                  scale = T,
                  scale.min = 0, scale.max = 1,
                  cols = pals::brewer.reds(20),
                  cluster.idents = T)
  pdf(file = file.path(plot.dir,"chromatin_modifiers_by_celltype.pdf"),
      width = 8, height = 3)
  gd + theme(strip.text = element_text(angle = 45, vjust=0, hjust=0))
  dev.off()

  f <- file.path(plot.dir,paste0('mammary_gland_refine_cluster_vs_others_dotplot_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 11*res, height = 3.1*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"mammary gland clusters"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  pdf(file = file.path(plot.dir,paste0('mammary_gland_refine_cluster_vs_others_heatmap_features_top',topn,'.pdf')), width = 11, height = 3.1)
  draw(ht,
       column_title = paste("TOP",topn,"mammary gland clusters"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()

  intersect(markers_list$HALLMARK_ESTROGEN_RESPONSE_LATE, markers_list$HALLMARK_ESTROGEN_RESPONSE_EARLY)

  DEsample_KO_vs_WT <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_sample_KO_vs_WT_wilcox_de.xlsx')))

  DEsample_KO_vs_WT[gene %in% markers_list$HALLMARK_ESTROGEN_RESPONSE_EARLY & p_val_adj < 0.05 &
                      cluster1=='KO' & avg_logFC > 0] %>% nrow

  DEsample_KO_vs_WT[gene %in% markers_list$HALLMARK_ESTROGEN_RESPONSE_EARLY & p_val_adj < 0.05 &
                      cluster1=='WT' & avg_logFC > 0] %>% nrow

  DEsample_KO_vs_WT[gene %in% markers_list$HALLMARK_ESTROGEN_RESPONSE_LATE & p_val_adj < 0.05]

  #### Heatmap ER Response Genes ####
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         cells = names(which(arid1a.multiome$mammary_gland_cluster_classification=="L-Hor")),
                         features = DEsample_KO_vs_WT[gene %in% markers_list$HALLMARK_ESTROGEN_RESPONSE_EARLY & p_val_adj < 0.05 &
                                                        avg_logFC > 0.1, unique(gene)],
                         group.by = 'seurat_clusters',
                         #group.by = 'condition',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.colors = parula(25),
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         group.fontsize = 12,
                         show.group.dendrogram = T,
                         cluster_group_slices = T,
                         clustering.within = T,
                         zscore_data = F,
                         drop.min.value = NA,
                         label.value.min = 0,
                         #show.cell.dendrogram = T,
                         group.colors = seurat_cluster_palette
                         #group.colors = condition_palette
                         )
  f <- file.path(plot.dir,paste0('HALLMARK_ESTROGEN_RESPONSE_EARLY_L_Hor.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 9*res, height = 3*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("HALLMARK_ESTROGEN_RESPONSE_EARLY_UP"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         cells = names(which(arid1a.multiome$mammary_gland_cluster_classification=="L-Hor")),
                         features = DEsample_KO_vs_WT[gene %in% markers_list$HALLMARK_ESTROGEN_RESPONSE_LATE & p_val_adj < 0.05 &
                                                        avg_logFC > 0.1,  unique(gene)],
                         group.by = 'seurat_clusters',
                         #group.by = 'condition',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.colors = parula(25),
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         group.fontsize = 12,
                         show.group.dendrogram = T,
                         cluster_group_slices = T,
                         clustering.within = T,
                         zscore_data = F,
                         drop.min.value = NA,
                         label.value.min = 0,
                         #show.cell.dendrogram = T,
                         group.colors = seurat_cluster_palette
                         #group.colors = condition_palette
                         )
  f <- file.path(plot.dir,paste0('HALLMARK_ESTROGEN_RESPONSE_LATE_L_Hor.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 9*res, height = 3*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("HALLMARK_ESTROGEN_RESPONSE_LATE_UP"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  #### Heatmap DE_KO_Lhor_vs_WT_Lhor ####
  topn <- 50
  ss_arid <- subset(arid1a.multiome, seurat_clusters %in% c(2,3,4,5) )
  x <- DE_KO_Lhor_vs_WT_Lhor[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = ss_arid,
                         features = x,
                         group.by = 'condition',
                         assay = 'RNA',
                         slot = 'data',
                         heatmap.orientation = 'horizontal',
                         feature.fontsize = 8,
                         show.group.dendrogram = T,
                         #show.cell.dendrogram = T,
                         group.colors = condition_palette)
  f <- file.path(plot.dir,paste0('DE_KO_Lhor_vs_WT_Lhor_heatmap_features_top',topn,'.png'))
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 12*res, height = 2.1*res, res = res, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"KO vs WT in L-Hor (clusters 2,3 vs 4,5)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()
  image_write(image = fig, path = f, format = 'png')

  pdf(file = file.path(plot.dir,paste0('DE_KO_Lhor_vs_WT_Lhor_heatmap_features_top',topn,'.pdf')), width = 12, height = 2.1, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"KO vs WT in L-Hor (clusters 2,3 vs 4,5)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1))
  dev.off()

  #### Signature DE_KO_Lhor_vs_WT_Lhor ####
  DE_KO_Lhor_vs_WT_Lhor <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'scRNA_KO_Lhor_vs_WT_Lhor_wilcox_de.xlsx')))
  KO_Genes <- DE_KO_Lhor_vs_WT_Lhor[cluster1=='KO' & p_val_adj < 0.05 & avg_logFC >= 0.5, head(gene,34)]
  WT_Genes <- DE_KO_Lhor_vs_WT_Lhor[cluster1=='WT' & p_val_adj < 0.05 & avg_logFC >= 0.5, head(gene,34)]

  library(UCell)
  arid1a.multiome <- AddModuleScore_UCell(obj = arid1a.multiome, features = list("KO_Genes"=KO_Genes), assay = "RNA")
  arid1a.multiome <- AddModuleScore_UCell(obj = arid1a.multiome, features = list("WT_Genes"=WT_Genes), assay = "RNA")

  dtplot1 <- data.table('cluster'=arid1a.multiome@meta.data$seurat_clusters,
                        'value'=arid1a.multiome@meta.data$KO_Genes_UCell,
                        'label'='KO up')
  dtplot2 <- data.table('cluster'=arid1a.multiome@meta.data$seurat_clusters,
                        'value'=arid1a.multiome@meta.data$WT_Genes_UCell,
                        'label'='WT up')
  dtplot <- rbindlist(list(dtplot1, dtplot2))

  pos <- position_dodge(width=.95)
  g1 <- ggplot(dtplot,
               aes(x=cluster, y=value, fill=label)) +
    geom_violin(trim = F, scale = 'width', position = pos) +
    geom_boxplot(alpha=0.5, width=.5, position = pos, color='black') +
    scale_y_continuous(name = "Signature scores") +
    ggtitle("Organoid L-Hor gene signature") +
    xlab("Organoid clusters") +
    scale_fill_manual("Gene signature",values = c("red","gray")) +
    theme_bw(16)
  pdf(file = file.path(plot.dir,"gene_sig_L_hor.pdf"), width=14, height=6)
  g1
  dev.off()

  arid1a.multiome <- AddModuleScore(object = arid1a.multiome,
                                    features = list(KO_Genes),
                                    assay = 'RNA',
                                    name = 'KO_Lhor_vs_WT_Up', # appends a 1
                                    search = F) # searches HUGO contains human gene names only
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='KO_Lhor_vs_WT_Up1')] <- 'KO_Lhor_vs_WT_Up'

  arid1a.multiome <- AddModuleScore(object = arid1a.multiome,
                                    features = list(WT_Genes),
                                    assay = 'RNA',
                                    name = 'KO_Lhor_vs_WT_Down', # appends a 1
                                    search = F) # searches HUGO contains human gene names only
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='KO_Lhor_vs_WT_Down1')] <- 'KO_Lhor_vs_WT_Down'

  #### KO None vs WT None ####
  topn <- 50
  x <- DEsample_KO_None_vs_WT_None[!gene %in% rem][order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1]
  x <- unique(x$V1)
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x,
                         group.by = 'bioNames',
                         assay = 'RNA',
                         slot = 'data',
                         show.dendrogram = T,
                         disp.heatmap.vertically = F,
                         feature.font.size = 8,
                         group.colors = bionames_palette)
  png(filename = file.path(plot.dir,paste0('Heatmap_KOnone_vs_WTnone_DE_topn_',topn,'.png')),
      width = 12, height = 3, units = 'in', res = 200, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"(No E2) KO vs E2\n(FDR ≤ 0.05, | Log2FC ≥0.05)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1),
       heatmap_legend_side = "right")
  dev.off()

  #### KO E2 vs WT E2 ####
  topn <- 50
  x <- DEsample_KO_E2_vs_WT_E2[order(rank, decreasing = T)][p_val_adj < 0.05 & abs(avg_logFC) >= 0.5, head(gene, topn), by=cluster1] %>% unique()
  ht <- DoHeatmapEnhance(object = arid1a.multiome,
                         features = x$V1,
                         group.by = 'bioNames',
                         assay = 'RNA',
                         slot = 'data',
                         show.dendrogram = T,
                         disp.heatmap.vertically = F,
                         feature.font.size = 8,
                         group.colors = bionames_palette)
  png(filename = file.path(plot.dir,paste0('Heatmap_E2KO_vs_E2WT_DE_topn_',topn,'.png')),
      width = 12, height = 3, units = 'in', res = 200, bg = 'transparent')
  draw(ht,
       column_title = paste("TOP",topn,"E2 KO vs E2 WT\n(FDR ≤ 0.05, | Log2FC ≥0.05)"),
       column_title_gp = gpar(fontsize = 10, fontface = "bold", y=1.1),
       heatmap_legend_side = "right")
  dev.off()

  #### Volcano plot ####
  # samples
  library(dplyr)
  library(ggrepel)
  gvs <- lapply(unique(DEsample_v_others$cluster1), function(x){
    dt <- DEsample_v_others[cluster1==x]
    topn <- 20
    top_markers <- rbindlist(
      list(arrange(dt, rev(statistic)) %>% #group_by(sign(avg_logFC)) %>%
             filter(rank( rev(statistic), ties.method = "first") <= topn) %>%
             ungroup() %>% setDT,
           arrange(dt, statistic) %>% #group_by(sign(avg_logFC)) %>%
             filter(rank( statistic, ties.method = "first") <= topn) %>%
             ungroup() %>% setDT))

    xlim <- c(-max(abs(dt$avg_logFC)), max(abs(dt$avg_logFC)))
    ylim <- c(0, ceiling(max(dt$auc)))
    p1 <- ggplot(dt, aes(avg_logFC, pmax(-0.5, auc)) ) +
      geom_point(aes(color = p_val_adj < 0.05), pch=20) +
      scale_color_manual("Significant", values = c("gray60","firebrick1")) +
      scale_y_continuous(limits = ylim) +
      scale_x_continuous(limits = xlim) +
      geom_point(data = top_markers, color = "navy", shape=1) +
      ggtitle(paste('Sample',x,'vs others')) +
      geom_text_repel(data = top_markers[auc >= .5], mapping = aes(label = gene),
                      size          =2.8,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = 0.1,
                      direction    = "x",
                      angle        = 90,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      max.time = 1) +
      geom_text_repel(data = top_markers[auc < .5], mapping = aes(label = gene),
                      size          =2.8,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = -0.1,
                      direction    = "x",
                      angle        = 90,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      max.time = 1) +
      theme_bw(9) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
      ylab("AUC") + xlab("mean log2 fold-change")
  })
  gvs[[1]]
  png(filename = file.path(plot.dir,'Sample_volcano.png'), width = 18, height = 9, units = 'in', res=150, bg='transparent')
  do.call(grid.arrange, c(gvs, nrow=2))
  dev.off()

  # clusters
  gvs <- lapply(sort(as.integer(unique(DEcluster_v_others$cluster1))), function(x){
    dt <- DEcluster_v_others[cluster1==x]
    topn <- 20
    top_markers <- rbindlist(
      list(arrange(dt, rev(statistic)) %>% #group_by(sign(avg_logFC)) %>%
             filter(rank( rev(statistic), ties.method = "first") <= topn) %>%
             ungroup() %>% setDT,
           arrange(dt, statistic) %>% #group_by(sign(avg_logFC)) %>%
             filter(rank( statistic, ties.method = "first") <= topn) %>%
             ungroup() %>% setDT))

    xlim <- c(-max(abs(dt$avg_logFC)), max(abs(dt$avg_logFC)))
    ylim <- c(0, ceiling(max(dt$auc)))
    p1 <- ggplot(dt, aes(avg_logFC, pmax(-0.5, auc)) ) +
      geom_point(aes(color = p_val_adj < 0.05), pch=20) +
      scale_color_manual("Significant", values = c("gray60","firebrick1")) +
      scale_y_continuous(limits = ylim) +
      scale_x_continuous(limits = xlim) +
      geom_point(data = top_markers, color = "navy", shape=1) +
      ggtitle(paste('Cluster',x,'vs others')) +
      geom_text_repel(data = top_markers[auc >= .5], mapping = aes(label = gene),
                      size          =2.8,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = 0.1,
                      direction    = "x",
                      angle        = 90,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      max.time = 1) +
      geom_text_repel(data = top_markers[auc < .5], mapping = aes(label = gene),
                      size          =2.8,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = -0.1,
                      direction    = "x",
                      angle        = 90,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      max.time = 1) +
      theme_bw(9) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
      ylab("AUC") + xlab("mean log2 fold-change")
  })
  gvs[[4]]
  png(filename = file.path(plot.dir,'Cluster_volcano.png'), width = 26, height = 16, units = 'in', res=150, bg='white')
  do.call(grid.arrange, c(gvs, nrow=4))
  dev.off()

  #### Correlation in sample gene expression ####
  library(corrplot)

  par(bg=NA)
  # average expression per sample
  Idents(arid1a.multiome.sct) <- 'bioNames'
  avexp <- AverageExpression(arid1a.multiome.sct, assays = 'SCT')$SCT
  png(filename = file.path(plot.dir,'sample_expression_correlation_matrix.png'), height = 8, width = 8, units = 'in', bg = 'transparent', res = 150)
  #plot(log1p(avexp), cex=.6)
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan",...)
  }
  pairs(log1p(avexp), #main='Avg. expression',
        panel = panel.smooth, upper.panel = NULL,
        cex = .6, pch = 1, bg = "transparent", horOdd=TRUE,
        diag.panel = panel.hist, cex.labels = 1, font.labels = 2)
  dev.off()

  old.par <- par('mar') #tlbr
  par(mar=c(12,4.1,4.1,2.1))
  cairo_pdf(filename = file.path(plot.dir,'sample_expression_correlation_matrix.pdf'), height = 7, width = 7,
            bg = 'transparent',
            onefile = T)
  corrplot(cor(avexp), type="upper", is.corr = F,
           order="hclust",
           method='color',
           tl.col = 'black',
           addCoef.col = 'black',
           bg='transparent',
           outline = T,
           #main = "Avg. expression",
           col=rev(pals::brewer.rdbu(8)),
           hclust.method = 'ward.D2')

  #### Correlation in sample PCA space ####
  pcaemb <- data.table(data.frame(arid1a.multiome.sct@reductions$pca@cell.embeddings), keep.rownames = 'sample')
  pcaemb$sample <- arid1a.multiome.sct@meta.data[pcaemb$sample, 'bioNames']
  d <- pcaemb[, lapply(.SD, mean), by=sample]
  mat <- as.matrix(d[,-1], rownames=d$sample)
  corrplot(cor(t(mat)), type="upper",
           order="hclust",
           method='color',
           tl.col = 'black',
           addCoef.col = 'black',
           outline = T,
           #tl.pos = "ld",
           bg='transparent',
           #main = "Avg. PCA embeddings",
           col=rev(pals::brewer.rdbu(8)),
           hclust.method = 'ward.D2')
  dev.off()
  par(mar=old.par)

  #### Sorted bulk populations correlate w/ single cells ####
  # First check which bulk samples to use
  dtcnt <- fread(file.path(data.dir, 'bulk_rna_seq_counts_normalized.tsv'))
  # wouter says the CD49 and CD24 are from mouse biological replicates
  dtcnt[,gene := Hmisc::capitalize(tolower(gene))]
  # First check which bulk samples to use
  # basal
  dtbasalnorm <- dtcnt[,grep('cd49', colnames(dtcnt), ignore.case = T, value = T), with=F]
  bct <- cor.test(dtbasalnorm$CD49f_pos_R_3, dtbasalnorm$MP_2CD49f_pos)

  # luminal
  dtluminalnorm <- dtcnt[,grep('cd24', colnames(dtcnt), ignore.case = T, value = T), with=F]
  lct <- cor.test(dtluminalnorm$MP_2CD24_pos, dtluminalnorm$MP_3CD24_pos)

  # plot - issue is poor correlation in basal.  Samir suggests `MP-2CD49f-pos`
  par(mfrow=c(2,1))
  plot(dtbasalnorm+1, log='xy', panel.first=grid()); abline(a=0, b=1, col='red');
  mtext(text = paste('Basal Bulk r=',round(bct$estimate, digits = 3), 'pval=', bct$p.value))
  plot(dtluminalnorm+1, log='xy', panel.first=grid()); abline(a=0, b=1, col='red');
  mtext(text = paste('Luminal Bulk r=',round(lct$estimate, digits = 3), 'pval=', lct$p.value))

  panel.cor <- function(x, y, ...)
  {
    par(usr = c(0, 1, 0, 1))
    txt <- as.character(format(cor(x, y, method = 'pearson'), digits=2))
    text(0.5, 0.5, txt, cex = abs(4 * cor(x, y)))
  }
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan",...)
  }
  lower.panel <- function(x, y,...){
    col = c("gray20", "red", "cornflowerblue")[group]
    shape = c(1, 19, 19)[group]
    points(x,y, pch=shape, col=col)
    #abline(a = 0,b = 1, col='black')
    fit <- loess(y ~ x);
    j <- order(x)
    lines(x[j],fit$fitted[j],col="red",lwd=1)
  }

  ## correlate with single cell by sample
  Idents(arid1a.multiome.sct) <- 'bioNames'
  bionames <- unique(arid1a.multiome.sct$bioNames)
  ##dtAvgExpSample <- AverageExpression(arid1a.multiome.sct, assays = 'SCT', slot = 'data')$SCT %>% data.table(., keep.rownames = 'gene')
  dtTotalExpSample <- TotalExpression(arid1a.multiome.sct, assays = 'SCT', slot = 'data')$SCT %>% data.table(., keep.rownames = 'gene')
  basal <- c('MP_2CD49f_pos','CD49f_pos_R_3')
  luminal <- c('MP_2CD24_pos','MP_3CD24_pos')
  DT <- dtTotalExpSample
  m <- merge(DT, dtcnt, by='gene')

  dtplot <- m[,colnames(DT), with=F]
  dtplot$LuminalBulk <- rowMeans(m[,luminal, with=F])
  dtplot$BasalBulk <- rowMeans(m[,basal, with=F])
  # Every gene must have minimum counts in any experiment
  min.cnt <- 40
  dtplot <- dtplot[dtplot[, all(.SD>= min.cnt), .SDcols=bionames, by=1:nrow(dtplot)]$V1]
  dtplot <- dtplot[dtplot[, all(.SD>= min.cnt), .SDcols=c('LuminalBulk','BasalBulk'), by=1:nrow(dtplot)]$V1]
  dtplot$group <- 1
  dtplot[gene %in% markers_list$`Prostate luminal mature`, group:=2]
  dtplot[gene %in% markers_list$`Prostate basal`, group:=3]
  dtplot <- dtplot[order(group)] #re-order to luminal and basal colors appear on top
  group <- dtplot$group
  png(filename = file.path(plot.dir,'sc_bulk_sample_expression_correlation_matrix.png'), height = 12, width = 12, units = 'in', res = 150, bg='transparent')
  #pairs(log1p(dtplot[,-c("gene","group"), with=F]), #main='Avg. expression',
  pairs(log1p(dtplot[,-c("gene","group"), with=F]), #main='Avg. expression',
        panel = panel.smooth, upper.panel = panel.cor, lower.panel = lower.panel,
        cex = .6, pch = 1, bg = "transparent", horOdd=TRUE, cex.labels = 1, font.labels = 2,
        diag.panel = panel.hist)
  dev.off()

  dtplot <- m[,colnames(DT), with=F]
  dtplot$LuminalBulk <- rowMeans(m[,luminal, with=F])
  dtplot$BasalBulk <- rowMeans(m[,basal, with=F])

  png(filename = file.path(plot.dir,'BulkCorrelated.png'), width = 18, height = 6, units = 'in', res = 150, bg = 'transparent')
  par(mfrow=c(2,length(bionames)))
  lapply(bionames, function(i){
    dt <- dtplot[,list(get(i), LuminalBulk)]
    dt <- log1p(dt)
    dt <- dt[V1 >= 1]
    dt <- dt[LuminalBulk >= 1]
    lim <- round(max(dt))
    lct <- cor.test(dt[[1]], dt[[2]], method = 'pearson')
    setnames(dt, 'V1', i)
    plot( dt , main=paste(i,'vs Luminal Bulk'),  panel.first=grid(), col='dodgerblue', pch=20, xlim=c(0,lim), ylim=c(0,lim));
    mtext(text = paste('Pearson r=',round(lct$estimate, digits = 3), 'pval=', lct$p.value))
  })
  lapply(bionames, function(i){
    dt <- dtplot[,list(get(i), BasalBulk)]
    dt <- log1p(dt)
    dt <- dt[V1 >= 1]
    dt <- dt[BasalBulk >= 1]
    lim <- round(max(dt))
    lct <- cor.test(dt[[1]], dt[[2]], method = 'pearson')
    setnames(dt, 'V1', i)
    plot( dt , main=paste(i,'vs Basal Bulk'),  panel.first=grid(), col='dodgerblue', pch=20, xlim=c(0,lim), ylim=c(0,lim));
    mtext(text = paste('Pearson r=',round(lct$estimate, digits = 3), 'pval=', lct$p.value))
  })
  dev.off()

  cairo_pdf(filename = file.path(plot.dir,'sample_expression_correlation_matrix_bulk.pdf'), height = 6.5, width = 6.5,
            bg = 'transparent',
            onefile = T)
  d <- as.matrix(log1p(dtplot[,-c("gene"), with=F]))
  d[which(d<1)] <- NA
  corrplot(cor(na.omit(d), method = 'pearson'),
           type="upper",
           order="hclust",
           method='color',
           tl.col = 'black',
           addCoef.col = 'white',
           addCoefasPercent = T,
           diag = T,
           outline = T,
           #tl.pos = "ld",
           bg='transparent',
           #main = "Avg. PCA embeddings",
           col=rev(pals::brewer.rdbu(100)),
           hclust.method = 'ward.D2')
  dev.off()

  ## correlate by cluster with bulk
  Idents(arid1a.multiome.sct) <- 'seurat_clusters'
  dtTotalExpCluster <- TotalExpression(arid1a.multiome.sct, assays = 'SCT', slot = 'data')$SCT
  clusters <- colnames(dtTotalExpCluster)
  dtTotalExpCluster <- data.table(dtTotalExpCluster, keep.rownames = 'gene')
  m <- merge(dtTotalExpCluster, dtcnt, by='gene')
  ssexp <- grep("gene", colnames(dtTotalExpCluster), invert = T, value = T)
  basal <- c('MP_2CD49f_pos','CD49f_pos_R_3')
  luminal <- c('MP_2CD24_pos','MP_3CD24_pos')
  min.cnt <- 40
  dtplot <- m[,colnames(dtTotalExpCluster), with=F]
  dtplot$LuminalBulk <- rowMeans(m[,luminal, with=F])
  dtplot$BasalBulk <- rowMeans(m[,basal, with=F])
  # Every gene must have 20 more counts in any experiment
  dtplot <- dtplot[dtplot[, any(.SD>= min.cnt), .SDcols=clusters, by=1:nrow(dtplot)]$V1]
  dtplot <- dtplot[dtplot[, all(.SD>= min.cnt), .SDcols=c('LuminalBulk','BasalBulk'), by=1:nrow(dtplot)]$V1]

  png(filename = file.path(plot.dir,'sc_bulk_cluster_expression_correlation_matrix_sep.png'), height = 24, width = 24, units = 'in', res = 150, bg='transparent')
  pairs(log1p(dtplot[,-c("gene"), with=F]), #main='Avg. expression',
        panel = panel.smooth, upper.panel = panel.cor, #lower.panel = lower.panel,
        cex = .6, pch = 1, bg = "transparent", horOdd=TRUE, cex.labels = 1, font.labels = 2,
        diag.panel = panel.hist)
  dev.off()

  cairo_pdf(filename = file.path(plot.dir,'cluster_expression_correlation_matrix_bulk.pdf'), height = 9.5, width = 9.5,
            bg = 'transparent',
            onefile = T)
  d <- log1p(dtplot[,-c("gene"), with=F])
  corrplot(cor(d), type="upper",
           order="hclust",
           method='color',
           tl.col = 'black',
           addCoef.col = 'white',
           addCoefasPercent = T,
           diag = T,
           outline = T,
           #tl.pos = "ld",
           bg='transparent',
           #main = "Avg. PCA embeddings",
           col=rev(pals::brewer.rdbu(100)),
           hclust.method = 'ward.D2')
  dev.off()

  ### transfer counts to DESEQ ###
  Idents(arid1a.multiome.sct) <- 'bioNames'
  dtSampleCounts <- TotalExpression(arid1a.multiome.sct, assays = 'RNA', slot = 'counts')$RNA
  dtSampleCounts <- setDT(dtSampleCounts, keep.rownames = 'gene')
  gz1 <- gzfile(file.path(data.dir, 'singlecell_rna_seq_counts_sample.csv.gz'), "w")
  write.csv(dtSampleCounts, gz1, row.names = F)
  close(gz1)

  Idents(arid1a.multiome.sct) <- 'seurat_clusters'
  dtClusterCounts <- TotalExpression(arid1a.multiome.sct, assays = 'RNA', slot = 'counts')$RNA
  dtClusterCounts <- setDT(dtClusterCounts, keep.rownames = 'gene')
  gz1 <- gzfile(file.path(data.dir, 'singlecell_rna_seq_counts_cluster.csv.gz'), "w")
  write.csv(dtClusterCounts, gz1, row.names = F)
  close(gz1)

  #### Mammary tissue signature L-HOR ####
  dt <- fread(file.path(data.dir, 'DEcluster_v_others_lhor.txt'))
  library(UCell)
  arid1a.multiome <- AddModuleScore_UCell(obj = arid1a.multiome, assay = "RNA",
                                          features = list("C2_mammary"=dt[cluster==2 & p_val_adj < 0.01 & myAUC > 0.65, head(gene, 15)]))
  arid1a.multiome <- AddModuleScore_UCell(obj = arid1a.multiome, assay = "RNA",
                                          features = list("C3_mammary"=dt[cluster==3 & p_val_adj < 0.01 & myAUC > 0.65, head(gene, 15)]))
  arid1a.multiome <- AddModuleScore_UCell(obj = arid1a.multiome, assay = "RNA",
                                          features = list("C4_mammary"=dt[cluster==4 & p_val_adj < 0.01 & myAUC > 0.65, head(gene, 15)]))
  dt <- data.table(arid1a.multiome@meta.data, keep.rownames = "cells")

  dt$mammary_lhor_UCell <- ''
  dt[,mammary_lhor_UCell := names(.SD)[max.col(.SD, ties.method = "first")], .SDcols=c('C2_mammary_UCell', 'C3_mammary_UCell', 'C4_mammary_UCell')]
  if(identical(Cells(arid1a.multiome), dt$cells)){
    arid1a.multiome$mammary_lhor_UCell <- dt$mammary_lhor_UCell
  } else {
    stop("Cells are not identical")
  }

  table(arid1a.multiome$mammary_lhor_UCell, arid1a.multiome$seurat_clusters)

  g <- DimPlot(arid1a.multiome, group.by = "mammary_lhor_UCell") %>%
    UmapPlotEnhance(obj = arid1a.multiome,
                    p1 = .,
                    axis.label = "mammary_lhor_UCell",
                    pt.size = 1,
                    order.random = F,
                    title = "mammary_lhor_UCell",
                    legend.title = "score",
                    pal = c("red",'blue','black'))
  g

  #### Save the arid1a object ####
  saveRDS(arid1a.multiome, arid1a.object)

  ### Score AR target markers ####
  arid1a.multiome <- AddModuleScore(object = arid1a.multiome,
                                    features = list(nelson_mouse_andr_up),
                                    assay = 'RNA',
                                    name = 'nelson_mouse_andr_up', # appends a 1
                                    search = F) # searches HUGO contains human gene names only
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='nelson_mouse_andr_up1')] <- 'nelson_mouse_andr_up'
  ggplot(data.table('cluster'=arid1a.multiome@meta.data$seurat_clusters,
                    'nelson_mouse_andr_up'=arid1a.multiome@meta.data$nelson_mouse_andr_up),
         aes(x=cluster, y=nelson_mouse_andr_up)) + geom_violin(trim = F, draw_quantiles = T, scale = 'width')

  arid1a.multiome <- AddModuleScore(object = arid1a.multiome,
                                   features = list(markers_list$HALLMARK_ANDROGEN_RESPONSE),
                                   assay = 'RNA',
                                   name = 'ar_targets', # appends a 1
                                   search = F) # searches HUGO contains human gene names only
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='ar_targets1')] <- 'ar_targets'

  ggplot(data.table('cluster'=arid1a.multiome@meta.data$seurat_clusters,
                    'Ar_targets'=arid1a.multiome@meta.data$ar_targets),
         aes(x=cluster, y=Ar_targets)) + geom_violin(trim = F, draw_quantiles = T, scale = 'width')

  colors <- colorRampPalette(c("royalblue4","royalblue2","deepskyblue",
                                           "grey90","orangered3","red3","darkred"), space="rgb")(1024);

  #### UMAP Ar ####
  glabels <- DimPlot(arid1a.multiome, group.by = 'seurat_clusters',
                     pt.size = 0.2, ncol = 1)
  glabels$data$seurat_clusters <- as.integer(glabels$data$seurat_clusters)
  glabels <- UmapPlotEnhance(obj = arid1a.multiome,
                             p1 = glabels,
                             addLabels = T,
                             label.size = 4,
                             pal=NA,
                             pt.size = 1)

  f <- 'Ar'
  g <- FeaturePlot(arid1a.multiome, features = f) %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    axis.label = "",
                    pt.size = 1,
                    scale.color=F,
                    title = "Ar",
                    legend.title = "Expression",
                    pal = rocket_palette, # pals::kovesi.linear_kryw_5_100_c67(10),
                    filename.prefix = 'Ar_Umap')
  g

  #### UMAP Ar targets ####
  g2 <- FeaturePlot(arid1a.multiome,
                    features = c('nelson_mouse_andr_up'),
                    slot = 'data',
                    min.cutoff = 'q1',
                    max.cutoff = 'q99'
  ) %>%
    UmapPlotEnhance(obj = arid1a.multiome,
                    p1 = .,
                    panel.border.size = 0,
                    break.unit = 0.1,
                    axis.label = NA,
                    scale.color = F,
                    label.size = 5,
                    pt.size = .8,
                    title.size = 15,
                    title = "Androgen up expression scores",
                    legend.title = "",
                    pal = rocket_palette,
                    filename.prefix = "Androgen_up_umap"
                    #pal= pals::ocean.thermal(20) # color_blind_parula_palette,
    )

  g2 <- FeaturePlot(arid1a.multiome,
                    features = c('ar_targets'),
                    slot = 'data',
                    min.cutoff = 'q1',
                    max.cutoff = 'q99'
  ) %>%
    UmapPlotEnhance(obj = arid1a.multiome,
                    p1 = .,
                    panel.border.size = 0,
                    break.unit = 0.1,
                    axis.label = NA,
                    scale.color = F,
                    label.size = 5,
                    pt.size = .5,
                    title.size = 15,
                    title = "Ar target expression scores",
                    legend.title = "",
                    pal = rocket_palette,
                    #pal= pals::ocean.thermal(20) # color_blind_parula_palette,
    )

  dt <- data.table(arid1a.multiome@meta.data[,c('bioNames',
                                                'condition',
                                                'seurat_clusters',
                                                'mammary_gland_cluster_classification',
                                                'ar_targets')])
  ### violin Ar target scores ####
  library(ggbeeswarm)
  #df$ident <- as.character(df$ident)
  gb <- ggplot(dt, aes(x=seurat_clusters, y=ar_targets, fill=seurat_clusters, color=mammary_gland_cluster_classification)) +
    geom_quasirandom(aes(x=seurat_clusters), cex = .05, varwidth = 1, width = .5, dodge.width = 2, pch = 20) +
    geom_boxplot(width=0.1, color="black", alpha=0.2) +
    #scale_x_discrete(breaks = c("Arid1a(-/-)", "Arid1a(+/+)", "Arid1a(-/-) E2", "Arid1a(+/+) E2"), expand = c(0,0)) +
    scale_fill_manual("", values = mammary_cluster_class_palette) +
    scale_color_manual("", values = mammary_cluster_class_palette) +
    xlab("Cluster") + ylab("Scores") + ggtitle("Ar Targets scores") +
    #scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
    theme_classic(14) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          legend.position = 'right',
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=9)))
  ### CDF Ar target scores ####
  gc <- ggplot(dt, aes(x=ar_targets, color=seurat_clusters)) +
    stat_ecdf(geom = "step", linewidth=1.5) +
    ggtitle("Ar Targets scores CDF") +
    ylab("Data (%)") + xlab("Ar target scores") +
    #scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(limits = c(-0.1, 0.3)) +
    geom_hline(yintercept = 0.5, linetype='dotted') +
    geom_vline(xintercept = 0.5, linetype='dotted') +
    scale_color_manual("Clusters", values = godsnot_102) +
    theme_bw(16) + theme(text = element_text(family='Arial'),
                         axis.text = element_text(size=12),
                         legend.text = element_text(size=12),
                         legend.key.height = unit(.5, 'cm'),
                         legend.key.size = unit(.5, 'cm'),
                         legend.spacing.y = unit(.02, 'cm'),
                         panel.grid = element_blank(),
                         plot.background = element_blank(),
                         panel.background = element_blank(),
                         legend.background = element_blank(),
                         strip.background = element_blank(),
                         strip.text = element_text(size=12, face='bold'))
  png(filename = file.path(plot.dir,'AR_target_scores.png'), width = 14, height = 3.5, units = 'in', res = 250, bg = 'transparent')
  grid.arrange(g + glabels$layers, g2 + glabels$layers, gb, widths=c(.3,.3,.4))
  dev.off()

  ### Score Hallmark ER early ####
  arid1a.multiome <- AddModuleScore(object = arid1a.multiome,
                                    features = list(markers_list$HALLMARK_ESTROGEN_RESPONSE_EARLY),
                                    assay = 'RNA',
                                    name = 'er_response_early', # appends a 1
                                    search = F) # searches HUGO contains human gene names only
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='er_response_early1')] <- 'er_response_early'

  ggplot(data.table('cluster'=arid1a.multiome@meta.data$seurat_clusters,
                    'er_response_early'=arid1a.multiome@meta.data$er_response_early),
         aes(x=cluster, y=er_response_early)) + geom_violin(trim = F, draw_quantiles = T, scale = 'width')

  ### Score Hallmark ER late ####
  arid1a.multiome <- AddModuleScore(object = arid1a.multiome,
                                    features = list(markers_list$HALLMARK_ESTROGEN_RESPONSE_LATE),
                                    assay = 'RNA',
                                    name = 'er_response_late', # appends a 1
                                    search = F) # searches HUGO contains human gene names only
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='er_response_late1')] <- 'er_response_late'

  ggplot(data.table('cluster'=arid1a.multiome@meta.data$seurat_clusters,
                    'er_response_late'=arid1a.multiome@meta.data$er_response_late),
         aes(x=cluster, y=er_response_late)) + geom_violin(trim = F, draw_quantiles = T, scale = 'width')

  #### UMAP Er ####
  glabels <- DimPlot(arid1a.multiome, group.by = 'seurat_clusters',
                     pt.size = 0.2, ncol = 1)
  glabels$data$seurat_clusters <- as.integer(glabels$data$seurat_clusters)
  glabels <- UmapPlotEnhance(obj = arid1a.multiome,
                             p1 = glabels,
                             addLabels = T,
                             label.size = 4,
                             pal=NA,
                             pt.size = 1)

  f <- 'Esr1'
  g <- FeaturePlot(arid1a.multiome, features = f) %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    axis.label = "",
                    pt.size = 1,
                    scale.color=F,
                    title = f,
                    legend.title = "Expression",
                    pal = rocket_palette, # pals::kovesi.linear_kryw_5_100_c67(10),
                    filename.prefix = paste0(f,'_Umap'))
  g
  #### UMAP ER targets ####
  g2 <- FeaturePlot(arid1a.multiome,
                    features = c('er_response_early'),
                    slot = 'data',
                    min.cutoff = 'q1',
                    max.cutoff = 'q99'
  ) %>%
    UmapPlotEnhance(obj = arid1a.multiome,
                    p1 = .,
                    panel.border.size = 0,
                    break.unit = 0.1,
                    axis.label = NA,
                    scale.color = F,
                    label.size = 5,
                    pt.size = .8,
                    title.size = 15,
                    title = "Er response early scores",
                    legend.title = "",
                    pal = rocket_palette,
                    filename.prefix = "UMAP_ER_early",
                    as.pdf = T
    )

  g3 <- FeaturePlot(arid1a.multiome,
                    features = c('er_response_late'),
                    slot = 'data',
                    min.cutoff = 'q1',
                    max.cutoff = 'q99'
  ) %>%
    UmapPlotEnhance(obj = arid1a.multiome,
                    p1 = .,
                    panel.border.size = 0,
                    break.unit = 0.1,
                    axis.label = NA,
                    scale.color = F,
                    label.size = 5,
                    pt.size = .8,
                    title.size = 15,
                    title = "Er response late scores",
                    legend.title = "",
                    pal = rocket_palette,
                    filename.prefix = "UMAP_ER_late",
                    as.pdf = T
                    #pal= pals::ocean.thermal(20) # color_blind_parula_palette,
    )

  require(ggrastr)
  library(ggbeeswarm)
  library(ggpubr)

  dt <- data.table(arid1a.multiome@meta.data[,c('bioNames',
                                                'condition',
                                                'seurat_clusters',
                                                'mammary_gland_cluster_classification',
                                                'er_response_early',
                                                "er_response_late")], keep.rownames = 'cells')
  ### violin esr1 ####
  Esr1 <- arid1a.multiome@assays$RNA@data['Esr1',]
  dt <- merge(dt, as.data.table(Esr1, keep.rownames = 'cells'), by='cells')
  gbesr1 <- ggplot(dt, aes(x=seurat_clusters, y=Esr1, fill=seurat_clusters, color=mammary_gland_cluster_classification)) +
    geom_quasirandom(aes(x=seurat_clusters), cex = .05, varwidth = 1, width = .5, dodge.width = 2, pch = 20) +
    geom_boxplot(width=0.1, color="black", alpha=0.2) +
    scale_fill_manual("", values = mammary_cluster_class_palette) +
    scale_color_manual("", values = mammary_cluster_class_palette) +
    xlab("Cluster") + ylab("Scores") + ggtitle("Esr1 expression") +
    theme_classic(14) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          legend.position = 'right',
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=9)))
  ### CDF er_response_early target scores ####
  gce <- ggplot(dt, aes(x=er_response_early, color=seurat_clusters)) +
    stat_ecdf(geom = "step", linewidth=1.5) +
    ggtitle("Er response early scores CDF") +
    ylab("Data (%)") + xlab("Ar target scores") +
    #scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(limits = c(-0.1, 0.3)) +
    geom_hline(yintercept = 0.5, linetype='dotted') +
    geom_vline(xintercept = 0.5, linetype='dotted') +
    scale_color_manual("Clusters", values = godsnot_102) +
    theme_bw(16) + theme(text = element_text(family='Arial'),
                         axis.text = element_text(size=12),
                         legend.text = element_text(size=12),
                         legend.key.height = unit(.5, 'cm'),
                         legend.key.size = unit(.5, 'cm'),
                         legend.spacing.y = unit(.02, 'cm'),
                         panel.grid = element_blank(),
                         plot.background = element_blank(),
                         panel.background = element_blank(),
                         legend.background = element_blank(),
                         strip.background = element_blank(),
                         strip.text = element_text(size=12, face='bold'))
  ### CDF er_response_early target scores ####
  gcl <- ggplot(dt, aes(x=er_response_late, color=seurat_clusters)) +
    stat_ecdf(geom = "step", linewidth=1.5) +
    ggtitle("Er response late scores CDF") +
    ylab("Data (%)") + xlab("Ar target scores") +
    #scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(limits = c(-0.1, 0.3)) +
    geom_hline(yintercept = 0.5, linetype='dotted') +
    geom_vline(xintercept = 0.5, linetype='dotted') +
    scale_color_manual("Clusters", values = godsnot_102) +
    theme_bw(16) + theme(text = element_text(family='Arial'),
                         axis.text = element_text(size=12),
                         legend.text = element_text(size=12),
                         legend.key.height = unit(.5, 'cm'),
                         legend.key.size = unit(.5, 'cm'),
                         legend.spacing.y = unit(.02, 'cm'),
                         panel.grid = element_blank(),
                         plot.background = element_blank(),
                         panel.background = element_blank(),
                         legend.background = element_blank(),
                         strip.background = element_blank(),
                         strip.text = element_text(size=12, face='bold'))

  png(filename = file.path(plot.dir,'er_response_target_scores.png'), width = 12, height = 6, units = 'in', res = 250, bg = 'transparent')
  grid.arrange(gbesr1, gce, gcl)
  dev.off()

  ### er_response_early target scores ####
  gbe <- ggplot(dt, aes(x=seurat_clusters, y=er_response_early, fill=seurat_clusters, color=mammary_gland_cluster_classification)) +
    geom_quasirandom(aes(x=seurat_clusters), cex = 1, varwidth = 1, width = .5, dodge.width = 2, pch = 20) +
    geom_boxplot(width=0.1, color="black", alpha=0.2) +
    #scale_x_discrete(breaks = c("Arid1a(-/-)", "Arid1a(+/+)", "Arid1a(-/-) E2", "Arid1a(+/+) E2"), expand = c(0,0)) +
    scale_fill_manual("", values = mammary_cluster_class_palette) +
    scale_color_manual("", values = mammary_cluster_class_palette) +
    xlab("Cluster") + ylab("Scores") + ggtitle("Er response early scores") +
    #scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
    theme_classic(14) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          legend.position = 'right',
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=9)))
  ### er_response_late target scores ####
  gbl <- ggplot(dt, aes(x=seurat_clusters, y=er_response_late, fill=seurat_clusters, color=mammary_gland_cluster_classification)) +
    geom_quasirandom(aes(x=seurat_clusters), cex = 1, varwidth = 1, width = .5, dodge.width = 2, pch = 20) +
    geom_boxplot(width=0.1, color="black", alpha=0.2) +
    #scale_x_discrete(breaks = c("Arid1a(-/-)", "Arid1a(+/+)", "Arid1a(-/-) E2", "Arid1a(+/+) E2"), expand = c(0,0)) +
    scale_fill_manual("", values = mammary_cluster_class_palette) +
    scale_color_manual("", values = mammary_cluster_class_palette) +
    xlab("Cluster") + ylab("Scores") + ggtitle("Er response late scores") +
    #scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
    theme_classic(14) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          legend.position = 'right',
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=9)))

  pdf(file = file.path(plot.dir,'er_response_target_scores.pdf'), width = 16, height = 5, bg = 'transparent')
  grid_arrange_list_shared_legend(plots = list(gbe,gbl), legend.position = "topright")
  dev.off()

  ### er_response_early target scores ####
  pos_dodge <- position_dodge(1)
  gbe <- ggplot(dt[mammary_gland_cluster_classification=='L-Hor'],
                aes(x=bioNames, y=er_response_early, fill=bioNames, color=mammary_gland_cluster_classification)) +
    #geom_quasirandom(aes(x=bioNames), cex = 1, varwidth = 1, width = .5, dodge.width = 3, pch = 20) +
    geom_violin(aes(x=bioNames), scale = "width",  cex = 1, width = .8, position = pos_dodge) +
    geom_boxplot(width=0.1, color="black", alpha=0.2, position = pos_dodge) +
    #scale_x_discrete(breaks = c("Arid1a(-/-)", "Arid1a(+/+)", "Arid1a(-/-) E2", "Arid1a(+/+) E2"), expand = c(0,0)) +
    scale_fill_manual("", values = bionames_palette) +
    scale_color_manual("", values = bionames_palette) +
    xlab("Cluster") + ylab("Scores") + ggtitle("Er response early scores in L-Hor cells") +
    stat_compare_means(method='anova', label.y=0.3) +
    stat_compare_means(label = "p.signif", method = "wilcox", ref.group = ".all.") +
    #scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
    theme_classic(14) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          legend.position = 'right',
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=9)))

  ### er_response_late target scores ####
  pos_dodge <- position_dodge(1)
  gbl <- ggplot(dt[mammary_gland_cluster_classification=='L-Hor'],
                aes(x=bioNames, y=er_response_late, fill=bioNames, color=mammary_gland_cluster_classification)) +
    #geom_quasirandom(aes(x=bioNames), cex = 1, varwidth = 1, width = .5, dodge.width = 3, pch = 20) +
    geom_violin(aes(x=bioNames), scale = "width",  cex = 1, width = .8, position = pos_dodge) +
    geom_boxplot(width=0.1, color="black", alpha=0.2, position = pos_dodge) +
    #scale_x_discrete(breaks = c("Arid1a(-/-)", "Arid1a(+/+)", "Arid1a(-/-) E2", "Arid1a(+/+) E2"), expand = c(0,0)) +
    scale_fill_manual("", values = bionames_palette) +
    scale_color_manual("", values = bionames_palette) +
    xlab("Cluster") + ylab("Scores") + ggtitle("Er response late scores in L-Hor cells") +
    stat_compare_means(method='anova', label.y=0.3) +
    stat_compare_means(label = "p.signif", method = "wilcox", ref.group = ".all.") +
    theme_classic(14) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          legend.position = 'right',
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=9)))

  pdf(file = file.path(plot.dir,'er_response_target_scores_L-Hor.pdf'), width = 12, height = 5, bg = 'transparent')
  grid_arrange_list_shared_legend(plots = list(gbe,gbl), legend.position = "topright")
  dev.off()

  #### DE Conserved markers ####
  # find if any genes are different in mutants vs empty vector
  FindConservedMarkers(seurat_obj,
                       ident.1 = cluster,
                       grouping.var = "group",
                       only.pos = TRUE)

  #### Define Luminal/Basal/Progenitors ####
  ## label clusters based on DE and GSEAs ##
  arid1a.multiome@meta.data$cluster_category <- ""
  arid1a.multiome@meta.data[which(arid1a.multiome@meta.data$seurat_clusters %in% c(1,2,5,6,7,8,13,15,16)),]$cluster_category <- 'Mammary Luminal mature'
  arid1a.multiome@meta.data[which(arid1a.multiome@meta.data$seurat_clusters %in% c(3,4,9,10,11,17,18)),]$cluster_category <- 'Mammary Luminal progenitor'
  arid1a.multiome@meta.data[which(arid1a.multiome@meta.data$seurat_clusters %in% c(9)),]$cluster_category <- 'Mammary Basal'
  arid1a.multiome@meta.data[which(arid1a.multiome@meta.data$seurat_clusters %in% c(12)),]$cluster_category <- 'cluster 12'
  table(arid1a.multiome@meta.data[,c("seurat_clusters","cluster_category")])

  #### DE Markers (luminal and basal cells) ####
  # https://github.com/immunogenomics/presto
  library(presto)
  # bionames
  DE_cell_category_v_others <- wilcoxauc(arid1a.multiome, group_by = 'cluster_category', seurat_assay = 'RNA', assay = 'data')
  DE_cell_category_v_others <- setDT(DE_cell_category_v_others)
  setnames(DE_cell_category_v_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DE_cell_category_v_others[,cluster2 := 'others']
  DE_cell_category_v_others <- DE_cell_category_v_others[order(statistic, decreasing = T)]
  DE_cell_category_v_others[, head(.SD, 10), by=cluster1]

  ## Add gene annotation and sort outputs
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description"))
  gene_annot <- gene_annot[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]

  DE_cell_category_v_others <- merge(DE_cell_category_v_others, gene_annot, by='gene', all.x=T)
  DE_cell_category_v_others <- DE_cell_category_v_others[order(cluster1, statistic, decreasing = T)]
  DE_cell_category_v_others[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  DE_cell_category_v_others <- DE_cell_category_v_others[order(rank, decreasing = T)]

  plot(DE_cell_category_v_others[cluster1=='Mammary Luminal mature' ,rank],
       DE_cell_category_v_others[cluster1=='Mammary Luminal mature' ,auc])

  plot(DE_cell_category_v_others[cluster1=='Mammary Luminal mature' ,auc],
       DE_cell_category_v_others[cluster1=='Mammary Luminal progenitor' ,auc])

  writexl::write_xlsx(x = DE_cell_category_v_others, path = file.path(spreadsheets.dir,'scRNA_cell_category_v_others_wilcox_de.xlsx'))
  DE_cell_category_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_cell_category_v_others_wilcox_de.xlsx')))

  table(DE_cell_category_v_others$cluster1)
  luminal_mature_genes <- DE_cell_category_v_others[cluster1=='Mammary Luminal mature' & p_val_adj < 0.01 & (auc >= .75), gene]
  luminal_progenitor_genes <- DE_cell_category_v_others[cluster1=='Mammary Luminal progenitor' & p_val_adj < 0.01 & (auc >= .75), gene]
  basal_genes <- DE_cell_category_v_others[cluster1=='Mammary Basal' & p_val_adj < 0.01 & (auc >= .75), gene]

  #### Single cell Basal ML LP classification ####
  refineLuminalBasalClassification()

  # How unique are luminal mature cells sample to sample?
  arid1a.multiome.luminalmature <- subset(arid1a.multiome, cell_classification=='Mammary Luminal mature')
  ( combos <- combn(x = unique(arid1a.multiome.luminalmature$bioNames), m = 2) )

  png(filename = file.path(plot.dir,"ML_UMAP.png"), width = 5, height = 7, units = 'in', bg = 'transparent', res = 200)
  FeaturePlot(object = arid1a.multiome.luminalmature, features = 'bioNames', slot = 'data') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., addLabels = F, label.size = 8,
                    pt.size = 0.5, axis.label = NA, legend.title = "Mature Luminal cells\nin samples",
                    show.legend = T, scale.color = F, pal = bionames_palette)
  dev.off()

  DEsample_one_v_all_luminalmature <- wilcoxauc(arid1a.multiome.luminalmature, group_by = 'bioNames', seurat_assay = 'RNA', assay = 'data')
  setDT(DEsample_one_v_all_luminalmature)
  setnames(DEsample_one_v_all_luminalmature,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEsample_one_v_all_luminalmature[,cluster2 := 'others']
  DEsample_one_v_all_luminalmature <- DEsample_one_v_all_luminalmature[order(statistic, decreasing = T)]
  DEsample_one_v_all_luminalmature[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  DEsample_one_v_all_luminalmature <- DEsample_one_v_all_luminalmature[order(rank, decreasing = T)]
  writexl::write_xlsx(x = DEsample_one_v_all_luminalmature, path = file.path(spreadsheets.dir,'scRNA_sample_luminal_mature_wilcox_de.xlsx'))

  aExp <- lapply(unique(arid1a.multiome.luminalmature$bioNames), function(b){
    DE_sample_ml <- DEsample_one_v_all_luminalmature[cluster1==b & p_val_adj < 0.05 & avg_logFC > 0.1,]
    unique(DE_sample_ml[,avgExpr, by=list(cluster1, gene)][order(avgExpr)])
  }) %>% rbindlist()
  aExp[,count:=c(1:.N), by=cluster1]
  aExp$cluster1Int <- factor(aExp$cluster1, levels = c(1:4))
  aExp$cluster1Int <- as.integer(aExp$cluster1)

  ggplot(aExp[,.(count=max(count)), by=cluster1], aes(x=count, y=cluster1, fill=cluster1)) +
    geom_col() + scale_fill_manual("", values = bionames_palette) +
    theme_classic(15) + theme(plot.background = element_blank(), panel.background = element_blank(), legend.position = 'none') +
    xlab("# D.E genes (FDR≤0.05 & logFC >0.1)") + ylab("")

  aExp$cluster1Int <- as.integer(aExp$cluster1)
  aExp$scale <- (aExp$avgExpr-min(aExp$avgExpr))/max(aExp$avgExpr)
  aExp <- aExp[order(cluster1, count)]
  png(filename = file.path(plot.dir,'Seurat_DE_Sample_vs_ALL_matureluminals.png'), width = 7, height = 3, bg = 'transparent', units = 'in', res = 300)
  ggplot(aExp, aes(y=cluster1Int-0.4, yend = cluster1Int+0.4, x=count, xend=count, color=cluster1, alpha=avgExpr, labels=cluster1)) +
    geom_segment(size = .2) +
    scale_alpha_continuous("µ expression",range = c(0.2, 1)) +
    guides(alpha=guide_legend(override.aes=list(size=c(2, 2, 2)))) +
    scale_y_continuous(expand = c(0,0), breaks = 1:4, labels = levels(aExp$cluster1)) +
    scale_color_manual("Samples", values = bionames_palette, guide='none') +
    theme_classic(18) + theme(text = element_text(colour = 'black'), plot.background = element_blank(),
                              panel.background = element_blank(),
                              legend.position = 'right', legend.background = element_blank()) +
    xlab("# D.E genes (FDR≤0.05 & logFC >0.1)") + ylab("")
  dev.off()

  topn <- 50
  degenes <- sapply(unique(arid1a.multiome.luminalmature$bioNames), function(b){
    DE_sample_ml <- DEsample_one_v_all_luminalmature[cluster1==b & p_val_adj < 0.05 & auc >= 0.7,]
    x <- DE_sample_ml[,head(gene, n=topn), by=cluster1]$V1
  })
  degenes <- unique(unlist(degenes))
  gsample <- DoHeatmapEnhance(object = arid1a.multiome.luminalmature,
                              features = degenes,
                              slot = 'scale.data',
                              assay = 'SCT',
                              group.by = 'bioNames',
                              disp.heatmap.vertically = F,
                              disp.identity.legend = F,
                              feature.font.size = 9,
                              legend.name = paste0("Mature luminal\ntop ",topn,"D.E."),
                              group.colors = bionames_palette
                              )

  pdf(file = file.path(plot.dir,'Seurat_DE_Sample_vs_ALL_matureluminals_label.pdf'), width = 12, height = 3, units = 'in', res = 150)
  gsample
  dev.off()

  # How unique are luminal cells sample to sample?
  arid1a.multiome.basal <- subset(arid1a.multiome, luminal_basal=='Basal')
  DEsample_one_v_all_basal <- wilcoxauc(arid1a.multiome.basal, group_by = 'bioNames', seurat_assay = 'RNA', assay = 'data')
  setDT(DEsample_one_v_all_basal)
  setnames(DEsample_one_v_all_basal,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEsample_one_v_all_basal[,cluster2 := 'others']
  DEsample_one_v_all_basal[order(statistic, decreasing = T)]
  plots <- lapply(unique(arid1a.multiome.basal$bioNames), function(b){
    DE_sample_basals <- DEsample_one_v_all_basal[cluster1==b & p_val_adj < 0.05 & auc >= 0.7,]
    x <- DE_sample_basals[,head(gene, n=50), by=cluster1]$V1
    gsample <- DoHeatmapEnhance(object = arid1a.multiome.basal,
                                features = x,
                                slot = 'scale.data',
                                assay = 'SCT',
                                group.by = 'bioNames',
                                size = 4,
                                angle = 90,
                                disp.identity.legend = F,
                                fill.legend.location = 'bottom',
                                group.colors = sample_palette, combine = T)
    gsample <- gsample + scale_fill_gradientn(colours = rev(pals::warmcool(10))) + ggtitle(paste(b,"basal up"))
  })
  png(filename = file.path(plot.dir,'Seurat_DE_Sample_vs_all_basals.png'),
      width = 4*6, height = 12, units = 'in', res = 150, bg = 'transparent')
  do.call(grid.arrange, c(plots, nrow=1))
  dev.off()

  DEsample_one_v_one_luminal <- lapply( seq_len(ncol(combos)) , function(i){
    dt <- wilcoxauc(arid1a.multiome.luminal, group_by = 'bioNames', groups_use = combos[,i], seurat_assay = 'RNA', assay = 'data') %>% setDT()
    setnames(dt,
             old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
             new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
    #dt[,cluster2 := 'others']
    label <- gsub(" ","_",paste(combos[,i], collapse = '_vs_'))
    dt[,label := label]
    dt[order(statistic, decreasing = T)]
  })
  names(DEsample_one_v_one_luminal) <- apply(combos, 2, function(x){ gsub(" ","_",paste(x[1],x[2],sep="_vs_")) })

  DEsample_one_v_one_luminal$Empty_Vector_vs_Foxa1_G275X[, head(.SD, 10), by=cluster1]
  DEsample_one_v_one_luminal$Foxa1_G275X_vs_Foxa1_WildType[, head(.SD, 10), by=cluster1]

  sapply( names(DEsample_one_v_one_luminal), function(x) { length(unique(DEsample_one_v_one_luminal[[x]][p_val_adj < 0.05, gene])) } )
  DE_sample_luminals <- rbindlist(lapply( DEsample_one_v_one_luminal, function(x) { x[p_val_adj < 0.05 & auc >= 0.7,] } ))

  plots <- lapply(unique(arid1a.multiome.luminal$bioNames), function(b){
    DE_sample_luminals <- rbindlist(lapply( DEsample_one_v_one_luminal, function(x) { x[cluster1==b & p_val_adj < 0.05 & auc >= 0.7,] } ))
    x <- DE_sample_luminals[,head(gene, n=20), by=label]$V1
    gsample <- DoHeatmapEnhance(object = arid1a.multiome.luminal,
                                features = x,
                                slot = 'scale.data',
                                assay = 'SCT',
                                group.by = 'bioNames',
                                size = 4,
                                angle = 90,
                                disp.identity.legend = F,
                                fill.legend.location = 'bottom',
                                group.colors = sample_palette, combine = T)
    gsample <- gsample + scale_fill_gradientn(colours = rev(pals::warmcool(10))) + ggtitle(paste(b,"luminal up"))
  })
  png(filename = file.path(plot.dir,'Seurat_DE_Sample_Luminals.png'),
      width = 4*6, height = 12, units = 'in', res = 150, bg = 'transparent')
  do.call(grid.arrange, c(plots, nrow=1))
  dev.off()

  # clusters
  grep("SCT_snn_res", colnames(arid1a.multiome@meta.data), value = T)
  Idents(arid1a.multiome) #SCT_snn_res_1.0
  identical(arid1a.multiome@meta.data$seurat_clusters, arid1a.multiome@meta.data[['SCT_snn_res.1']])
  DEcluster_v_others <- wilcoxauc(arid1a.multiome, group_by = 'seurat_clusters', seurat_assay = 'RNA', assay = 'data')
  DEcluster_v_others <- setDT(DEcluster_v_others)
  setnames(DEcluster_v_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEcluster_v_others[,cluster2 := 'others']
  DEcluster_v_others <- DEcluster_v_others[order(statistic, decreasing = T)]
  DEcluster_v_others[p_val_adj < 0.05, .N, by=cluster1]

  DEcluster_v_others[, head(.SD, 15), by=cluster1]

  #### DE Conserved Luminal vs Basal regardless cluster ####
  ## Regardless of cluster which genes are D.E. in luminal vs basal 'like' cells
  # markers that are conserved i.e. genes that are d.e. between cluster and others irrespective of luminal or basal identity in cluster1
  table(arid1a.multiome@meta.data$cluster_luminal_basal, arid1a.multiome@meta.data$seurat_clusters)

  lb.markers <- lapply(unique(arid1a.multiome$luminal_basal), function(i){
    dt <- FindConservedMarkers(object = arid1a.multiome, ident.1 = i, grouping.var = 'seurat_clusters')
    setDT(dt, keep.rownames = 'genes')
    dt[,luminal_basal := i]
  })
  head(lb.markers, 10)

  dt.markers <- rbindlist(lb.markers)
  setnames(dt.markers, 'cluster', 'cluster1')
  dt.markers[, cluster2 := 'others']
  features <- c(head(dt.markers[cluster==1], 10)$genes,'Sox9')
  FeaturePlot(arid1a.multiome, features = features, min.cutoff = "q9")
  writexl::write_xlsx(x = dt.markers, path = file.path(spreadsheets.dir,'scRNA_conserved_cluster_vs_others_regardless_luminal_basal_de.xlsx'))

  #### DE Conserved regardless Luminal vs Basal ####
  ## Regardless of luminal or basal cell label which genes are D.E. in clusters
  # markers that are conserved i.e. genes that are d.e. between cluster and others irrespective of luminal or basal identity in cluster1
  table(arid1a.multiome@meta.data$cluster_luminal_basal, arid1a.multiome@meta.data$seurat_clusters)

  lb.markers <- lapply(unique(arid1a.multiome$seurat_clusters), function(i){
    dt <- FindConservedMarkers(object = arid1a.multiome, ident.1 = i, grouping.var = 'luminal_basal')
    setDT(dt, keep.rownames = 'genes')
    dt[,cluster := i]
  }) %>% set_names(unique(arid1a.multiome$seurat_clusters))
  head(lb.markers, 10)

  dt.markers <- rbindlist(lb.markers)
  setnames(dt.markers, 'cluster', 'cluster1')
  dt.markers[, cluster2 := 'others']
  features <- c(head(dt.markers[cluster==1], 10)$genes,'Sox9')
  FeaturePlot(arid1a.multiome, features = features, min.cutoff = "q9")
  writexl::write_xlsx(x = dt.markers, path = file.path(spreadsheets.dir,'scRNA_conserved_cluster_vs_others_regardless_luminal_basal_de.xlsx'))

  # clusters
  library(dplyr)
  gvs <- lapply(sort(as.integer(unique(dt.markers$cluster1))), function(x){
    dt <- dt.markers[cluster1==x][order(Luminal_p_val, decreasing = F)]
    topn <- 20
    top_markers <- rbindlist(
      list(arrange(dt, rev(minimump_p_val)) %>% # Log2FC Up
             filter(rank( Luminal_avg_logFC, ties.method = "first") <= topn) %>%
             ungroup() %>% setDT,
           arrange(dt, minimump_p_val) %>% # Log2FC Down
             filter(rank( -Luminal_avg_logFC, ties.method = "first") <= topn) %>%
             ungroup() %>% setDT))

    xlim <- c(-max(abs(dt$avg_logFC)), max(abs(dt$avg_logFC)))
    ylim <- c(0, ceiling(max(dt$auc)))
    p1 <- ggplot(dt, aes(avg_logFC, pmax(-0.5, auc)) ) +
      geom_point(aes(color = p_val_adj < 0.05), pch=20) +
      scale_color_manual("Significant", values = c("gray60","firebrick1")) +
      scale_y_continuous(limits = ylim) +
      scale_x_continuous(limits = xlim) +
      geom_point(data = top_markers, color = "navy", shape=1) +
      ggtitle(paste('Cluster',x,'vs others')) +
      geom_text_repel(data = top_markers[auc >= .5], mapping = aes(label = gene),
                      size          =2.8,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = 0.1,
                      direction    = "x",
                      angle        = 90,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      max.time = 1) +
      geom_text_repel(data = top_markers[auc < .5], mapping = aes(label = gene),
                      size          =2.8,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = -0.1,
                      direction    = "x",
                      angle        = 90,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      max.time = 1) +
      theme_bw(9) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
      ylab("AUC") + xlab("mean log2 fold-change")
  })
  gvs[[4]]
  png(filename = file.path(plot.dir,'Cluster_volcano.png'), width = 26, height = 16, units = 'in', res=150, bg='transparent')
  do.call(grid.arrange, c(gvs, nrow=4))
  dev.off()

  #### Porportion of Luminal/Basal proportion vs gene expression ####
  genesTested <- c("Tcf7l1","Hnf4a","Hif3a","Gata3","Foxj1","Arid5a","Pou2f1","Trp63")

  p <- lapply(genesTested, function(g){
    dt <- data.table(arid1a.multiome@meta.data[,c("luminal_basal","bioNames")], keep.rownames = 'cells')
    ds <- as.data.table(arid1a.multiome@assays$RNA@data[g,], keep.rownames = 'cells')
    setnames(ds, c('cells',g))
    m <- merge(dt, ds, by='cells')

    g1 <- ggplot(m, aes(x=bioNames, y=get(g), color=luminal_basal)) +
      #geom_histogram(aes(y=..density..), colour="grey80", fill="white")+
      #geom_density(alpha=.6, size=1) +
      #geom_violin(trim = F, scale = 'width') +
      geom_boxplot() + ylab(paste("expression level")) + ggtitle(paste(g,"expression")) +
      xlab("bioNames") +
      scale_color_manual("Luminal/Basal", values = c('#1f77b4','#FFD966') ) + theme_classic(10)
    g1
  })
  do.call(grid.arrange, c(p, nrow=8))

  msfilt <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'BioNames_Motif_Heatmap_Vierstra_Filt.xlsx')))
  gplots <- lapply(sort(unique(DEluminal_v_basal_sample$bioNames)), function(i){
    genes <- unique(msfilt[bioNames==i & plot==T, geneCorrelated])
    genesTested <- ""
    if(i=='Foxa1 G275X'){
      genes <- unique(c(genesTested, genes))
      genesTested <- c("Tcf7l1","Hnf4a","Hif3a","Gata3","Foxj1","Arid5a","Pou2f1","Trp63")
    }
    toplot <- DEluminal_v_basal_sample[cluster1=='Luminal' & bioNames==i]
    toplot[,percent_diff := percent_cluster1 - percent_cluster2]
    p_cut_high <- quantile(toplot$percent_diff, .999)
    p_cut_low <- quantile(toplot$percent_diff, .001)
    ggplot(toplot, aes(x=percent_cluster1, y=percent_cluster2)) +
      geom_point(data = toplot[p_val_adj > 0.01], color='gray80', size=1) +
      geom_point(data = toplot[p_val_adj < 0.01], color='indianred', size=.9, shape=21) +
      geom_point(data = toplot[p_val_adj < 0.01 & abs(percent_diff) > p_cut_high], fill='red', shape=21, size=.9) +
      geom_point(data = toplot[p_val_adj < 0.01 & gene %in% genes & !gene %in% genesTested], color='blue', shape=21, size=.9) +
      geom_point(data = toplot[p_val_adj < 0.01 & gene %in% genesTested], color='darkgreen', shape=21, size=.9) +
      geom_text_repel(data = toplot[p_val_adj < 0.01 & percent_diff > p_cut_high & !gene %in% genes],
                      mapping = aes(label = gene),
                      size          =2.2,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = -10,
                      direction    = "y",
                      angle        = 0,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      force        = 0.5,
                      max.time = 1) +
      geom_text_repel(data = toplot[p_val_adj < 0.01 & percent_diff < p_cut_low & !gene %in% genes],
                      mapping = aes(label = gene),
                      size          =2.2,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = 10,
                      direction    = "y",
                      angle        = 0,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      force        = 0.5,
                      max.time = 1) +
      geom_text_repel(data = toplot[gene %in% genes & !gene %in% genesTested & percent_diff > 0],
                      mapping = aes(label = gene),
                      color        ='blue',
                      size          =2.5,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = -10,
                      direction    = "y",
                      angle        = 0,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      force        = 0.5,
                      max.time = 1) +
      geom_text_repel(data = toplot[gene %in% genes & !gene %in% genesTested & percent_diff < 0],
                      mapping = aes(label = gene),
                      color        ='blue',
                      size          =2.5,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = 10,
                      direction    = "y",
                      angle        = 0,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      force        = 0.5,
                      max.time = 1) +
      geom_text_repel(data = toplot[gene %in% genesTested & percent_diff > 0],
                      mapping = aes(label = gene),
                      color        ='darkgreen',
                      size          =2.5,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = -10,
                      direction    = "y",
                      angle        = 0,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      force        = 0.5,
                      max.time = 1) +
      geom_text_repel(data = toplot[gene %in% genesTested & percent_diff < 0],
                      mapping = aes(label = gene),
                      color        ='darkgreen',
                      size          =2.5,
                      force_pull   = 0, # do not pull toward data points
                      nudge_y      = 10,
                      direction    = "y",
                      angle        = 0,
                      hjust        = 0,
                      segment.size = 0.1,
                      max.overlaps = Inf,
                      max.iter = 1e4,
                      force        = 0.5,
                      max.time = 1) +
      geom_abline(size=.1) +
      xlab("Luminal cells (%)") + ylab("Basal cells (%)") + theme_bw() + facet_wrap(~bioNames) +
      theme(strip.background = element_blank(), strip.text = element_text(face = 'bold'), panel.grid = element_line(size=.1))
  })
  png(filename = file.path(plot.dir,"scRNA_multiome_luminal_basal_scatter_candidates.png"),
      width = 12, height = 8, units = 'in', bg = 'transparent', res = 150)
  do.call(grid.arrange, c(gplots, nrow=2))
  dev.off()

  #### Export object ####
  ExportToSingleCell(obj = arid1a.multiome)

  #### ATAC metadata ####
  dt <- data.table(data.frame(proj@cellColData), keep.rownames = "cells")
  arid1a.multiome@meta.data$atac_clusters <- ""
  arid1a.multiome@meta.data[dt$gex_cellid,]$atac_clusters <- dt$Clusters

  glabel <- DimPlot(arid1a.multiome, label = TRUE, group.by = 'seurat_clusters')
  glabel <- UmapPlotEnhance(obj = arid1a.multiome, p1 = glabel, addLabels = T,
                            label.size = 4, axis.label = NA,
                            pal=NA, title = "clusters", legend.title = '')

  g1 <- DimPlot(arid1a.multiome, label = F, group.by = 'atac_clusters') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., label.size = 10, pt.size = 0.9,
                    axis.label = NA, legend.title = '', title='snRNA samples',
                    pal=cluster_palette, panel.border.size = 0)
  pdf(file = file.path(plot.dir,'ClustersATAC_umap_basic.pdf'), width = 5, height = 2.6, bg = 'transparent')
  g1 + glabel$layers
  dev.off()

}

#### Analysis functions ####

#' Export a Seurat object to a SingleCellExperiment for use in scanpy
#'
#' Stores SCT counts/data and the RNA counts/data as assays together with the
#' PCA/UMAP embeddings, so the object can be converted to AnnData (h5ad).
#'
#' @param obj Seurat object (SCT + RNA assays, with `pca` and `umap` reductions).
#' @param out.file Path of the .rds file to write.
#' @param test If TRUE, inspect the object instead of writing `out.file`.
#' @return Nothing; the SingleCellExperiment is saved to `out.file`.
ExportToSingleCell <- function(obj,
                               out.file=file.path(robj.dir,'arid1a_sctransform_multiome_regress_sce.rds'),
                               test=F){
  # Just tests what a single cell experiment looks like
  # Will import Seurat object directly to scanpy
  # library(SeuratDisk) - does not work no longer developed
  # SaveH5Seurat(erg.scRNA, filename = file.path(robj.dir,'erg_sctransform_multiome_regress_adata.h5Seurat'))
  # Convert(file.path(robj.dir,'erg_sctransform_multiome_regress_adata.h5Seurat'), dest = "h5ad")

  DefaultAssay(obj) <- "RNA"

  #### Seurat way ####

  #### My way ####
  library(SingleCellExperiment)

  # SCT@counts -> adata.X corrected UMI counts
  # adata.layers['data'] -> SCT@data : log-normalized corrected UMI counts (Feature/Violin)
  # adata.layers['normcounts'] -> SCT@scale.data : pearson residuals (Heatmaps)
  # adata.layers['rnacounts'] -> RNA@counts : raw counts
  # adata.layers['rnadata'] -> RNA@data : log-normalized counts (non SCT) (Feature/Violin)
  # Corrected Pearson residuals for joint analysis, like PCA and clustering.
  # For DE gene calling, original SCT assay (or even the original RNA assay).
  # Heat maps: scale.data slot
  # Feature/Violin plots data slot
  # Never use integrated/corrected values when performing DE.
  # https://github.com/satijalab/seurat/issues/2899

  # Get the total variance:
  total_variance <- sum(matrixStats::rowVars(
    Seurat::GetAssayData(obj, assay = "RNA", slot = "scale.data")
  ))
  eigValues = (obj[["pca"]]@stdev)^2  ## EigenValues
  varExplained = eigValues / total_variance

  adata <- SingleCellExperiment(assays = list(counts = GetAssayData(object = obj, layer = "counts"),
                                              data = GetAssayData(object = obj, layer = "data"), #logcounts
                                              #normcounts = GetAssayData(object = obj, slot = "scale.data"),
                                              rnacounts = GetAssayData(object = obj, assay='RNA', layer = "counts"),
                                              rnadata = GetAssayData(object = obj, assay='RNA', layer = "data")), #rnalogcounts
                                colData=DataFrame(obj@meta.data)) #,rowData=DataFrame(GetAssayData(object = obj, slot='meta.features')))
  # creates obsm
  reducedDims(adata) <- list(PCA=obj@reductions[['pca']]@cell.embeddings,
                             UMAP=obj@reductions[['umap']]@cell.embeddings,
                             #KNNDistances=Distances(obj@graphs[["SCT_nn"]]),
                             neighbors=as(obj@graphs[['SCT_nn']], "dgCMatrix"), # should go into UNS
                             snn=as(obj@graphs[['SCT_snn']], "dgCMatrix"))
                             #knn=as(obj@graphs[["SCT_nn"]]@nn.idx, "dgCMatrix")) # Gives cell with self
  # creates uns
  metadata(adata) <- list(PCs=Loadings(obj[['pca']]),
                          highly_variable = DataFrame('highly_variable'=obj@assays$RNA@var.features),
                          variance=eigValues,
                          variance_ratio=varExplained)

  (factors <- which(sapply(adata@colData, is.factor)))
  for(i in factors){
    adata@colData[,i] <- as.numeric(adata@colData[,i])
  }
  which(sapply(adata@colData, is.factor))
  if(!is.null(out.file)){
    message("exporting data to sce format: ",out.file)
    saveRDS(object = adata, out.file)
  }

  ## Test to Read in
  if(test){
    obj <- readRDS(out.file)

    adata <- SingleCellExperiment(assays = list(counts = assays(obj)$counts,
                                                data = assays(obj)$data,
                                                rnacounts = assays(obj)$rnacounts,
                                                rnadata = assays(obj)$rnadata),
                                  colData=DataFrame(obj@colData))
    # creates obsm
    reducedDims(adata) <- list(PCA=reducedDims(obj)$PCA,
                               UMAP=reducedDims(obj)$UMAP,
                               neighbors=as(reducedDims(obj)$neighbors, "dgCMatrix"), # should go into UNS
                               snn=as(reducedDims(obj)$snn, "dgCMatrix")) # should go into UNS

    # creates uns
    metadata(adata) <- list(PCs=metadata(obj)$PCs,
                            highly_variable = DataFrame('highly_variable'=metadata(obj)$highly_variable),
                            variance=metadata(obj)$variance,
                            variance_ratio=metadata(obj)$variance_ratio)

  }

  adata

}

#' Import scanpy/palantir reductions back into the Seurat object
#'
#' Reads the diffusion map and force-directed layout computed in scanpy
#' (`<data.dir>/scanpy/scRNA_diffmap.txt`, `scRNA_fa.txt`), adds them as the
#' `dc` and `fa` reductions of `arid1a.multiome`, plots them and saves the object.
ScanpyImportReductions <- function(){

  #### Diffusion Map ####
  # to process data in scanpy with palantir run ArchRProjectToSCE
  diffmap <- fread(file.path(data.dir, 'scanpy', 'scRNA_diffmap.txt'),
                   header = T, sep = "\t")
  diffmap <- data.frame(diffmap, row.names = diffmap$V1)
  diffmap$V1 <- NULL
  message("Removing 0 DC")
  diffmap$X0 <- NULL
  colnames(diffmap) <- gsub("X","IterativeLSI#DC_Dimension_",colnames(diffmap))
  ## outlier.cells <- c('EV-day5#AATTGCCAGTAGGCGG-1','WT-day1#GCGAAGTAGAGGAAGG-1')
  ## diffmap[rownames(diffmap) %in% outlier.cells,] <- NA
  dr <- CreateDimReducObject(embeddings = data.matrix(diffmap), key='DC')
  arid1a.multiome[["dc"]] <- dr

  #### Force Directed Map ####
  # to process data in scanpy with palantir run ArchRProjectToSCE
  fa <- fread(file.path(data.dir, 'scanpy', 'scRNA_fa.txt'),
              header = T, sep = "\t")
  fa <- data.frame(fa, row.names = fa$V1)
  fa$V1 <- NULL
  colnames(fa) <- gsub("X0","PCA#FA_Dimension_1",colnames(fa))
  colnames(fa) <- gsub("X1","PCA#FA_Dimension_2",colnames(fa))
  ## outlier.cells <- c('EV-day5#AATTGCCAGTAGGCGG-1','WT-day1#GCGAAGTAGAGGAAGG-1')
  ## fa <- fa[!rownames(fa) %in% outlier.cells,]
  ## fa[rownames(fa) %in% outlier.cells,] <- NA
  dr <- CreateDimReducObject(embeddings = data.matrix(fa), key='FA')
  arid1a.multiome[["fa"]] <- dr

  #### Plot Reductions ####
  reduction <- 'dc'
  reduction <- 'fa'
  coord1 <- paste0(toupper(reduction),'_1')
  coord2 <- paste0(toupper(reduction),'_2')

  p <- DimPlot(arid1a.multiome, label = F, reduction = reduction, group.by = 'bioNames') + NoLegend()
  g0 <- UmapPlotEnhance(obj = arid1a.multiome, p1 = p, label.size = 10,
                       pt.size = 1.1,
                       pal=bionames_palette, title = '',
                       coord1 = coord1, coord2 = coord2,
                       fixed_coord = T)

  df <- data.frame(arid1a.multiome@reductions[[reduction]]@cell.embeddings, arid1a.multiome$seurat_clusters)
  g1 <- UmapPlotEnhance(obj = arid1a.multiome, p1 = df,
                        label.size = 5, pt.size = 1,
                        coord1 = coord1, coord2 = coord2, addLabels = T,
                        pal=NA, title = "") + theme(legend.position = 'none')

  pdf(file = file.path(plot.dir,paste0(reduction,'_bionames.pdf')), width = 6, height = 6, bg = 'transparent')
  g0 + g1$layers
  dev.off()

  message("Saving seurat object to: ",arid1a.object)
  saveRDS(arid1a.multiome, file = arid1a.object)

}

#' Classify cells as basal / luminal progenitor / luminal mature and refine clusters
#'
#' Scores cells (AddModuleScore) with gene sets derived from the one-vs-rest
#' differential expression between luminal mature, luminal progenitor and basal
#' cells, clusters the scores (kmeans, then a multivariate Gaussian model),
#' draws ternary plots, and saves the classification in the object.
#' Uses the `arid1a.multiome` and `DE_cell_category_v_others` objects from main().
refineLuminalBasalClassification <- function(){
  library(MASS)
  library(plotly)
  #### Score cells from DE ####
  luminal_mature_genes <- DE_cell_category_v_others[cluster1=='Mammary Luminal mature' & p_val_adj < 0.01 & (auc >= .75), gene]
  luminal_progenitor_genes <- DE_cell_category_v_others[cluster1=='Mammary Luminal progenitor' & p_val_adj < 0.01 & (auc >= .75), gene]
  basal_genes <- DE_cell_category_v_others[cluster1=='Mammary Basal' & p_val_adj < 0.01 & (auc >= .75), gene]

  # Choose confident genes not overlapping betweens sets
  rem <- c(intersect(luminal_mature_genes, basal_genes),
           intersect(luminal_progenitor_genes, basal_genes),
           intersect(luminal_mature_genes, luminal_progenitor_genes)) %>% unique
  luminal_mature_genes <- setdiff(luminal_mature_genes, rem)
  luminal_progenitor_genes <- setdiff(luminal_progenitor_genes, rem)
  basal_genes <- setdiff(basal_genes, rem)

  arid1a.multiome <- AddModuleScore(object = arid1a.multiome,
                                   features = list(luminal_mature_genes, luminal_progenitor_genes, basal_genes),
                                   assay = 'RNA',
                                   name = 'Luminal_vs_basal',
                                   search = F) # searches HUGO containing human gene names only
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='Luminal_vs_basal1')] <- 'luminal_mature_gene_scores'
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='Luminal_vs_basal2')] <- 'luminal_progenitor_gene_scores'
  colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='Luminal_vs_basal3')] <- 'basal_gene_scores'

  dt <- data.table(arid1a.multiome@meta.data[,grep("luminal_mature_gene_scores|luminal_progenitor_gene_scores|basal_gene_scores|bioNames", colnames(arid1a.multiome@meta.data))])
  X <- as.matrix(dt[,list(luminal_mature_gene_scores,luminal_progenitor_gene_scores,basal_gene_scores)])
  rownames(X) <- rownames(arid1a.multiome@meta.data)
  ## Look at density
  library(MASS)

  #### Kmeans ####
  library(factoextra)
  nbc <- fviz_nbclust(x = X, FUNcluster = kmeans)
  nbc
  (cl <- kmeans(x = X, centers = 3))
  plot(X, col=cl$cluster)
  points(cl$centers, col=1:3, pch=8, cex=2)
  print(cl)
  # totss = total variance in data
  cl$totss
  # size of each cluster
  cl$size
  fviz_cluster(object = cl, data = X)

  #### MVGaussian ####
  library(mclust)
  mod <- densityMclust(X)
  summary(mod)
  plot(mod, what = "BIC")
  plot(mod, what = "density")

  BIC <- mclustBIC(X)
  mod1 <- Mclust(X, x=BIC)
  plot(mod1, what = "BIC")
  plot(mod1, what = "density")
  plot(mod1, what = "classification", col=as.vector(pals::alphabet()))
  arid1a.multiome$mod_classification <- as.character(mod1$classification)
  g1 <- FeaturePlot(arid1a.multiome, features = 'mod_classification', ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., pal = pals::alphabet(10), pt.size = 1, axis.label = "", addLabels = F, legend.title = 'class')

  mod9 <- Mclust(X, G = 9)
  plot(mod9, what = "classification", col=pals::alphabet())
  plot(mod9, what = "uncertainty", col=pals::alphabet())
  arid1a.multiome@meta.data$BIC_cell_classification <- ""
  arid1a.multiome@meta.data[names(mod9$classification),]$BIC_cell_classification <- mod9$classification
  arid1a.multiome$cell_classification_uncertainty <- mod9$uncertainty

  mod5 <- Mclust(X, G = 5)
  plot(mod5, what = "density")
  plot(mod5, what = "uncertainty", col=pals::alphabet())
  plot(mod5, what = "classification")
  uncerPlot(z = mod5$z)
  arid1a.multiome@meta.data[names(mod5$classification),]$BIC_cell_classification <- mod5$classification
  arid1a.multiome$cell_classification_uncertainty <- mod5$uncertainty

  mod4 <- Mclust(X, G = 4)
  plot(mod4, what = "density")
  plot(mod4, what = "uncertainty", col=pals::alphabet())
  plot(mod4, what = "classification")
  uncerPlot(z = mod4$z)
  arid1a.multiome@meta.data[names(mod4$classification),]$BIC_cell_classification <- mod4$classification

  mod3 <- Mclust(X, G = 3)
  summary(mod3)
  plot(mod3, what = "density")
  plot(mod3, what = "uncertainty", col=pals::alphabet())
  uncerPlot(z = mod3$z)
  arid1a.multiome@meta.data[names(mod3$classification),]$BIC_cell_classification <- mod3$classification
  arid1a.multiome$cell_classification_uncertainty <- mod3$uncertainty

  plot(mod3, col = topo.colors(50))
  plot(mod3, what = "classification")
  table(mod3$classification)
  plot(mod3, what = "uncertainty")
  plot(mod3, what = "density")

  dens <- densityMclust(X)
  table(arid1a.multiome$seurat_clusters, mod3$classification)
  # 1 luminal mature
  # 2 luminal progenitor
  # 3 basal
  dtclass <- data.table(data.frame(cell_classification_num=mod3$classification), keep.rownames = 'cells')
  dtclass[cell_classification_num==1, cell_classification := 'Mammary Luminal mature']
  dtclass[cell_classification_num==2, cell_classification := 'Mammary Luminal progenitor']
  dtclass[cell_classification_num==3, cell_classification := 'Mammary basal']
  table(dtclass[,list(cell_classification_num, cell_classification)])
  arid1a.multiome@meta.data$cell_classification <- ""
  arid1a.multiome@meta.data[dtclass$cells,]$cell_classification <- dtclass$cell_classification
  table(arid1a.multiome@meta.data$cell_classification)
  dtclassprob <- data.table(data.frame(mod3$z), keep.rownames = 'cells')
  setnames(dtclassprob, c("X1","X2","X3"), c("cell_class_prob_luminal_mature","cell_class_prob_luminal_progenitor","cell_class_prob_basal"))
  arid1a.multiome@meta.data$cell_class_prob_luminal_mature <- 0
  arid1a.multiome@meta.data[dtclassprob$cells,]$cell_class_prob_luminal_mature <- dtclassprob$cell_class_prob_luminal_mature
  arid1a.multiome@meta.data$cell_class_prob_luminal_progenitor <- 0
  arid1a.multiome@meta.data[dtclassprob$cells,]$cell_class_prob_luminal_progenitor <- dtclassprob$cell_class_prob_luminal_progenitor
  arid1a.multiome@meta.data$cell_class_prob_basal <- 0
  arid1a.multiome@meta.data[dtclassprob$cells,]$cell_class_prob_basal <- dtclassprob$cell_class_prob_basal
  if(identical(names(mod3$uncertainty), rownames(arid1a.multiome@meta.data))){
    arid1a.multiome$cell_classification_uncertainty <- mod3$uncertainty
  } else {
    warning("Problem with uncertainty")
  }

  dt <- data.table(arid1a.multiome@meta.data, keep.rownames = 'cells')

  # Bar classifications
  PlotBarVariableInClusters(seuratObj = arid1a.multiome, variable = 'cell_classification', palette = luminal_basal_palette,
                            legend.title = "Basal/LP/ML class", filename = "Basal_Luminal_BarPlot")

  # UMAP classified
  g1 <- FeaturePlot(arid1a.multiome, features = c('cell_classification'),
                    pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., addLabels = F, axis.label = "", pt.size = .8,
                    legend.title = 'Basal/LP/ML class', pal=luminal_basal_palette)
  # classification uncertainty
  g2 <- FeaturePlot(arid1a.multiome, features = c('cell_classification_uncertainty'),
                    pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = ., addLabels = F, axis.label = "", pt.size = .8,
                    scale.color = F, pal = pals::brewer.orrd(10),
                    legend.title = 'Basal/LP/ML uncertainty')
  # BIC
  g3 <- FeaturePlot(arid1a.multiome, features = c('BIC_cell_classification'),
                    pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    addLabels = F, label.size = 5,
                    axis.label = "", pt.size = .8,
                    legend.title = 'BIC',pal=as.vector(pals::alphabet(10)))
  png(filename = file.path(plot.dir,'Luminal_Basal_Classification_Umap.png'),
      width = 9, height = 9, units = 'in',
      bg = 'transparent', res = 200)
  #grid::grid.draw(rbind( ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3))) # left align plots
  grid.arrange(g1, g2, g3, ncol=1)
  dev.off()

  library(plotly)
  plot_ly(x=dt$luminal_progenitor_gene_scores, y=dt$luminal_mature_gene_scores, z=dt$basal_gene_scores, size = 5,
          type="scatter3d", mode="markers", color=dt$cell_classification_uncertainty) %>%
    layout(scene = list(xaxis=list(title = 'lp gene scores'),
                        yaxis=list(title = 'ml gene scores'),
                        zaxis=list(title = 'basal gene scores'),
                        legendgrouptitle=list(size=12)),
           font=list(
             family="Arial",
             size=18,
             color="RebeccaPurple"
           ),
           paper_bgcolor='rgba(0,0,0,0)',
           plot_bgcolor='rgba(0,0,0,0)')

  ### density classify by sample
  dt <- data.table(arid1a.multiome@meta.data)
  #dt$density <- dens$density
  minlim <- min(dt$luminal_mature_gene_scores, dt$luminal_progenitor_gene_scores, dt$basal_gene_scores)
  maxlim <- max(dt$luminal_mature_gene_scores, dt$luminal_progenitor_gene_scores, dt$basal_gene_scores)
  gp <- lapply(sort(unique(dt$bioNames)), function(x){
    idx <- which(dt$bioNames==x)
    set <- dt[idx]
    g0 <- ggplot(dt,
           aes(x=luminal_mature_gene_scores, y=luminal_progenitor_gene_scores)) +
      # fill density
      geom_density_2d_filled(data = set, contour_var = "ndensity") +
      # show rings
      geom_density_2d(bins=15, size=.2, color='gray50') +
      # fill colors
      scale_fill_manual(values = rev(pals::kovesi.linear_kryw_5_100_c64(10))) +
      ggtitle(x) +
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(16, units='points'),
            legend.text = element_text(size = 12), legend.position = 'none') +
      guides(colour = guide_legend(override.aes = list(size=4, alpha=1)))
    g1 <- ggplot(dt,
           aes(x=luminal_mature_gene_scores, y=basal_gene_scores)) +
      # fill density
      geom_density_2d_filled(data = set, contour_var = "ndensity") +
      # show rings
      geom_density_2d(bins=15, size=.2, color='gray50') +
      # fill colors
      scale_fill_manual(values = rev(pals::kovesi.linear_kryw_5_100_c64(10))) +
      ggtitle(x) +
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(16, units='points'),
            legend.text = element_text(size = 12), legend.position = 'none') +
      guides(colour = guide_legend(override.aes = list(size=4, alpha=1)))
    g2 <- ggplot(dt,
           aes(x=luminal_progenitor_gene_scores, y=basal_gene_scores)) +
      # fill density
      geom_density_2d_filled(data = set, contour_var = "ndensity") +
      # show rings
      geom_density_2d(bins=15, size=.2, color='gray50') +
      # fill colors
      scale_fill_manual(values = rev(pals::kovesi.linear_kryw_5_100_c64(10))) +
      ggtitle(x) +
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(16, units='points'),
            legend.text = element_text(size = 12), legend.position = 'none') +
      guides(colour = guide_legend(override.aes = list(size=4, alpha=1)))
    list(g0,g1,g2)
  })

  pdf(file = file.path(plot.dir,'Luminal_Basal_MVG_classification_by_sample.pdf'),
      width = 10, height = 3, bg = 'transparent', onefile = T)
  do.call(grid.arrange, c(gp[[1]], nrow=1))
  do.call(grid.arrange, c(gp[[2]], nrow=1))
  do.call(grid.arrange, c(gp[[3]], nrow=1))
  do.call(grid.arrange, c(gp[[4]], nrow=1))
  dev.off()

  #### Ternary plot ####
  library(Ternary)
  ternaryDat <- data.table(arid1a.multiome@meta.data[,c("bioNames",
                                                        "cell_class_prob_luminal_progenitor",
                                                        "cell_class_prob_luminal_mature",
                                                        "cell_class_prob_basal")], keep.rownames = 'cells')
  png(filename = file.path(plot.dir,'TernaryPlot.png'), width = 6, height = 6, units = 'in', bg = 'transparent', res = 200)
  TernaryPlot(alab = "Basal \u2192", blab = "Mature \u2192", clab = " \u2190 Progenitor",
              lab.col = luminal_basal_palette,
              point = 'up', grid.lines = 5,
              lab.cex = 1.2, grid.minor.lines = 0,
              grid.lty = 'dotted', #col = rgb(0.95, 0.95, 0.95), grid.col = 'white',
              axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
              #axis.rotate = FALSE,
              padding = 0.08)
  # Colour the background:
  #cols <- TernaryPointValues(rgb)
  #ColourTernary(cols, spectrum = NULL)
  # make points
  dat_col <- sapply(ternaryDat$bioNames, function(x) {bionames_palette[x]})
  ternaryDat <- round(ternaryDat[, .(cell_class_prob_basal, cell_class_prob_luminal_mature, cell_class_prob_luminal_progenitor)], digits = 8)
  #ColourTernary(TernaryDensity(ternaryDat, resolution = 10L))
  TernaryPoints(ternaryDat,
                cex = 0.8,
                col = 'black', #dat_col,
                pch = 16)
  #TernaryDensityContour(ternaryDat, resolution = 30L)
  dev.off()

  # Ternary of Arid1a +/+ only!
  ternaryDat <- data.table(arid1a.multiome@meta.data[,c("bioNames","condition",
                                                        "cell_class_prob_luminal_progenitor",
                                                        "cell_class_prob_luminal_mature",
                                                        "cell_class_prob_basal")], keep.rownames = 'cells')
  ternaryDat <- ternaryDat[condition=='WT']
  png(filename = file.path(plot.dir,'TernaryPlotWT.png'), width = 6, height = 6, units = 'in', bg = 'transparent', res = 200)
  TernaryPlot(alab = "Basal \u2192", blab = "Mature \u2192", clab = " \u2190 Progenitor",
              lab.col = luminal_basal_palette,
              point = 'up', grid.lines = 5,
              lab.cex = 1.2, grid.minor.lines = 0,
              grid.lty = 'dotted', #col = rgb(0.95, 0.95, 0.95), grid.col = 'white',
              axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
              #axis.rotate = FALSE,
              padding = 0.08)
  # Colour the background:
  #cols <- TernaryPointValues(rgb)
  #ColourTernary(cols, spectrum = NULL)
  # make points
  dat_col <- sapply(ternaryDat$bioNames, function(x) {bionames_palette[x]})
  ternaryDat <- round(ternaryDat[, .(cell_class_prob_basal, cell_class_prob_luminal_mature, cell_class_prob_luminal_progenitor)], digits = 8)
  #ColourTernary(TernaryDensity(ternaryDat, resolution = 10L))
  TernaryPoints(ternaryDat,
                cex = 0.8,
                col = 'black', #dat_col,
                pch = 16)
  #TernaryDensityContour(ternaryDat, resolution = 30L)
  dev.off()

  ### classify by cluster
  gp <- lapply(unique(dt$seurat_clusters), function(x){
    idx <- which(dt$seurat_clusters==x)
    set <- dt[idx]
    Density <- dens$density[idx]
    ggplot(dt,
           aes(x=Cd24.Scores, y=Cd49.Scores)) +
      geom_point(data = set, aes(color=Density)) +
      #stat_density_2d(data = dt[bioNames==x], aes(color=luminal_basal_classification), size=.5, alpha=0.6) +
      geom_density_2d(bins=15, size=.4, color='black') +
      #scale_color_viridis() +
      scale_color_gradientn(colours = pals::parula(50)) +
      #scale_color_manual(values = bionames_palette_list[[x]]) +
      ggtitle(x) +
      theme_bw(12) +
      theme(panel.grid = element_blank(), plot.background = element_blank(),
            legend.key.size = unit(24, units='points'),
            legend.text = element_text(size = 12))
  })
  png(filename = file.path(plot.dir,'Luminal_Basal_MVG_classification_by_sample.png'), width = 12, height = 8, units = 'in', res=200, bg = 'transparent')
  do.call(grid.arrange, c(gp, nrow=4))
  dev.off()

  dt <- data.table(foxa1.multiome@meta.data[,c('bioNames', 'seurat_clusters', 'luminal_basal_classification_uncertainty', 'luminal_basal_classification')])
  position=position_dodge(1)
  ggplot(dt, aes(x=luminal_basal_classification_uncertainty, color=bioNames)) +
    #geom_histogram(aes(y=..density..), position='identity', alpha=0.5) + #,  binwidth=.01) +
    geom_density(alpha=0.2) +
    #scale_fill_manual("", values = sample_palette) +
    scale_color_manual("", values = sample_palette) +
    theme_bw() + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(expand = F)

  ggplot(dt, aes(y=luminal_basal_classification_uncertainty, x=seurat_clusters, color=seurat_clusters)) +
    geom_col() +
    scale_fill_manual("", values = seurat_cluster_palette) +
    theme_bw() + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(ylim = c(0,1), expand = F)

  ggplot(dt, aes(x=luminal_basal_classification_uncertainty, color=luminal_basal_classification, fill=luminal_basal_classification)) +
    geom_density(alpha=0.5) +
    xlab("Classification uncertainty") +
    scale_fill_manual("", values = luminal_basal_palette) +
    scale_color_manual("", values = luminal_basal_palette) +
    theme_bw() + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(expand = F)

  message("Saving seurat object to: ",arid1a.object)
  saveRDS(arid1a.multiome, file = arid1a.object)

}

#' Convert the Saeki et al. mammary cluster markers to a GMT file
#'
#' Reads every sheet of Mammary_Significant_GeneMarkers_by_Cluster_Saeki_et_al.xlsx
#' (in `data.dir`) and writes the gene sets to `saeki.gmt` in `annot.dir`. Then
#' scores the cells on each marker set (AddModuleScore), derives a Saeki-based
#' mammary cell classification, refines it, and saves the object.
SaekiMammaryMarkers <- function(){

  s <- readxl::excel_sheets(file.path(data.dir,'Mammary_Significant_GeneMarkers_by_Cluster_Saeki_et_al.xlsx'))
  geneLists <- lapply(s, function(i){
    wb <- readxl::read_xlsx(file.path(data.dir,'Mammary_Significant_GeneMarkers_by_Cluster_Saeki_et_al.xlsx'),
                            sheet = i, col_names = T)
    setDT(wb)
  }) %>% set_names(s)

  writeGMT(object = geneLists, fname = file.path(annot.dir, 'saeki.gmt'))

  colnames(arid1a.multiome@meta.data)
  x <- grep("gene", invert = T, value=T, grep("mammary_scores",colnames(arid1a.multiome@meta.data), value=T))
  for(i in x)
    arid1a.multiome@meta.data[[i]] <- NULL
  colnames(arid1a.multiome@meta.data)

  for(i in names(geneLists)){
    message("Scoring list: ",i)
    arid1a.multiome <- AddModuleScore(object = arid1a.multiome,
                                      features = head(geneLists[[i]],100),
                                      assay = 'RNA',
                                      name = 'Saeki_scores',
                                      search = F) # searches HUGO containing human gene names only
    colnames(arid1a.multiome@meta.data)[which(colnames(arid1a.multiome@meta.data)=='Saeki_scores1')] <- paste0(i,'_mammary_scores')
  }
  grep("Klk", rownames(arid1a.multiome), value = T)

  colnames(arid1a.multiome@meta.data)

  gp <- lapply( names(saeki_mammary_scores), function(i){
    p <- DimPlot(arid1a.multiome, label = F, group.by = i) + NoLegend()
    p$data[[i]] <- as.numeric(as.character(p$data[[i]]))
    g0 <- UmapPlotEnhance(obj = arid1a.multiome,
                          p1 = p,
                          panel_width_in = 2.5,
                          label.size = 10,
                          pt.size = 1.1,
                          legend.height.factor = .8,
                          pal = blue_gray_purple_diverging_palette, # pals::magma(10),
                          title = saeki_mammary_scores[[i]],
                          break.unit = 0.2,
                          fixed_coord = T)
    df <- data.frame(arid1a.multiome@reductions$umap@cell.embeddings,
                     seurat_clusters=arid1a.multiome$seurat_clusters)
    df$seurat_clusters <- as.numeric(as.character(df$seurat_clusters))
    g1 <- UmapPlotEnhance(obj = arid1a.multiome,
                          p1 = df,
                          panel_width_in = 2.5,
                          label.size = 3,
                          pt.size = 1,
                          addLabels = T,
                          pal=NA, title = saeki_mammary_scores[[i]]) + theme(legend.position = 'none')
    g0 + g1$layers
  }) %>% set_names(names(saeki_mammary_scores))

  png(filename = file.path(plot.dir,paste0('umap_mammary_scores_grid.png')),
      width = 9,
      height = 8,
      units = 'in',
      res = 150)
  do.call(grid.arrange, c(gp, nrow=3, top='snRNA'))
  dev.off()

  #### Violins of Saeki_scores ####
  f <- c("S.Score","G2M.Score", names(saeki_mammary_scores))
  g <- ViolinPlotEnhance(arid1a.multiome,
                         features = f,
                         assay = 'RNA',
                         slot = data,
                         pt.size = 0,
                         flip = T,
                         group.by = 'seurat_clusters',
                         stack = T,
                         fill.by = 'ident',
                         title = "Cell cycle and mammary scores",
                         y_axis_label = "Normalized scores",
                         x_axis_text_angle = 90,
                         strip.text.size = 14,
                         pal = pals::brewer.pubu(25),
                         #pal = pals::coolwarm(25)
                         ) + xlab("")
  pdf(file = file.path(plot.dir,"cell_cycle_saeki_mammary_violin.pdf"), width = 6.5, height = 6)
  g
  dev.off()

  dt <- data.table(arid1a.multiome@meta.data, keep.rownames = "cells")
  dt[, lapply(.SD, mean), .SDcols = names(saeki_mammary_scores), by=seurat_clusters]

  #### Create Saeki classification ####

  arid1a.multiome@meta.data$mammary_gland_cluster_classification <- ""
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==1, 'mammary_gland_cluster_classification'] <- 'L-Hor'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==2, 'mammary_gland_cluster_classification'] <- 'L-Hor'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==3, 'mammary_gland_cluster_classification'] <- 'L-Alv'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==4, 'mammary_gland_cluster_classification'] <- 'L-Alv'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==5, 'mammary_gland_cluster_classification'] <- 'L-Hor'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==6, 'mammary_gland_cluster_classification'] <- 'L-Hor'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==7, 'mammary_gland_cluster_classification'] <- 'L-Alv'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==8, 'mammary_gland_cluster_classification'] <- 'L-Alv'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==9, 'mammary_gland_cluster_classification'] <- 'L-Hor'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==10, 'mammary_gland_cluster_classification'] <- 'L-Hor'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==11, 'mammary_gland_cluster_classification'] <- 'LH-Pro'
  arid1a.multiome@meta.data[arid1a.multiome@meta.data$seurat_clusters==12, 'mammary_gland_cluster_classification'] <- 'Basal/MaSC'

  #### Refine Saeki classes ####
  arid1a.multiome@meta.data[arid1a.multiome$bioNames %in% c("Arid1a(-/-)","Arid1a(-/-) E2") &
                              arid1a.multiome$mammary_gland_cluster_classification=="L-Hor" ,]$mammary_gland_refine_clusters <- "L-Hor-Arid1a-"

  table(arid1a.multiome$mammary_gland_cluster_classification, arid1a.multiome$seurat_clusters)
  ## MaSC is difficult - density of MaSC scores in clusters ##
  dt <- data.table(arid1a.multiome@meta.data, keep.rownames = 'cells')
  ggplot(dt, aes(x=masc_mammary_scores, y=basal_mammary_scores, color=seurat_clusters)) +
    geom_point_rast() + theme_bw() + facet_wrap(~seurat_clusters, nrow=3)

  dt[order(masc_mammary_scores, decreasing=T)][1:10]
  # MMO1_Control_multiome#CGGAATCGTTAGTACG-1

  coi <- grep("mammary_scores|seurat_clusters", colnames(dt), value=F)
  dmelt <- melt(dt[,coi, with=F], id.vars = 'seurat_clusters')
  png(filename = file.path(plot.dir,"cluster_mammaryscores_violin.png"), width = 8, height = 6,
      bg = 'transparent', units = 'in', res=150)
  ggplot(dmelt, aes(x=seurat_clusters, y=value, fill=variable)) +
    #facet_wrap(~seurat_clusters, scales = 'free_x', nrow=3) +
    facet_wrap(~variable, scales = 'free_x', nrow=3) +
    ylab("Scores") + xlab("Clusters") +
    scale_fill_manual("mammary_score", values = unname(mammary_cluster_class_palette)) +
    geom_violin(scale = 'width') + theme_bw(16) + theme(legend.position = 'none')
  dev.off()

  reduction <- 'umap'
  coord1 <- paste0(toupper(reduction),'_1')
  coord2 <- paste0(toupper(reduction),'_2')
  p <- DimPlot(arid1a.multiome, label = F, reduction = reduction,
               group.by = 'mammary_gland_cluster_classification')
  g0 <- UmapPlotEnhance(obj = arid1a.multiome,
                        p1 = p,
                        label.size = 10,
                        pt.size = .5,
                        pal = mammary_cluster_class_palette_list,
                        title = '',
                        legend.title = 'snRNA',
                        legend.title.rotation = 0,
                        coord1 = coord1,
                        coord2 = coord2,
                        fixed_coord = T)

  df <- data.frame(arid1a.multiome@reductions[[reduction]]@cell.embeddings, arid1a.multiome$seurat_clusters)
  g1 <- UmapPlotEnhance(obj = arid1a.multiome,
                        p1 = df,
                        label.size = 3,
                        coord1 = coord1, coord2 = coord2,
                        addLabels = T,
                        pal=NA,
                        title="")

  png(filename = file.path(plot.dir,paste0(reduction,'_mammary_cluster_classification_single.png')),
      width = 5,
      height = 3,
      units = 'in', res = 250)
  g0 + g1$layers
  dev.off()

  message("Saving seurat object to: ",arid1a.object)
  saveRDS(arid1a.multiome, file = arid1a.object)
}

#' Sum expression per identity class
#'
#' Modelled on Seurat's AverageExpression() but sums (rather than averages) the
#' values of the chosen slot, which is what is needed to create pseudo-bulk
#' count tables (e.g. for DESeq2).
#'
#' @param object Seurat object.
#' @param assays Assays to use (default: all).
#' @param features Features to include (default: all).
#' @param slot Slot to sum ("counts", "data" or "scale.data").
#' @return A list with one matrix (genes x identities) per assay.
TotalExpression <- function (object, assays = NULL, features = NULL,
                             return.seurat = FALSE,
                             add.ident = NULL, slot = "data",
                             use.scale = FALSE, use.counts = FALSE,
                             verbose = TRUE, ...)
{
  library(rlang)
  Seurat:::CheckDots(..., fxns = "CreateSeuratObject")
  if (use.scale) {
    .Deprecated(msg = "'use.scale' is a deprecated argument, please use the 'slot' argument instead")
    slot <- "scale.data"
  }
  if (use.counts) {
    .Deprecated(msg = "'use.counts' is a deprecated argument, please use the 'slot' argument instead")
    if (use.scale) {
      warning("Both 'use.scale' and 'use.counts' were set; using counts",
              call. = FALSE, immediate. = TRUE)
    }
    slot <- "counts"
  }
  fxn.total <- switch(EXPR = slot, data = function(x) {
    rowSums(x = expm1(x = x))
  }, rowSums)
  object.assays <- Seurat:::FilterObjects(object = object, classes.keep = "Assay")
  assays <- assays %||% object.assays
  ident.orig <- Idents(object = object)
  orig.levels <- levels(x = Idents(object = object))
  ident.new <- c()
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(assays) == 0) {
      stop("None of the requested assays are present in the object")
    }
    else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (!is.null(x = add.ident)) {
    new.data <- Seurat:::FetchData(object = object, vars = add.ident)
    new.ident <- paste(Idents(object)[rownames(x = new.data)],
                       new.data[, 1], sep = "_")
    Idents(object, cells = rownames(new.data)) <- new.ident
  }
  data.return <- list()
  for (i in 1:length(x = assays)) {
    data.use <- GetAssayData(object = object, assay = assays[i],
                             slot = slot)
    features.assay <- features
    if (length(x = intersect(x = features, y = rownames(x = data.use))) <
        1) {
      features.assay <- rownames(x = data.use)
    }
    data.all <- list(data.frame(row.names = features.assay))
    for (j in levels(x = Idents(object))) {
      temp.cells <- WhichCells(object = object, idents = j)
      features.assay <- unique(x = intersect(x = features.assay,
                                             y = rownames(x = data.use)))
      if (length(x = temp.cells) == 1) {
        data.temp <- (data.use[features.assay, temp.cells])
        if (slot == "data") {
          data.temp <- expm1(x = data.temp)
        }
      }
      if (length(x = temp.cells) > 1) {
        data.temp <- fxn.total(data.use[features.assay,
                                        temp.cells, drop = FALSE])
      }
      data.all[[j]] <- data.temp
      if (verbose) {
        message(paste("Finished totalling", assays[i],
                      "for cluster", j))
      }
      if (i == 1) {
        ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
      }
    }
    names(x = ident.new) <- levels(x = Idents(object))
    data.return[[i]] <- do.call(cbind, data.all)
    names(x = data.return)[i] <- assays[[i]]
  }
  if (return.seurat) {
    toRet <- CreateSeuratObject(counts = data.return[[1]],
                                project = "Total", assay = names(x = data.return)[1],
                                ...)
    toRet <- SetAssayData(object = toRet, assay = names(x = data.return)[1],
                          slot = "data", new.data = log1p(x = as.matrix(x = data.return[[1]])))
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]])
        toRet <- SetAssayData(object = toRet, assay = names(x = data.return)[i],
                              slot = "data", new.data = log1p(x = as.matrix(x = data.return[[i]])))
      }
    }
    if (DefaultAssay(object = object) %in% names(x = data.return)) {
      DefaultAssay(object = toRet) <- DefaultAssay(object = object)
    }
    Idents(toRet, cells = colnames(x = toRet)) <- ident.new[colnames(x = toRet)]
    Idents(object = toRet) <- factor(x = Idents(object = toRet),
                                     levels = as.character(x = orig.levels), ordered = TRUE)
    toRet <- ScaleData(object = toRet, verbose = verbose)
    return(toRet)
  }
  else {
    return(data.return)
  }
}

#' Differential abundance analysis with miloR
#'
#' Builds a Milo neighbourhood graph on the PCA embedding. Because there are no
#' biological replicates, pseudo-replicates are created and p-values are only
#' indicative. Tests differential abundance by genotype, by condition and by
#' condition within cell type, plots the results as beeswarm plots and saves the
#' Milo object and results as an .rds.
miloAnalysis <- function(){

  library(miloR)
  require(patchwork)

  arid1a.multiome <- readRDS(file=arid1a.object)
  # since we dont have replicates I will create pseudo-replicates
  # we will ignore p-values after this point as they will no longer be valid
  sce <- as.SingleCellExperiment(arid1a.multiome)
  reducedDim(x = sce, "PCA") <- arid1a.multiome@reductions[['pca']]@cell.embeddings
  reducedDim(x = sce, "UMAP") <- arid1a.multiome@reductions[['umap']]@cell.embeddings
  reducedDim(x = sce, "FA") <- arid1a.multiome@reductions[['fa']]@cell.embeddings
  reducedDimNames(sce)

  #### Milo by genotype ####
  df <- colData(sce)
  set.seed(seed = 13)
  samples <- unique(sce$orig.ident)
  reps <- sapply(samples, function(x){
    dfs <- df[df$orig.ident==x,]
    sample(x = rownames(dfs),
           size = ceiling(nrow(dfs)/2), replace = F)
  }) %>% unlist(use.names = F)
  sce$orig.ident <- paste0(sce$orig.ident,'_rep1')
  sce@colData[reps,"orig.ident"] <- gsub("rep1","rep2",sce@colData[reps,"orig.ident"])
  table(sce$orig.ident)
  milo.obj <- Milo(sce)
  miloR::graph(milo.obj) <- miloR::graph(buildFromAdjacency(arid1a.multiome@graphs$SCT_snn, k=30))

  # As per the guidelines in the original Milo paper - the recommendation is to
  # set k such that you have a minimum of Sx5 counts per nhood,
  # where S=total number of samples, so in your case this would be 6*5=30.
  k <- 30

  milo.obj <- makeNhoods(milo.obj,
                         refinement_scheme = "graph",
                         prop = 0.1,
                         k = k,
                         d=30,
                         refined = TRUE,
                         reduced_dims = "PCA")

  plotNhoodSizeHist(milo.obj)

  # Counting cells in neighborhood
  milo.obj <- countCells(milo.obj, meta.data = as.data.frame(colData(milo.obj)), sample="orig.ident")
  head(nhoodCounts(milo.obj))

  # Experimental design
  milo_design <- data.frame(colData(milo.obj))[,c("orig.ident","bioNames")]

  milo_design$bioNames <- gsub("\\(\\+\\/\\+\\)","",milo_design$bioNames)
  milo_design$bioNames <- gsub("\\(\\-\\/\\-\\)","ko",milo_design$bioNames)

  milo_design <- dplyr::distinct(milo_design)
  milo_design$genotype <- gsub(" ","_",milo_design$bioNames)
  rownames(milo_design) <- milo_design$orig.ident
  ## Reorder rownames to match columns of nhoodCounts(milo)
  milo_design <- milo_design[colnames(nhoodCounts(milo.obj)), , drop=FALSE]
  milo_design

  milo.obj <- buildNhoodGraph(milo.obj)

  milo_design

  ## Avg difference of all
  model <- model.matrix(~ 0 + genotype, data=milo_design)
  make.names(colnames(model))
  ave.contrast <- c("(genotypeArid1a_E2 + genotypeArid1ako + genotypeArid1ako_E2)/3 - genotypeArid1a")
  #ave.contrast <- c("(genotypeArid1a + genotypeArid1a_E2 + genotypeArid1ako + genotypeArid1ako_E2)/4")
  mod.contrast <- makeContrasts(contrasts=ave.contrast, levels=model)
  mod.contrast

  da_results_avg <- testNhoods(milo.obj,
                               design = ~0 + genotype,
                               design.df = milo_design,
                               model.contrasts = ave.contrast,
                               fdr.weighting = "graph-overlap")

  message("!Ignore all p-values as we do not have true replicates!")

  library(scales)
  # base FA
  g0 <- DimPlot(arid1a.multiome, label = F, group.by = 'bioNames', reduction = "fa") %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    pt.size = 0.3,
                    panel.border.size = 0,
                    pal=bionames_palette,
                    title = "Prostate organoid genotypes")

  ggplot(da_results_avg, aes(logFC, -log10(SpatialFDR))) +
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

  da_results_avg <- annotateNhoods(milo.obj, da_results_avg, coldata_col = "bioNames")
  head(da_results_avg)

  da_results_avg$bioNames <- factor(x = da_results_avg$bioNames,
                                    levels = c('Arid1a(+/+)','Arid1a(+/+) E2','Arid1a(-/-)','Arid1a(-/-) E2'),
                                    ordered = T )

  # beeswarm
  gb <- plotDAbeeswarm(da_results_avg, group.by = "bioNames", alpha = 1)
  gb <- gb + scale_color_gradient2(low=muted("blue"), mid="white", high=muted("red"),midpoint=0) + xlab("") +
    scale_y_continuous(limits = c(-10,10)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'gray90', linewidth = .2))

  pdf(file = file.path(plot.dir,'Beeswarm_genotype_avg_vs_Arid1a_milo.pdf'), width = 5, height = 3, bg = 'transparent')
  gb
  dev.off()

  # FA
  g1 <- plotNhoodGraphDA(milo.obj, da_results_avg, alpha=1, layout="FA") +
    plot_layout(guides="collect") + ggtitle("Milo avg(Genotype)") +
    scale_fill_gradient2(low=muted("blue"), mid="white", high=muted("red"), midpoint=0) + # reverse colors
    ylab("") +
    coord_fixed()

  pdf(file = file.path(plot.dir,'UMAP_genotype_avg_vs_Arid1a_milo.pdf'), width = 16, height = 8, bg = 'transparent')
  g0 + g1 + gb
  dev.off()

  ## cell type
  da_results <- annotateNhoods(milo.obj, da_results_avg, coldata_col = "mammary_gland_refine_clusters")
  head(da_results)
  gb <- plotDAbeeswarm(da_results, group.by = "mammary_gland_refine_clusters", alpha = 1)
  gb <- gb + scale_color_gradient2(low=muted("blue"), mid="white", high=muted("red"),midpoint=0) + xlab("") +
    scale_y_continuous(limits = c(-10,10)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'gray90', linewidth = .2))

  pdf(file = file.path(plot.dir,'Beeswarm_celltype_avg_vs_Arid1a_milo.pdf'), width = 5, height = 3, bg = 'transparent')
  gb
  dev.off()

  #### Milo by condition ####
  df <- colData(sce)
  set.seed(seed = 13)
  samples <- unique(sce$orig.ident)
  reps <- sapply(samples, function(x){
    dfs <- df[df$orig.ident==x,]
    sample(x = rownames(dfs),
           size = ceiling(nrow(dfs)/2), replace = F)
  }) %>% unlist(use.names = F)
  sce$orig.ident <- paste0(sce$orig.ident,'_rep1')
  sce@colData[reps,"orig.ident"] <- gsub("rep1","rep2",sce@colData[reps,"orig.ident"])
  table(sce$orig.ident)
  milo.obj <- Milo(sce)
  miloR::graph(milo.obj) <- miloR::graph(buildFromAdjacency(arid1a.multiome@graphs$SCT_snn, k=30))

  # As per the guidelines in the original Milo paper - the recommendation is to
  # set k such that you have a minimum of Sx5 counts per nhood,
  # where S=total number of samples, so in your case this would be 6*5=30.
  k <- 30
  milo.obj@.k=k
  milo.obj <- makeNhoods(milo.obj,
                         refinement_scheme = "graph",
                         prop = 0.1,
                         k = k,
                         d=30,
                         refined = TRUE,
                         reduced_dims = "PCA")

  plotNhoodSizeHist(milo.obj)

  #**
  # milo.obj <- calcNhoodDistance(milo.obj, d=30)
  # milo.obj <- countCells(milo.obj, samples="Sample", meta.data=milo.meta)
  #
  # milo.design <- as.data.frame(xtabs(~ Condition + Sample, data=milo.meta))
  # milo.design <- milo.design[milo.design$Freq > 0, ]
  # rownames(milo.design) <- milo.design$Sample
  # milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]
  #
  # milo.res <- testNhoods(milo.obj, design=~Condition, design.df=milo.design)
  # head(milo.res)
  #*

  milo.meta <- arid1a.multiome@meta.data
  milo.obj <- countCells(milo.obj, samples="orig.ident", meta.data=milo.meta)

  milo.design <- as.data.frame(xtabs(~ condition + orig.ident, data=milo.meta))
  milo.design <- milo.design[milo.design$Freq > 0, ]
  rownames(milo.design) <- milo.design$orig.ident
  milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]

  milo.res <- testNhoods(milo.obj, design=~condition, design.df=milo.design)
  head(milo.res)

  da_results <- annotateNhoods(milo.obj, milo.res, coldata_col = "bioNames")
  head(da_results)

  da_results$bioNames <- factor(x = da_results$bioNames,
                                    levels = c('Arid1a(+/+)','Arid1a(+/+) E2','Arid1a(-/-)','Arid1a(-/-) E2'),
                                    ordered = T )

  # beeswarm
  gb <- plotDAbeeswarm(da_results, group.by = "bioNames", alpha = 1)
  gb <- gb + scale_color_gradient2(low=muted("blue"), mid="white", high=muted("red"),midpoint=0) +
    xlab("") + scale_y_continuous(limits = c(-10,10)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'gray90', linewidth = .2))
  gb

  #### Milo by Condition x Cell Type ####
  library(miloR)
  require(patchwork)
  require(dplyr)

  # since we dont have replicates I will create pseudo-replicates
  # we will ignore p-values after this point as they will no longer be valid
  sce <- as.SingleCellExperiment(arid1a.multiome)
  reducedDim(x = sce, "PCA")  <- arid1a.multiome@reductions[['pca']]@cell.embeddings
  reducedDim(x = sce, "UMAP") <- arid1a.multiome@reductions[['umap']]@cell.embeddings
  reducedDim(x = sce, "FA") <- arid1a.multiome@reductions[['fa']]@cell.embeddings
  reducedDimNames(sce)
  df <- colData(sce)
  set.seed(seed = 13)
  (samples <- unique(sce$orig.ident))
  reps <- sapply(samples, function(x){
    dfs <- df[df$orig.ident==x,]
    sample(x = rownames(dfs),
           size = ceiling(nrow(dfs)/2), replace = F)
  }) %>% unlist(use.names = F)
  sce$orig.ident <- paste0(sce$orig.ident,'_rep1')
  sce@colData[reps,"orig.ident"] <- gsub("rep1","rep2",sce@colData[reps,"orig.ident"])
  table(sce$orig.ident)
  milo.obj <- Milo(sce)

  ## Version 1.
  miloR::graph(milo.obj) <- miloR::graph(buildFromAdjacency(arid1a.multiome@graphs$SCT_snn, k=30))

  # As per the guidelines in the original Milo paper - the recommendation is to
  # set k such that you have a minimum of Sx5 counts per nhood,
  # where S=total number of samples, so in your case this would be 6*5=30.
  k <- 30

  milo.obj <- makeNhoods(milo.obj,
                         refinement_scheme = "graph",
                         prop = 0.1,
                         k = k,
                         d=30,
                         refined = TRUE,
                         reduced_dims = "PCA")

  plotNhoodSizeHist(milo.obj)

  ## Version 2.
  milo.obj <- buildGraph(milo.obj, k = 30, d = 30, reduced.dim = "PCA")
  milo.obj <- makeNhoods(milo.obj, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")
  plotNhoodSizeHist(milo.obj)

  # Counting cells in neighborhood
  milo.obj <- countCells(milo.obj,
                         meta.data = as.data.frame(colData(milo.obj)),
                         sample="orig.ident")
  head(nhoodCounts(milo.obj))

  # Experimental design
  milo_design <- data.frame(colData(milo.obj))[,c("orig.ident","bioNames","condition")]
  head(milo_design)
  milo_design$bioNames <- gsub("\\(\\+\\/\\+\\)","_wt",milo_design$bioNames)
  milo_design$bioNames <- gsub("\\(\\-\\/\\-\\)","_ko",milo_design$bioNames)

  head(milo_design)

  milo_design <- dplyr::distinct(milo_design)
  rownames(milo_design) <- milo_design$orig.ident
  ## Reorder rownames to match columns of nhoodCounts(milo)
  milo_design <- milo_design[colnames(nhoodCounts(milo.obj)), ,
                             drop=FALSE]
  milo_design

  ## Neighborhood connectivity ##
  milo.obj <- calcNhoodDistance(milo.obj, d=30, reduced.dim = "PCA")

  da_results <- testNhoods(milo.obj,
                           design = ~condition,
                           design.df = milo_design,
                           fdr.weighting = "graph-overlap")
  message("!Ignore all p-values as we do not have true replicates!")

  head(da_results)

  da_results %>%
    arrange(SpatialFDR) %>%
    head()

  ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

  ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

  library(scales)
  milo.obj <- buildNhoodGraph(milo.obj)
  # UMAP
  umap_pl <- DimPlot(arid1a.multiome, label = F, group.by = 'bioNames', reduction = "umap") %>%
    UmapPlotEnhance(obj = arid1a.multiome, p1 = .,
                    pt.size = 0.3,
                    panel.border.size = 0,
                    pal=bionames_palette,
                    title = "Mammary genotypes")
  umap_pl
  ## Plot neighbourhood graph
  nh_graph_pl <- plotNhoodGraphDA(milo.obj, da_results, alpha=0.1) +
    guides(fill = guide_colorbar(direction = "horizontal"))

  umap_pl + nh_graph_pl + plot_layout(guides="collect")

  pdf(file = file.path(plot.dir,"milo_umap_condition.pdf"), width = 6, height = 5)
  umap_pl / nh_graph_pl + plot_layout(guides="collect")
  dev.off()

  #### beeswarm mammary_gland_refine_clusters ####

  da_results <- annotateNhoods(milo.obj, da_results, coldata_col = c("mammary_gland_refine_clusters"))
  head(da_results)

  gb <- plotDAbeeswarm(da_results, group.by = "mammary_gland_refine_clusters", alpha = 1)
  gb <- gb + scale_color_gradient2(low=muted("red"), mid="white", high=muted("blue"),midpoint=0) +
    xlab("") +
    scale_y_continuous(limits = c(-10,10)) +
    ggtitle(label = "Arid1a WT vs Arid1a KO", subtitle = "KO ↔ WT") +
    theme(title = element_text(hjust = .5), plot.subtitle = element_text(hjust = .5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'gray90', linewidth = .2))

  cairo_pdf(filename = file.path(plot.dir,'Beeswarm_condition_vs_mammary_gland_refine_clusters_Arid1a_milo.pdf'),
            width = 7, height = 6, bg = 'transparent')
  gb
  dev.off()

  saveRDS(object = list('milo.obj'=milo.obj, 'da_results'=da_results),
          file = file.path(robj.dir,"MiloResults.rds"))

  #### beeswarm bionames ####

  da_results <- annotateNhoods(milo.obj, da_results, coldata_col = c("bioNames"))
  head(da_results)
  da_results$bioNames <- factor(x = da_results$bioNames,
                                levels = c("Arid1a(-/-)","Arid1a(-/-) E2",
                                           "Arid1a(+/+)","Arid1a(+/+) E2"))

  gb <- plotDAbeeswarm(da_results, group.by = "bioNames", alpha = 1)
  gb <- gb + scale_color_gradient2(low=muted("red"), mid="white", high=muted("blue"),midpoint=0) +
    xlab("") +
    scale_y_continuous(limits = c(-10,10)) +
    ggtitle(label = "Arid1a WT vs Arid1a KO", subtitle = "KO ↔ WT") +
    theme(title = element_text(hjust = .5), plot.subtitle = element_text(hjust = .5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'gray90', linewidth = .2))

  cairo_pdf(filename = file.path(plot.dir,'Beeswarm_condition_vs_bioNames_Arid1a_milo.pdf'),
            width = 7, height = 5, bg = 'transparent')
  gb
  dev.off()

  #### beeswarm clusters ####

  da_results <- annotateNhoods(milo.obj, da_results, coldata_col = c("seurat_clusters"))
  head(da_results)
  if(is.factor(da_results$seurat_clusters))
    da_results$seurat_clusters <- as.numeric(levels(da_results$seurat_clusters))[da_results$seurat_clusters]
  else
    da_results$seurat_clusters <- as.integer(da_results$seurat_clusters)
  da_results$seurat_clusters <- factor(x = da_results$seurat_clusters,
                                       levels = rev(sort(unique(da_results$seurat_clusters))))
  head(da_results)

  gb <- plotDAbeeswarm(da_results, group.by = "seurat_clusters", alpha = 1)
  gb <- gb + scale_color_gradient2(low=muted("red"), mid="white", high=muted("blue"),midpoint=0) +
    xlab("") +
    scale_y_continuous(limits = c(-10,10)) +
    ggtitle(label = "Arid1a WT vs Arid1a KO", subtitle = "KO ↔ WT") +
    theme(title = element_text(hjust = .5), plot.subtitle = element_text(hjust = .5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'gray90', linewidth = .2))

  cairo_pdf(filename = file.path(plot.dir,'Beeswarm_condition_vs_seurat_clusters_Arid1a_milo.pdf'),
            width = 7, height = 6.5, bg = 'transparent')
  gb
  dev.off()

}

#' Functional annotation of Hotspot gene modules
#'
#' Reads the gene modules computed with Hotspot in scanpy
#' (`<data.dir>/scanpy/hotspot_gene_modules.txt`), writes them to a spreadsheet
#' and runs gProfiler enrichment per module (GO, Hallmark and GSEA G5 sets via
#' custom gProfiler organisms), writing the results to `spreadsheets.dir`.
hotspot_go <- function(){

  dt <- fread(file.path(data.dir, 'scanpy', 'hotspot_gene_modules.txt'))

  writexl::write_xlsx(x = dt, path = file.path(spreadsheets.dir,'hotspot_gene_modules.xlsx'))

  library(gprofiler2)
  library(ggrepel)
  library(magrittr)

  totExp <- TotalExpression(object = arid1a.multiome)
  mu <- rowMeans(totExp$RNA)
  quantile(mu, probs=c(seq(0,1,by=.1)))
  (minexp <- quantile(mu, probs=.25))
  expgenes <- mu[mu >= minexp]

  (modules <- unique(dt$Module))
  gostres <- lapply(modules, function(i){
    query <- list(dt[Module==i][FDR < 0.05, Gene])
    names(query) <- paste0("Module_",i)
    message("Gathering results for ",names(query))
    gost(query = query,
         sources = c("GO", "KEGG", "REAC"),
         organism = "mmusculus",
         ordered_query = TRUE, # Items are ordered by Hotspot
         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
         measure_underrepresentation = FALSE, evcodes = FALSE,
         user_threshold = 0.05, correction_method = "g_SCS",
         #domain_scope = "annotated", custom_bg = NULL,
         #domain_scope = "known", custom_bg = NULL,
         #domain_scope = "custom_annotated", custom_bg = foxa1.multiome@assays$RNA@var.features,
         #domain_scope = "custom",
         #custom_bg = names(expgenes),
         custom_bg = rownames(arid1a.multiome),
         numeric_ns = "", as_short_link = FALSE)
  }) %>% set_names(paste0("m",modules))

  dtgostres <- lapply(gostres, function(x){ data.table(x$result) }) %>% rbindlist()

  writexl::write_xlsx(x = dtgostres, path = file.path(spreadsheets.dir,'go_results_of_hotspot_gene_modules.xlsx'))

  ## If grid port error then expand the grid for publish_gostplot
  plot.width <- 15
  plot.height <- 15
  gp <- lapply(names(gostres), function(i){
    if(!is.null(gostres[[i]])){
      p <- gostplot(gostres[[i]], capped = T, interactive = F)
      dt <- data.table(p$data)
      dt <- dt[order(p_value, decreasing = F),head(.SD, 10), by=.group] # label top n by group
      maxy <- ceiling(max(dt$logpval))
      p <- p + geom_text_repel(data = dt, aes(x=order,
                                              y=logpval,
                                              label=term_name),
                               size = 2,
                               min.segment.length = 0, # draw all line segments
                               max.overlaps = 15) + scale_y_continuous(limits = c(NA, maxy))
      terms <- dt[,head(.SD, n=6*length(unique(source))), by=list(source)]
      p <- publish_gostplot(p = p, highlight_terms = terms, width = plot.width, height = plot.height)
    }
  }) %>% set_names(names(gostres))

  # remove empty entries
  gp <- gp[!sapply(gp,is.null)]

  pdf(file = file.path(plot.dir,"hotspot_go_expressed.pdf"), width = plot.width, height = plot.height, bg = "transparent")
  gp
  dev.off()

  ## Table of results
  gt <- lapply(names(gostres), function(i){
    p <- publish_gosttable(gostres[[i]],
                           highlight_terms = head(gostres[[i]]$result, 80),
                           use_colors = TRUE,
                           show_columns = c("source", "term_name", "term_size", "intersection_size"))
    p <- p + ggtitle(i) + theme(title = element_text(face = "bold", size = 16, hjust = 0.5))
  }) %>% set_names(names(gostres))

  pdf(file = file.path(plot.dir,"hotspot_go_annotated_tables.pdf"), width = 10, height = 22, bg = "transparent")
  gt
  dev.off()

  #### custom organism ####
  # Hallmark
  #upload_GMT_file(gmtfile = file.path(annot.dir, 'mh.all.v2023.2.Mm.symbols.gmt'))
  # gp__VECM_xsGW_iuI
  #
  # GO terms
  # upload_GMT_file(gmtfile = file.path(annot.dir, 'm5.all.v2023.2.Mm.symbols.gmt'))
  # gp__XGL1_P9f8_rl4

  library(gprofiler2)
  library(ggrepel)
  library(magrittr)

  dt <- fread(file.path(data.dir, 'scanpy', 'hotspot_gene_modules.txt'))
  custom_gostres_organism <- Sys.getenv("GPROFILER_HALLMARK_ID")
  stopifnot(nzchar(custom_gostres_organism))

  (modules <- unique(dt$Module))
  gostres <- lapply(modules, function(i){
    query <- list(dt[Module==i][FDR < 0.05, Gene])
    names(query) <- paste0("Module_",i)
    message("Gathering results for ",names(query))
    gost(query = query,
         organism = custom_gostres_organism,
         ordered_query = TRUE, # Items are ordered by Hotspot
         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
         measure_underrepresentation = FALSE, evcodes = FALSE,
         user_threshold = 0.05, correction_method = "g_SCS",
         domain_scope = "annotated", custom_bg = NULL,
         #domain_scope = "known", custom_bg = NULL,
         #domain_scope = "custom_annotated", custom_bg = foxa1.multiome@assays$RNA@var.features,
         #domain_scope = "custom_annotated", custom_bg = rownames(foxa1.multiome),
         numeric_ns = "", as_short_link = FALSE)
  }) %>% set_names(paste0("m",modules))

  dtgostres <- lapply(gostres, function(x){ data.table(x$result) }) %>% rbindlist()

  writexl::write_xlsx(x = dtgostres, path = file.path(spreadsheets.dir,'hallmark_results_of_hotspot_gene_modules.xlsx'))

  gp <- lapply(names(gostres), function(i){
    if(!is.null(gostres[[i]])){
      p <- gostplot(gostres[[i]], capped = T, interactive = F)
      pd <- p$data
      pd$term_name <- gsub("https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/","",pd$term_name)
      p$data <- pd
      dt <- data.table(p$data)
      dt <- dt[,head(.SD, 10), by=.group] # label top n by group

      p <- p + geom_text_repel(data = dt, aes(x=order,
                                              y=logpval,
                                              label=term_name),
                               size = 2,
                               min.segment.length = 0, # draw all line segments
                               max.overlaps = 15)
      p <- publish_gostplot(p = p, highlight_terms = head(dt, 30))
    }
  }) %>% set_names(names(gostres))

  pdf(file = file.path(plot.dir,"hotspot_hallmark_annotated.pdf"), width = 12, height = 8, bg = "transparent")
  gp
  dev.off()

  # GO terms from GSEA
  dt <- fread(file.path(data.dir, 'scanpy', 'hotspot_gene_modules.txt'))
  custom_gostres_organism <- Sys.getenv("GPROFILER_G5_ID")
  stopifnot(nzchar(custom_gostres_organism))

  (modules <- unique(dt$Module))
  gostres <- lapply(modules, function(i){
    query <- list(dt[Module==i][FDR < 0.05, Gene])
    names(query) <- paste0("Module_",i)
    message("Gathering results for ",names(query))
    gost(query = query,
         organism = custom_gostres_organism,
         ordered_query = TRUE, # Items are ordered by Hotspot
         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
         measure_underrepresentation = FALSE, evcodes = FALSE,
         user_threshold = 0.05, correction_method = "g_SCS",
         domain_scope = "annotated", custom_bg = NULL,
         #domain_scope = "known", custom_bg = NULL,
         #domain_scope = "custom_annotated", custom_bg = foxa1.multiome@assays$RNA@var.features,
         #domain_scope = "custom_annotated", custom_bg = rownames(foxa1.multiome),
         numeric_ns = "", as_short_link = FALSE)
  }) %>% set_names(paste0("m",modules))

  dtgostres <- lapply(gostres, function(x){ data.table(x$result) }) %>% rbindlist()

  writexl::write_xlsx(x = dtgostres, path = file.path(spreadsheets.dir,'GSEA_G5_results_of_hotspot_gene_modules.xlsx'))

  gp <- lapply(names(gostres), function(i){
    if(!is.null(gostres[[i]])){
      p <- gostplot(gostres[[i]], capped = T, interactive = F)
      pd <- p$data
      pd$term_name <- gsub("https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/","",pd$term_name)
      p$data <- pd
      dt <- data.table(p$data)
      dt <- dt[,head(.SD, 10), by=.group] # label top n by group

      p <- p + geom_text_repel(data = dt, aes(x=order,
                                              y=logpval,
                                              label=term_name),
                               size = 2,
                               min.segment.length = 0, # draw all line segments
                               max.overlaps = 15)
      p <- publish_gostplot(p = p, highlight_terms = head(dt, 20))
    }
  }) %>% set_names(names(gostres))
  # keep non null entries
  gp <- gp[!sapply(gp,is.null)]

  pdf(file = file.path(plot.dir,"hotspot_GSEA_G5_annotated.pdf"), width = 8, height = 16, bg = "transparent")
  gp
  dev.off()

}

#' Silhouette widths of a clustering in a reduced space
#'
#' @param obj Seurat object.
#' @param reduction Name of the reduction whose embedding is used for distances.
#' @param clusters Name of the meta.data column holding the cluster labels.
#' @param dims Embedding dimensions to use.
#' @return A `silhouette` object (cluster package).
silhouetteScore <- function(obj, reduction='umap', clusters='seurat_clusters', dims=c(1:30)){
  library(cluster)

  dist.matrix <- dist(x = Embeddings(object = obj[[reduction]])[, dims])
  clusters <- obj@meta.data[[clusters]]
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
}

#### Plots ####
#' Publication-style UMAP plot
#'
#' Wraps a Seurat DimPlot/FeaturePlot (`p1`) or builds one from `group.by`, and
#' redraws it with custom palettes, rasterised points, tidy legends (optionally
#' with per-group cell counts), a fixed panel size and optional faceting. Writes a
#' png (or pdf if `as.pdf = TRUE`) named `filename.prefix` in `plot.dir`.
#' Adapted from Seurat's DimPlot/SingleDimPlot.
#'
#' @param obj Seurat object.
#' @param p1 Optional Seurat plot to take the data from.
#' @param group.by Meta.data column used for colouring when `p1` is NULL.
#' @param pal Colour palette (named vector for discrete groups).
#' @param filename.prefix File name (without extension) of the saved figure.
#' @param as.pdf Write a pdf instead of a png.
#' @return The ggplot object.
UmapPlotEnhance <- function(obj,
                            group.by=NULL,
                            p1=NULL,
                            cells=NULL,
                            addLabels=NULL,
                            axis.label=NA, #"orig.ident",
                            outline.size=0,
                            outline.pt.size=0,
                            outline.color='black', #'grey80',
                            pal=NA,
                            panel.border.size=0,
                            show.legend=T,
                            label.size=4,
                            title.size=14,
                            scale.color=F,
                            pt.size=0.1,
                            pt.alpha=1,
                            title=NA,
                            legend.title="",
                            legend.title.rotation=90,
                            legend.title.size=12,
                            legend.height.factor=.4,
                            legend.position="right",
                            legend.justification="top",
                            legend.direction="vertical",
                            legend.text.size=12,
                            legend.show.na=T,
                            legend.ncol=NA,
                            legend.size=5,
                            break.unit=NULL,
                            facet.text.size=10,
                            fixed_coord=T,
                            panel_width_in=3,
                            add.legend.counts=T,
                            order.random=F, # ordered for higher expression on top
                            vmin=NA,
                            vmax=NA,
                            na.color="gray80",
                            filename.prefix,
                            nrow=NULL,
                            reduction=NULL,
                            coord1='UMAP_1',
                            coord2='UMAP_2',
                            small_axis_title=T,
                            resolution=150,
                            as.pdf=F){
  # Plot the integration UMAP
  require(shadowtext)
  require(ggrastr)
  require(ggh4x)
  require(gtable)

  dt <- NULL

  if(is.null(p1)){
    dims = c(1, 2)
    reduction <- reduction %||% DefaultDimReduc(object = obj)
    cells <- cells %||% Cells(x = obj, assay = DefaultAssay(object = obj[[reduction]]))
    dims <- paste0(Key(object = obj[[reduction]]), dims)
    embed <- Embeddings(obj[[reduction]])[cells, dims, drop = FALSE]
    embed <- as.data.frame(embed)
    orig.groups <- group.by
    group.by <- group.by %||% "ident"
    meta <- FetchData(object = obj,
                      layer = "data",
                      vars = group.by,
                      cells = cells,
                      clean = "none")
    # 44346
    message("Dim embed:",dim(embed))
    message("Dim meta:",dim(meta))
    data <- cbind(embed, meta)
    group.by <- colnames(x = data)[3:ncol(x = data)]
    for (group in group.by) {
      # Want this for discrete colors
      if (!is.factor(x = data[, group]) & group %in% colnames(obj@meta.data) ){
        data[, group] <- factor(x = data[, group])
      }
    }
    dt <- data.table(data, keep.rownames = 'rn')
  }
  else if(is.data.frame(p1))
    dt <- data.table(p1, keep.rownames = 'rn')
  else
    dt <- data.table(p1$data, keep.rownames = 'rn')
  cols <- colnames(dt)

  panel_width <- unit(panel_width_in, "in")
  if(length(title)==0){
    title <- data.table::last(cols)
  }
  # handle alt reductions
  if(!coord1 %in% colnames(dt)){
    coord1 <- colnames(dt)[2]
    coord2 <- colnames(dt)[3]
  }

  setnames(dt, old = c(coord1,coord2), new = c("x","y"))
  #dt[,color:=data.table::last(colnames(dt))]
  setnames(dt, data.table::last(colnames(dt)), new = c("color"))

  if(add.legend.counts){
    message("axis label")
    dt[,count:=.N,by=color]
  }

  ## This creates 0 valued points not in faceting window
  if(axis.label %in% colnames(obj@meta.data)){
    dt[,axis.label := obj@meta.data[dt$rn,eval(axis.label)]]
    dt <- lapply(unique(dt$axis.label), function(i){
      ds <- copy(dt[axis.label!=i])
      dt_zero <- copy(ds)
      dt_zero$color <- NA
      dt_zero$axis.label <- i
      rbindlist(list(dt_zero, ds))
    }) %>% rbindlist
  }
  # build the data taable
  legend.continuous <- ifelse(is.factor(dt$color), F, T)
  if(legend.continuous){
    message("Continuous data [",paste(range(dt$color), collapse = ":"),"]")
  } else {
    message("Discrete data range [",paste(levels(dt$color), collapse = ":"),"]")
  }
  if(order.random)
    dt <- dt[order(sample(color), na.last = F)]
  else{
    message("Order by color (ascending)")
    dt <- dt[order(color, decreasing = F, na.last = F)]
  }
  dt[,rn := factor(rn, levels = unique(rn), ordered = T)]

  # build the plot
  p1 <- NULL
  if(all(is.na(pal))){
    p1 <- ggplot(dt, aes(x=x, y=y)) #, fill=color
    outline.size <- 0
  }
  else{
    p1 <- ggplot(dt, aes(x=x, y=y, fill=color))
  }

  if(is.null(outline.size))
    outline.size = pt.size

  if(outline.size > 0){
    if(is.null(outline.pt.size))
      outline.pt.size <- pt.size
    p1 <- p1 + geom_point_rast(size=outline.pt.size, shape=1,
                               stroke=outline.size, col=outline.color) # standard stroke #'grey40'
  }

  p1 <- p1 + geom_point_rast(size=pt.size, shape=21, stroke=0, alpha=pt.alpha) +
    theme_bw()

  p1 <- p1 + xlab( gsub("_","",coord1) ) + ylab( gsub("_","",coord2) ) +
    theme(text = element_text(size = 15),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=title.size, hjust = 0.5),
          legend.text = element_text(size = legend.text.size),
          strip.text = element_text(size = facet.text.size),
          panel.border = element_rect(linewidth = panel.border.size, color='black'),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          plot.margin = margin(t=0,r=0,b=0,l=0),
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.direction = legend.direction,
          legend.position = legend.position,
          legend.justification = legend.justification,
          legend.title = element_text(angle = legend.title.rotation, size = legend.title.size),
          legend.background = element_blank(), # get rid of legend bg element_rect(fill = 'transparent')
          legend.box.background = element_blank(), # get rid o
          legend.key=element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(t=0,r=0,b=0,l=0)) #trbl
  if(fixed_coord)
    p1 <- p1 + coord_fixed(clip = 'off')

  if(!is.null(addLabels)){
    if(addLabels %in% colnames(obj@meta.data)){
      ds <- setDT(DimPlot(obj, label = F, group.by = addLabels)$data)
      setnames(ds, c("x","y","label"))
      #names(pdfFonts())
      median.labels <- ds[,.(x=median(x),y=median(y), color=I("black")), by=label]
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

  #  get dimensions of plot
  gt <- ggplotGrob(p1)
  panel <- gtable_filter(gt, "panel")
  # height is scaled compared to width
  panel_height_in <- round(panel_width_in*as.numeric(gsub("null","",panel$heights)), 2)
  # set the panel height based on width
  panel_height <- unit(panel_height_in, "in") # reduce for title
  p1 <- p1 + force_panelsizes(cols = panel_width, rows = panel_height)
  p1 <- p1 + ggtitle(label = title) + theme(plot.title = element_text(size = title.size))

  # Do we want color fill?
  if(all(is.na(pal))){
    ##nothing
  }
  else {
    if(legend.continuous){
      message("continous color")
      if(missing(pal))
        pal <- pals::brewer.purples(10)

      lim <- range(dt$color, na.rm = T)
      lim <- c(floor(1000*lim[1])/1000, ceiling(1000*lim[2])/1000)

      if(!is.na(vmin))
        #lim[1] <- min(dt[color>=vmin]$color, na.rm = T)
        lim[1] <- min(vmin, dt$color, na.rm = T)

      if(!is.na(vmax))
        lim[2] <- max(vmax, dt$color, na.rm = T)

      if(scale.color)
        lim <- c(-max(abs(dt$color), na.rm = T), max(abs(dt$color), na.rm = T))
      else if(lim[1] > 0)
        lim <- c(0,max(abs(dt$color), na.rm = T))

      if(is.null(break.unit)){
        if(ceiling(lim[2]) > 1)
          break.unit <- 1
        else
          break.unit <- .1
        if( identical( range(c(0,1)), lim ) )
          break.unit <- 0.2
      }

      if(is.na(legend.ncol))
        legend.ncol <- 1

      p1 <- p1 + scale_fill_gradientn(legend.title,
                                      guide = "colourbar",
                                      colours = pal,
                                      na.value = na.color,
                                      limits=lim,
                                      breaks=seq(from=floor(lim[1]), to=ceiling(lim[2]), by=break.unit))

      legend.title.position <- ifelse(legend.title=="","right","left")
      p1 <- p1 + guides(fill = guide_colourbar(title.position = legend.title.position,
                                               ticks = T,
                                               ticks.linewidth = .5,
                                               ticks.colour = 'black',
                                               frame.colour = 'black',
                                               frame.linewidth = 0.2,
                                               nbin = 1000, # high to set extreme ends
                                               barwidth = .5,
                                               barheight = (panel_height / ifelse(is.null(nrow), 1, nrow))*legend.height.factor))
    }
    #discrete
    else{
      legend.title.position <- ifelse(legend.title.rotation=="0","top","left")
      n <- length(unique(na.omit( dt$color ) ))
      message("discrete color : ", n)
      if(missing(pal))
        pal <- pals::brewer.accent(n = n)
      if(is.na(legend.ncol))
        legend.ncol <- ifelse(n>10,ifelse(n>20,3,2),1)
      if(is.null(names(pal)))
        names(pal) <- unique(dt$color)
      if(add.legend.counts & is.na(axis.label)){
        cnts <- unique(dt[, count,by=color])
        names(pal) <- sapply(names(pal), function(x){ cnts[color==x, paste0(x," (n=",count,")")] } )
        dt[,color :=                                      paste0(color," (n=",.N,")"), by=color]
      }
      p1 <- p1 + scale_fill_manual(legend.title, values = pal, na.value = na.color)
      p1 <- p1 + guides(fill = guide_legend(title.position=legend.title.position,
                                            override.aes = list(size=legend.size),
                                            ncol=legend.ncol))
      if(legend.height.factor == 1)
        p1 <- p1 + theme(legend.box.margin=margin(-10,10,-10,-10), legend.key.height = unit(1, "null")) #trbl
      # accomodate labels
      #View(dt)
    }
  }

  if(!show.legend){
    p1 <- p1 + theme(legend.position = 'none')
    panel_width_in <- panel_width_in+.5
    panel_height_in <- panel_height_in+.5
  }
  else{
    gt <- ggplotGrob(p1)
    known_wd = sum(grid::convertWidth(gt$widths, 'in', valueOnly=TRUE))
    #message("known wd ",known_wd)
    known_ht = sum(grid::convertHeight(gt$heights, 'in', valueOnly=TRUE))
    #message("known ht ",known_ht)
    panel_width_in <- round(known_wd + .1, 2)
    panel_height_in <- round(known_ht + .1, 2)
    #message("panel_width_in ",panel_width_in)
  }

  if(legend.title!=""){
    p1 <- p1 + theme(legend.box.margin=margin(-10,10,-10,-10),
                     legend.title = element_text(angle = legend.title.rotation)) #trbl
  }

  if("axis.label" %in% colnames(dt)){
    if(!is.null(nrow))
      p1 <- p1 + facet_wrap(~ axis.label, nrow = nrow)
    else
      p1 <- p1 + facet_wrap(~ axis.label)
  }

  if(small_axis_title){
    gm <- max(panel_width_in, panel_height_in)
    xmin <- floor(min(dt[,x]))
    xmax <- xmin + gm*.3
    ymin <- floor(min(dt[,y]))
    ymax <- ymin + gm*.3

    p1 <- p1 +
      # vertical arrow
      annotate("segment", x=xmin, y=ymin, xend=xmin, yend=ymax, linewidth=.5,
               arrow = arrow(type = "closed", length = unit(0.15, "cm")), color="black") +
      annotate("segment", x=xmin, y=ymin, xend=xmax, yend=ymin, linewidth=.5,
               arrow = arrow(type = "closed", length = unit(0.15, "cm")), color="black") +
      theme(axis.title.x =  element_text(hjust = 0, vjust = 2, size = 12),
            axis.title.y = element_text(hjust = 0, vjust = 0, size = 12))

  }

  if(!missing(filename.prefix)){
    if(as.pdf){
      f <- file.path(plot.dir,paste0(filename.prefix,".pdf"))
      message("Plotting [",panel_width_in,"x",panel_height_in,"] to:",f)
      pdf(file = f, width = panel_width_in, height = panel_height_in, bg = 'transparent')
      print(p1)
      dev.off()
    } else{
      f <- file.path(plot.dir,paste0(filename.prefix,".png"))
      message("Plotting [",panel_width_in,"x",panel_height_in,"] to:",f)
      png(file = f,
          width = panel_width_in, height = panel_height_in, units = 'in', bg = 'transparent', res = resolution)
      print(p1)
      dev.off()
    }
  }
  p1
}

#' UMAP coloured by local point density, split by `feature`
#'
#' @param obj Seurat object.
#' @param feature Meta.data column to facet by.
#' @param pal Colour palette for the density.
#' @param filename.prefix File name (without extension) of the saved figure.
#' @return The ggplot object.
UmapPlotDensity <- function(obj,
                            feature='Phase',
                            title.size=16,
                            pt.size=0.8,
                            panel_width_in=3.5,
                            pal=rev(pals::kovesi.linear_kryw_5_100_c64(10)), #pal=pals::viridis(10), #pal=pals::parula(10)
                            na.color="gray80",
                            contour.color = NA,
                            filename.prefix,
                            fixed_coord=T,
                            embedding='umap',
                            coord1='UMAP_1',
                            coord2='UMAP_2',
                            legend.height.factor=0.5,
                            small_axis_title=T,
                            nrow=1,
                            as.pdf=F){

  .get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }

  dt <- merge( data.table( obj@reductions[[embedding]]@cell.embeddings, keep.rownames = 'cellid' ),
               data.table( obj@meta.data, keep.rownames = 'cellid' ), by='cellid')
  dt[,temp_features := get(feature)]

  dt <- dt[!is.na(get(coord1))]

  require(ggh4x)
  require(gtable)

  panel_height_in <<- 0

  p <- lapply(sort(unique(dt[,temp_features])), function(s){
    # select data of interest
    toplot <- dt[temp_features == s]
    if(nrow(toplot) > 1){
      toplot[,density := .get_density(x = get(coord1), y= get(coord2), n=500)]
      toplot[,density := density/max(density)]
      palx <- pal
      if(is.list(pal)){
        palx <- pal[[s]]
      }

      g <- ggplot(toplot, aes(x=get(coord1), y=get(coord2), color=density)) +
        geom_point_rast(data = dt, size = pt.size, fill = na.color, shape=21, color='black', stroke=0.5) + # include all points and color gray
        geom_point_rast(data = dt, size = pt.size, color = na.color, shape=1) + # include all points and color gray
        geom_point_rast(size = pt.size, stroke=0.5) +
        xlab( gsub("_","",coord1) ) + ylab( gsub("_","",coord2) ) +
        scale_color_gradientn('', guide = "colourbar", colours = palx,
                              limits=c(0,1),
                              breaks=seq(from=0, to=1, by=.2)) +
        ggtitle(paste("Density", gsub("_"," ",s) )) + xlab(coord1) + ylab(coord2) +
        theme_bw(12) +
        theme(panel.grid = element_blank(),
              plot.title = element_text(size=title.size, hjust = 0.5, vjust = -1),
              plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
              plot.margin = margin(0,0,0,0),
              strip.background = element_blank(),
              panel.border = element_blank(),
              panel.background = element_rect(fill = "transparent"), # bg of the panel
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.background = element_blank(),
              legend.key.size = unit(16, units='points'),
              legend.title = element_text(angle = 90, hjust = 0.5),
              legend.text = element_text(size = 12),
              legend.position = 'right',
              legend.justification = "top",
              legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(t=0,r=0,b=0,l=0)) #trbl

              # legend.margin=margin(0,0,0,0),
              # legend.box.margin=margin(-10,10,-10,-30)) #trbl

      if(!is.na(contour.color))
        g <- g +  geom_density_2d(data=toplot, linewidth=.4, color=contour.color)

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

      if(small_axis_title){
        xmin <- dt[,min(get(coord1))]
        #xmin <- xmin -(abs(xmin)*.25)
        xmax <- xmin + abs(xmin)*.3
        xlen <- abs(xmax-xmin)

        ymin <- dt[,min(get(coord2))]
        #ymin <- ymin -(abs(ymin)*.25)
        ymax <- ymin + abs(ymin)*.3
        ylen <- abs(ymax-ymin)

        message("xmax",xmax)
        message("xmin",xmin)
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

      }
    }

  })

  p <- p[lengths(p) > 0]

  panel_width_in <- (panel_width_in) * (length(p)/nrow)
  panel_height_in <- (panel_height_in+.5) * nrow

  if(!missing(filename.prefix)){
    if(as.pdf){
      f <- file.path(plot.dir,paste0(filename.prefix,".pdf"))
      message("Plotting to:",f)
      pdf(file = f, width = panel_width_in, height = panel_height_in, bg = 'transparent')
      #do.call(grid.arrange, c(p, nrow=nrow))
      grid_arrange_list_shared_legend(p, nrow = nrow, legend.position = "topright")
      dev.off()
    } else{
      f <- file.path(plot.dir,paste0(filename.prefix,".png"))
      message("Plotting [",panel_width_in,"x",panel_height_in,"] to:",f)
      png(filename = f, width = panel_width_in, height = panel_height_in, bg = 'transparent', units = 'in', res=200)
      #do.call(grid.arrange, c(p, nrow=nrow))
      grid_arrange_list_shared_legend(p, nrow = nrow, legend.position = "topright")
      dev.off()
    }
  }

  p
}

#' ComplexHeatmap-based replacement for Seurat's DoHeatmap
#'
#' Draws a cell x feature heatmap (z-scored by default) with a group annotation
#' bar, optional clustering of features and groups, and dendrograms.
#'
#' @param object Seurat object.
#' @param features Genes to show.
#' @param group.by Meta.data column used to group cells.
#' @param group.colors Named colours for the groups.
#' @param heatmap.colors Colour ramp of the heatmap body.
#' @return A ComplexHeatmap object.
DoHeatmapEnhance <- function(object,
                             features = NULL,
                             cells = NULL,
                             order.cells.by.desc=NULL,
                             heatmap.colors = heatmap_palette,
                             group.by = "ident",
                             group.bar = TRUE,
                             group.colors = NULL, disp.min = -2.5,
                             disp.max = NULL, layer = "data", assay = "RNA",
                             show.feature.names = TRUE,
                             show.group.dendrogram = TRUE,
                             show.cell.dendrogram = FALSE,
                             feature.fontsize = 8,
                             group.fontsize = 14,
                             legend.fontsize = 9,
                             raster = TRUE, draw.lines = TRUE,
                             lines.width = NULL, disp.identity.legend = TRUE,
                             heatmap.orientation='vertical',
                             legend.direction="vertical",
                             title="",
                             legend.title='Expression',
                             sample.cell.ratio = 1,
                             clustering.within = T,
                             cluster.features = T,
                             no.clustering = F,
                             zscore_data = T,
                             drop.min.value = NA,
                             label.value.min = NA,
                             cluster_group_slices = T,
                             group.bar.height = 0.02)
{
  require(Matrix)
  require(ComplexHeatmap)
  require(grid)
  require(circlize)

  if(is.null(cells)){
    cells <- cells %||% colnames(x = object)
    if (is.numeric(x = cells)) {
      cells <- colnames(x = object)[cells]
    }
  }
  message("Getting features")
  if(is.null(features)){
    features <- features %||% VariableFeatures(object = object)
    features <- rev(x = unique(x = features))
  }
  else{
    message("Retrieving requested features")
  }
  disp.max <- disp.max %||% ifelse(test = layer == "scale.data",
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object,
                                                 layer = layer))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", layer,
           " layer for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ",
            layer, " layer for the ", assay, " assay: ", paste(bad.features,
                                                             collapse = ", "))
    #features <- setdiff(features, bad.features)
  }
  cat(paste(features, collapse = ","))

  if(!is.null(order.cells.by.desc)){
    dt <- data.table(object@meta.data, keep.rownames = 'cells')
    cells <- dt[order(get(order.cells.by.desc), decreasing = T), .SD, by=group.by]$cells
  }

  message("Getting data from assay ",assay," in layer ",layer)
  # matrix[cells, features]
  mat <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object,
                                                            assay = assay,
                                                            layer = layer)[features, cells, drop = FALSE]))) %>% as.matrix()

  group_anno <- factor(object@meta.data[cells,group.by])
  ngroup <- length(unique(group_anno))

  bionames_anno <- factor(object@meta.data[cells,'bioNames'], levels = bionames_order)
  nbionames <- length(unique(bionames_anno))

  if(!is.na(drop.min.value)){
    nc <- ncol(mat)
    # Keep any columns that are > drop.min.value
    # Drops columns in which all values are ≤ drop.min.value
    mat <- mat[,apply(mat,2,function(x) any(x > drop.min.value))]
    ncm <- ncol(mat)
    if(ncm < nc){
      message("removed ",nc-ncm," features with all values ≤ ",drop.min.value)
    }
  }

  if(sample.cell.ratio<1){
    message("Sampling ",100*sample.cell.ratio,"% of cells in each group")
    d <- dim(mat)
    # r <- sample(x = rownames(mat), size = d[1]*sample.cell.ratio)
    r <- sapply( levels(group_anno), function(i){
      sample(x = which(group_anno == i), size = length(which(group_anno == i))*sample.cell.ratio)
    })
    print(table(group_anno))
    print(lengths(r))
    idx <- unlist(r, use.names = F)
    # subset matrix to subsampled cells
    mat <- mat[idx, ]
    # subset cluster annotation to subsampled cells
    group_anno<- factor(object@meta.data[idx,group.by])
    ngroup <- length(unique(group_anno))
    bionames_anno <- factor(object@meta.data[idx,'bioNames'])
    nbionames <- length(unique(bionames_anno))
  } else {
    message("Using all cells")
  }
  message("groups: ", paste(levels(group_anno), collapse = ", "))
  message("bioNames: ", paste(levels(bionames_anno), collapse = ", "))

  fontsizes <- rep(feature.fontsize, ncol(mat))
  fontcolors <- rep('black', ncol(mat))
  fontfaces <- rep('plain',ncol(mat))
  colAnno <- NULL
  if(show.feature.names){
    colAnno <- columnAnnotation(columns = anno_text(colnames(mat), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))
    if(!is.na(label.value.min)){
      col_idx <- which(apply(mat,2,function(x) any(x >= label.value.min)))
      fontsizes[col_idx] <- 8
      fontcolors[col_idx] <- 'black'
      fontfaces[col_idx] <- 'bold'
      # Create text annotation object for displaying row names
      colAnno <- columnAnnotation(columns = anno_text(colnames(mat), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))
    }
  }

  q <- c(min(mat), max(mat))
  #mat <- t(mat)
  if(layer=='scale.data'){
    q <- quantile(mat, c(0.1, 0.95))
    # what's the value range in the matrix
    q <- c(-max(abs(q)),max(abs(q)))
  }
  else if(zscore_data){
    message("Scaling data")
    mat <- scale(mat)
    # what's the value range in the matrix
    q <- quantile(mat, c(0.1, 0.95))
    q <- c(-max(abs(q)), max(abs(q)))
    legend.title <- paste0("scaled\n",legend.title)
  }
  mat <- t(mat)
  message("mat range:",paste(range(mat), collapse = ":"))
  message("lower-upper limit:",paste(round(q,2), collapse = ":"))

  if(is.null(heatmap.colors)){
    heatmap.colors <- colorRampPalette(c("royalblue4","royalblue2","deepskyblue",
                                         "grey90","orangered3","red3","darkred"), space="rgb")(1024);
  }

  col_fun = circlize::colorRamp2(seq(floor(q[1]), ceiling(q[2]), length = length(heatmap.colors)), heatmap.colors)

  dendc <- NULL
  column_title <- NULL
  column_split <- NULL
  show_column_dend <- FALSE
  cluster_rows <- FALSE
  cluster_column_slices <- cluster_group_slices
  show.group.bar.legend <- F
  ## We don't input single cell so clustering within is same as group clustering
  if( (nrow(mat) >= 65000)){
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
      dendc <- cluster_within_group(mat,  group_anno) # both groups and within groups are reordered, better for discerning clusters, but looks rougher
      # For large heatmaps clustering rows with dendrogram works best
      cluster_rows <- as.dendrogram(hclust(dist(mat, method = 'binary')))
      column_split <- ngroup
      message("Finished clustering")
    }
    else{
      message("Performing between group clustering")
      dendc <- cluster_between_groups(mat, group_anno) # only reordered between groups, looks more smooth
      # For smaller heatmaps (with labels) clustering rows based on dendc works better
      cluster_rows <- as.dendrogram(hclust(dist(mat)))
      show_column_names <- FALSE # dendrogram shows these
      column_split <- ngroup
    }
  }
  if(!cluster.features)
    cluster_rows = F

  if(show.group.dendrogram & !show.cell.dendrogram){
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
    column_split <- group_anno
    dendc <- TRUE # Let's heatmap cluster and adds dendrogram for each cell
  }

  ht <- NULL
  message("Generating heatmap in orientation ",heatmap.orientation)
  message("mat[",paste(dim(mat),collapse = ":"),"]")
  if(heatmap.orientation=='vertical'){
    top_anno <- NULL
    if(group.bar)
      if(group.by == 'bioNames'){
        top_anno = HeatmapAnnotation(df = data.frame(group = group_anno),
                                     col = list(group = c(group.colors)),
                                     height = unit(x = group.bar.height, units = 'cm'),
                                     show_annotation_name = F,
                                     show_legend = show.group.bar.legend,
                                     border=T)
      } else {
        top_anno = HeatmapAnnotation(df = data.frame(group = group_anno, bionames = bionames_anno),
                                     col = list(group = c(group.colors), bionames = bionames_palette),
                                     height = unit(x = group.bar.height, units = 'cm'),
                                     show_annotation_name = F,
                                     show_legend = show.group.bar.legend,
                                     border=T)
      }

    ht <- Heatmap(mat,
                  column_title = title,
                  name = legend.title,
                  cluster_columns = dendc,
                  column_split = column_split,
                  show_column_dend = show_column_dend,
                  column_dend_gp = gpar(lwd=0.5),
                  cluster_column_slices = cluster_column_slices,
                  column_title_gp = gpar(fontsize = group.fontsize, fontface='bold'),
                  column_gap = unit(0.5, "mm"),
                  cluster_rows = cluster_rows,
                  show_row_dend = FALSE,
                  col = col_fun,
                  row_names_gp = gpar(fontsize = feature.fontsize),
                  column_names_gp = gpar(fontsize = feature.fontsize),
                  column_title_rot = 90,
                  top_annotation = top_anno,
                  heatmap_legend_param = list(direction = legend.direction, legend_gp = gpar(fontsize = legend.fontsize)),
                  show_column_names = FALSE,
                  show_row_names = show.feature.names,
                  use_raster = raster,
                  border = T,
                  raster_quality = 8)
  }
  else{
    left_anno <- NULL
    if(group.bar){
      if(group.by == 'bioNames'){
        left_anno = HeatmapAnnotation(df = data.frame(group = group_anno),
                                      which='row',
                                      col = list(group = c(group.colors)),
                                      width = unit(x = group.bar.height, units = 'cm'),
                                      show_annotation_name = F,
                                      show_legend = show.group.bar.legend,
                                      border=T)
      } else {
        left_anno = HeatmapAnnotation(df = data.frame(group = group_anno, bionames = bionames_anno),
                                      which='row',
                                      col = list(group = c(group.colors), bionames = bionames_palette),
                                      width = unit(x = group.bar.height, units = 'cm'),
                                      show_annotation_name = F,
                                      show_legend = show.group.bar.legend,
                                      border=T)
      }
    }

    ht <- Heatmap(t(mat),
                  column_title = title,
                  name = legend.title,
                  cluster_rows = dendc,
                  row_split = column_split,
                  show_row_dend = show_column_dend,
                  row_dend_gp = gpar(lwd=0.5),
                  cluster_row_slices = cluster_column_slices,
                  row_title_gp = gpar(fontsize = group.fontsize, fontface='bold'),
                  row_gap = unit(0.5, "mm"),
                  cluster_columns = cluster_rows,
                  show_column_dend = FALSE,
                  col = col_fun,
                  #column_names_gp = gpar(fontsize = feature.fontsize),
                  row_title_rot = 0,
                  left_annotation = left_anno,
                  heatmap_legend_param = list(direction = legend.direction, legend_gp = gpar(fontsize = legend.fontsize)),
                  show_row_names = FALSE,
                  show_column_names = FALSE, #show.feature.names,
                  bottom_annotation = colAnno,
                  use_raster = raster,
                  border = T,
                  raster_quality = 8)
  }
  ht

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

#' Customised version of Seurat's VlnPlot
#'
#' Adds control over palette, axis and strip text sizes, sorting and stacking
#' (via MultiViolinPlotEnhance()). Adapted from Seurat::VlnPlot.
ViolinPlotEnhance <- function (object, features, type = "violin", idents = NULL, ncol = NULL,
                               sort = FALSE, assay = NULL, y.max = NULL, same.y.lims = FALSE,
                               adjust = 1, cols = NULL, pt.size = 0, group.by = NULL, split.by = NULL,
                               log = FALSE, layer = "data", stack = FALSE, combine = TRUE,
                               pal = rocket_palette,
                               y_axis_label="Expression Level (log normalized)",
                               title = NULL,
                               x.axis.text.size=16, y.axis.text.size=12, strip.text.size=18,
                               x_axis_text_angle=0, y_text_hjust=0.5,
                               fill.by = NULL, flip = FALSE, raster = NULL)
{
  require(Seurat)
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  if (isTRUE(x = stack)) {
    if (!is.null(x = ncol)) {
      warning("'ncol' is ignored with 'stack' is TRUE",
              call. = FALSE, immediate. = TRUE)
    }
    if (!is.null(x = y.max)) {
      warning("'y.max' is ignored when 'stack' is TRUE",
              call. = FALSE, immediate. = TRUE)
    }
  }
  else {
    ncol <- ncol %||% ifelse(test = length(x = features) >
                               9, yes = 4, no = min(length(x = features), 3))
  }

  data <- NULL
  # features are genes or other slot
  if(!all(features %in% colnames(object@meta.data)) ){
    data <- FetchData(object = object, vars = features, layer = layer)
  }

  # check if features are in meta.data
  mdata <- NULL
  if( any(features %in% colnames(object@meta.data)) ){
    mfeatures <- features[which(features %in% colnames(object@meta.data))]
    mdata <- subset(object@meta.data, select = mfeatures)
    if(!is.null(data))
      data <- cbind(data, mdata)
    else data <- mdata
  }

  pt.size <- pt.size %||% AutoPointSize(data = object)
  features <- colnames(x = data)
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  }
  else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in%
                                                 idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  }
  else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  }
  else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    else if (length(x = cols) == 1 && cols == "interaction") {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    }
    else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- levels(x = split)
    if ((length(x = cols) > 2) & (type == "splitViolin")) {
      warning("Split violin is only supported for <3 groups, using multi-violin.")
      type <- "violin"
    }
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  if (isTRUE(x = stack)) {
    return(MultiViolinPlotEnhance(type = type, data = data, idents = idents,
                                  split = split, sort = sort, same.y.lims = same.y.lims,
                                  adjust = adjust, cols = cols, pt.size = pt.size,
                                  pal = pal, title = title,
                                  y_axis_label = y_axis_label,
                                  x.axis.text.size = x.axis.text.size,
                                  y.axis.text.size=y.axis.text.size,
                                  x_axis_text_angle=x_axis_text_angle,
                                  y_text_hjust=y_text_hjust,
                                  strip.text.size = strip.text.size,
                                  log = log, fill.by = fill.by, flip = flip))
  }
  plots <- lapply(X = features, FUN = function(x) {
    return(SingleExIPlot(type = type, data = data[, x, drop = FALSE],
                         idents = idents, split = split, sort = sort, y.max = y.max,
                         adjust = adjust, cols = cols, pt.size = pt.size,
                         log = log, raster = raster))
  })
  label.fxn <- switch(EXPR = type, violin = if (stack) {
    xlab
  } else {
    ylab
  }, splitViolin = if (stack) {
    xlab
  } else {
    ylab
  }, ridge = xlab, stop("Unknown ExIPlot type ", type, call. = FALSE))
  for (i in 1:length(x = plots)) {
    key <- paste0(unlist(x = strsplit(x = features[i], split = "_"))[1],
                  "_")
    obj <- names(x = which(x = Key(object = object) == key))
    if (length(x = obj) == 1) {
      if (inherits(x = object[[obj]], what = "DimReduc")) {
        plots[[i]] <- plots[[i]] + label.fxn(label = "Embeddings Value")
      }
      else if (inherits(x = object[[obj]], what = "Assay")) {
        next
      }
      else {
        warning("Unknown object type ", class(x = object),
                immediate. = TRUE, call. = FALSE)
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    }
    else if (!features[i] %in% rownames(x = object)) {
      plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
    }
  }
  if (combine) {
    require(patchwork)
    plots <- wrap_plots(plots, ncol = ncol)
    if (length(x = features) > 1) {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

#' Stacked violin plot of several features (helper of ViolinPlotEnhance)
#'
#' Adapted from Seurat's internal MultiExIPlot.
MultiViolinPlotEnhance <- function (data, idents, split = NULL, type = "violin", sort = FALSE,
                                    same.y.lims = same.y.lims, adjust = 1, pt.size = 0, cols = NULL,
                                    pal = rocket_palette, title = NULL, y_axis_label=NULL,
                                    x.axis.text.size=9, y.axis.text.size=10,
                                    strip.text.size=24, x_axis_text_angle=45,
                                    y_text_hjust=0.5,
                                    seed.use = 42, log = FALSE, fill.by = NULL, flip = NULL)
{
  require(dplyr)
  require(cowplot)
  #theme_set(theme_cowplot())
  if (!(fill.by %in% c("feature", "ident"))) {
    stop("`fill.by` must be either `feature` or `ident`")
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) < 2) {
    stop("MultiExIPlot requires a data frame with >1 column")
  }
  #data <- Melt(x = data)
  data <- Seurat:::Melt(x = data)
  data <- data.frame(feature = data$cols, expression = data$vals,
                     ident = rep_len(x = idents, length.out = nrow(x = data)))
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$feature <- as.vector(x = data$feature)
    data$ident <- as.vector(x = data$ident)
    avgs.matrix <- sapply(X = split(x = data, f = data$ident),
                          FUN = function(df) {
                            return(tapply(X = df$expression, INDEX = df$feature,
                                          FUN = mean))
                          })
    idents.order <- hclust(d = dist(x = t(x = L2Norm(mat = avgs.matrix,
                                                     MARGIN = 2))))$order
    avgs.matrix <- avgs.matrix[, idents.order]
    avgs.matrix <- L2Norm(mat = avgs.matrix, MARGIN = 1)
    position <- apply(X = avgs.matrix, MARGIN = 1, FUN = which.max)
    mat <- hclust(d = dist(x = avgs.matrix))$merge
    orderings <- list()
    for (i in 1:nrow(mat)) {
      x <- if (mat[i, 1] < 0)
        -mat[i, 1]
      else orderings[[mat[i, 1]]]
      y <- if (mat[i, 2] < 0)
        -mat[i, 2]
      else orderings[[mat[i, 2]]]
      x.pos <- min(x = position[x])
      y.pos <- min(x = position[y])
      orderings[[i]] <- if (x.pos < y.pos) {
        c(x, y)
      }
      else {
        c(y, x)
      }
    }
    features.order <- orderings[[length(x = orderings)]]
    data$feature <- factor(x = data$feature, levels = unique(x = sort(x = data$feature))[features.order])
    data$ident <- factor(x = data$ident, levels = unique(x = sort(x = data$ident))[rev(x = idents.order)])
  }
  else {
    data$feature <- factor(x = data$feature, levels = unique(x = data$feature))
  }
  if (log) {
    noise <- rnorm(n = nrow(x = data))/200
    data$expression <- data$expression + 1
  }
  else {
    noise <- rnorm(n = nrow(x = data))/1e+05
  }
  for (f in unique(x = data$feature)) {
    if (all(data$expression[(data$feature == f)] == data$expression[(data$feature ==f)][1])) {
      warning("All cells have the same value of ", f, call. = FALSE,
              immediate. = TRUE)
    }
    else {
      data$expression[(data$feature == f)] <- data$expression[(data$feature ==f)] + noise[(data$feature == f)]
    }
  }
  if (type == "violin" && !is.null(x = split)) {
    data$split <- rep_len(x = split, length.out = nrow(data))
    vln.geom <- geom_violin
    fill.by <- "split"
  }
  else if (type == "splitViolin" && !is.null(x = split)) {
    data$split <- rep_len(x = split, length.out = nrow(data))
    vln.geom <- geom_split_violin
    fill.by <- "split"
    type <- "violin"
  }
  else {
    vln.geom <- geom_violin
  }
  switch(EXPR = type,
         violin = {
           geom <- list(vln.geom(scale = "width", adjust = adjust, trim = TRUE, linewidth = .1)) # size is border thickness of violin
         }, ridge = {
           geom <- list(geom_density_ridges(scale = 4), theme_ridges(),
                        scale_y_discrete(expand = c(0.01, 0)))
         }, stop("Unknown plot type: ", type))
  if (flip) {
    x <- "ident"
    x.label <- "Identity"
    y <- "expression"
    y.label <- y_axis_label #"Expression Level (log normalized)"
  }
  else {
    y <- "ident"
    y.label <- "Identity"
    x <- "expression"
    x.label <- y_axis_label #"Expression Level (log normalized)"
  }

  ### New here start
  setDT(data)
  data[,mean_expression:=mean(expression), by=list(feature, ident)]
  #View(data)
  plot <- ggplot(data = data, mapping = aes(x = data[[x]], y = as.numeric(data[[y]]), fill = mean_expression)[c(2, 3, 1)]) +
    labs(x = x.label, y = y.label, fill = NULL) + theme_cowplot()

  plot <- do.call(what = "+", args = list(plot, geom))
  ### New here end

  if (flip) {
    hjust <- 0.5
    if(x_axis_text_angle != 0)
      hjust <- 1
    plot <- plot + scale_y_continuous(expand = c(0, 0), labels = function(x) c(rep(x = "",times = length(x) - 2), x[length(x) - 1], "")) +
      facet_grid(feature ~ ., scales = (if (same.y.lims)
        "fixed"
        else "free")) + Seurat:::FacetTheme(panel.spacing = unit(0,"lines"),
                                            panel.background = element_rect(fill = NA,color = "black"),
                                            axis.text.y = element_text(size = y.axis.text.size, angle=0),
                                            axis.text.x = element_text(angle = x_axis_text_angle, vjust=0.5, hjust = hjust, size=x.axis.text.size),
                                            plot.title = element_text(hjust=0.5),
                                            strip.text.y.right = element_text(angle = 0, size=strip.text.size, hjust = y_text_hjust))
  }
  else {
    plot <- plot + scale_x_continuous(expand = c(0, 0), labels = function(x) c(rep(x = "",times = length(x) - 2), x[length(x) - 1], "")) +
      facet_grid(. ~ feature, scales = (if (same.y.lims)
        "fixed"
        else "free")) + Seurat:::FacetTheme(panel.spacing = unit(0,"lines"),
                                            panel.background = element_rect(fill = NA, color = "black"),
                                            axis.text.y = element_text(size = y.axis.text.size),
                                            axis.text.x = element_text(size = x.axis.text.size),
                                            plot.title = element_text(hjust=0.5),
                                            strip.text.x = element_text(angle = -90, size=strip.text.size))
  }
  if (log) {
    plot <- plot + scale_x_log10()
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      }
      else {
        unlist(x = lapply(X = idents, FUN = function(pattern,
                                                     x) {
          x.mod <- gsub(pattern = paste0(pattern, "."),
                        replacement = paste0(pattern, ": "), x = x,
                        fixed = TRUE)
          x.keep <- grep(pattern = ": ", x = x.mod, fixed = TRUE)
          x.return <- x.mod[x.keep]
          names(x = x.return) <- x[x.keep]
          return(x.return)
        }, x = unique(x = as.vector(x = data$split))))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    }
    else {
      #labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  } else {

    sum_func <- list(Min = min, Mean = mean, Max = max)

    ls_val <- unlist(lapply(sum_func, function(f) f(data$mean_expression)))

    # colors the violins across levels
    if(ls_val[1] > 0)
      ls_val[1] <- 0

    plot <- plot + scale_fill_gradientn(
      colours =  pal,
      limits=c(ls_val[1], ls_val[3]),
      breaks = c(ls_val[1],ls_val[3])
    )

  }

  plot <-  plot + ggtitle(title)

  #print(plot)

  #return(data)

  return(plot)
}

#' Stacked bar chart of a meta.data variable per group
#'
#' Shows, for each level of `by.group` (e.g. cluster), the percentage of cells
#' belonging to each level of `variable` (e.g. sample), optionally annotated
#' with cell counts.
#'
#' @param seuratObj Seurat object.
#' @param variable Meta.data column plotted as stacked fill.
#' @param by.group Meta.data column giving the bars.
#' @param palette Named colours for `variable`.
#' @param filename Output file; if empty the plot is only returned.
PlotBarVariableInClusters <- function(seuratObj,
                                      variable='orig.ident',
                                      by.group='seurat_clusters',
                                      legend.title="",
                                      legend.column=1,
                                      palette=sample_palette,
                                      addcounts=T,
                                      orientation='horizontal',
                                      width_multiplier=1,
                                      height_multiplier=1,
                                      text.angle=0,
                                      x.label,
                                      y.label,
                                      filename=""){
  dt <- as.data.table(seuratObj@meta.data[,c(eval(by.group),eval(variable))])
  setnames(dt, by.group, 'grouping')

  dt[,nCluster := .N, by=grouping]
  dt[,nClustervar := .N, by=list(grouping, get(variable))]
  dt <- unique(dt)
  dt <- dt[order(grouping, decreasing = F)]
  dt[,percClustervar := 100 * (nClustervar / nCluster) ]
  dtlab <- unique(dt[,list(grouping, nCluster)])

  fill <- scale_fill_manual(legend.title, values = palette)

  g <- NULL
  width <- 9
  height <- 4
  if(missing(x.label))
    x.label <- by.group
  if(missing(y.label))
    y.label <- "Percentage of cells"
  if(orientation=='horizontal'){
    g <- ggplot(dt, aes(x=grouping, y=percClustervar)) +
      geom_bar(stat='identity', aes(fill=get(variable)), linewidth = .95) +
      fill +
      xlab(x.label) +
      ylab(y.label) +
      theme_bw(14) +
      coord_cartesian(ylim = c(0,100), expand = F, clip = 'off') +
      theme(panel.ontop = T,
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major.y = element_line(colour = 'black', size=.2, linetype = 'longdash'),
            legend.direction = "vertical",
            legend.position = "right",
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = text.angle, hjust=0.5, vjust=0.5),
            axis.title = element_text(size = 12),
            text = element_text(size = 12),
            legend.text = element_text(size = 11),
            legend.background = element_blank(),
            legend.margin = margin(0,0,0,1, unit = "cm"), #trbl
            plot.margin = unit(x = c(2,1,1,1), units = 'cm'),
            plot.background = element_blank())
    g <- g + guides(fill = guide_legend(title = legend.title, title.position = "top", ncol = legend.column))

    if(addcounts)
      g <- g + geom_text(data = dtlab, aes(label=paste0("n=",nCluster)), y=100, color='black',
                         vjust=0, hjust=0, angle=90, size=3.5)
  } else {
    width <- 4
    height <- 4
    if(missing(x.label))
      x.label <- "Percentage of cells"
    if(missing(y.label))
      y.label <- by.group
    dt <- dt[order(grouping, decreasing = T)]
    dt[,grouping:=factor(grouping, unique(grouping))]
    dtlab <- unique(dt[,list(grouping, nCluster)])
    g <- ggplot(dt, aes(y=grouping, x=percClustervar)) +
      geom_bar(stat='identity', aes(fill=get(variable)), linewidth = .95) +
      fill +
      xlab(x.label) +
      ylab(y.label) +
      theme_bw(14) +
      coord_cartesian(xlim = c(0,100), expand = F, clip = 'off') +
      theme(panel.ontop = T,
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major.x = element_line(colour = 'black', linewidth=.2, linetype = 'longdash'),
            legend.direction = "horizontal",
            legend.position = "right",
            axis.text.x = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 11),
            legend.margin = margin(0,0,0,1, unit = "cm"), #trbl
            legend.background = element_blank(),
            plot.margin = unit(x = c(1,2,1,1), units = 'cm'),
            plot.background = element_blank())
    g <- g + guides(fill = guide_legend(title = legend.title, title.position = "top", ncol = legend.column))

    if(addcounts)
      g <- g + geom_text(data = dtlab, aes(label=paste0("n=",nCluster)), x=101, color='black',
                         vjust=0.5, hjust=0, angle=0, size=3.5)
  }

  if(filename !=""){
    fout <- file.path(plot.dir,paste0(filename,'.pdf'))
    message("Plotting to:",fout)
    pdf(file = fout, width = width * width_multiplier, height = height * height_multiplier, bg = 'transparent')

    print(g)
    dev.off()
  }
  g

}

#' Customised version of Seurat's DotPlot
#'
#' Supports named feature groups, horizontal/vertical orientation and tunable
#' font and legend sizes. Adapted from Seurat::DotPlot.
MyDotPlot <- function (object,
                       features,
                       assay = NULL,
                       cols = c("lightgrey", "blue"),
                       font.base.size = 14,
                       legend.title.size = 11,
                       plot.orientation = "Horizontal",
                       legend.box = "horizontal",
                       col.min = -2.5, col.max = 2.5, dot.min = 0.0, dot.scale = 6,
                       pt.stroke.size = 0.1,
                       panel.grid.stroke.size = 0.05,
                       idents = NULL, group.by = NULL, cluster.idents = FALSE,
                       scale = TRUE, scale.min = 0, scale.max = 1, # scale color bar
                       scale.by = "radius", scale.min.pct = NA, scale.max.pct = NA # scale point size
)
{
  #@object is a seurat object

  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features),
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE,
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }

  dt.feature.groups <-   data.table(do.call(rbind,as.list(feature.groups)), keep.rownames = "feature")
  setnames(dt.feature.groups,c("features.plot","feature.groups"))

  cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]),
                                        idents = idents))
  data.features <- FetchData(object = object,
                             vars = features,
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)

  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot,
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled",
            call. = FALSE, immediate. = TRUE)
  }
  else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results",
            call. = FALSE, immediate. = TRUE)
  }
  library(scales)
  message("Setting data to log1p")
  legend.exp.title <- ifelse(scale,"Mean Expression\nper group scaled","Average\nExpression")
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
                             # scales expression of gene gene across clusters (samples, etc)
                             if (scale) {
                               data.use.log <- log1p(data.use)
                               if(!is.na(scale.min) | !is.na(scale.max)){
                                 col.min <- ifelse(!is.na(scale.min), scale.min, min(data.use.log))
                                 col.max <- ifelse(!is.na(scale.max), scale.max, max(data.use.log))
                                 data.use <- rescale(data.use.log, to=c(col.min, col.max))
                               } else {
                                 data.use <- scale(x = data.use.log)
                                 #data.use <- MinMax(data = data.use, min = col.min, max = col.max)
                               }
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))

  data.plot$avg.exp.scaled <- avg.exp.scaled
  dups <- which(duplicated(features))
  if(length(dups)==0)
    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features)
  else message("There are duplicated features:",paste(features[dups], collapse = ","))

  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100

  color.by <- "avg.exp.scaled"

  if (!is.na(x = scale.min.pct)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min.pct
  }
  if (!is.na(x = scale.max.pct)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max.pct
  } else {
    scale.max <- ceiling(max(data.plot[,"pct.exp"]))
  }
  if (!is.null(x = feature.groups)) {
    feature.groups.sort <- unique(feature.groups)

    # This allows for genes to exist in more than one category
    setDT(data.plot, key = "features.plot")
    setkey(x = dt.feature.groups,"features.plot")
    data.plot <- merge(data.plot, dt.feature.groups, allow.cartesian=TRUE)
    data.plot[, feature.groups := factor(x = feature.groups, levels = feature.groups.sort)]
  }
  data.plot <- data.plot[order(-data.plot$avg.exp),]
  data.plot$rn <- 1:nrow(data.plot)

  if(tolower(plot.orientation) %in% c('h','horizontal'))
    plot.orientation <- 'horizontal'
  else
    plot.orientation <- 'vertical'

  plot <- NULL
  if(plot.orientation=='horizontal')
    plot <- ggplot(data = data.plot, aes(x = reorder(features.plot, rn), y = id))
  else
    plot <- ggplot(data = data.plot, aes(x = reorder(features.plot, -rn), y = id))

  plot <- plot +
    geom_point(mapping = aes_string(size = "pct.exp", fill = color.by),
               color='grey50', stroke = pt.stroke.size, shape = 21) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min.pct, scale.max.pct)) +
    guides(size = guide_legend(title = "(%) of Cells\nExpressing"), theme = theme(fill = "black")) +
    labs(x = "Features", y = "Identity") +
    theme_bw(base_size = font.base.size) +
    theme(legend.title = element_text(size=legend.title.size),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=90))

  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(rows = ~feature.groups,
                              scales = "free_x",
                              space = "free_x",
                              switch = "y")
    plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                         panel.spacing = unit(0, "pt"),
                         legend.key.height = unit(1, "pt"),
                         panel.grid = element_line(linewidth = panel.grid.stroke.size),
                         panel.border = element_rect(color = "black", fill = NA, size = .5),
                         strip.background = element_blank(),
                         plot.background = element_blank(),
                         panel.background = element_blank())
  }

  if (length(x = cols) == 1) {
    plot <- plot + scale_fill_distiller(palette = cols)
  }
  else if (length(x = cols) == 2) {
    plot <- plot + scale_fill_gradient(low = cols[1], high = cols[2])
  }
  else{
    plot <- plot + scale_fill_gradientn(colours = cols)
  }
  plot <- plot + guides(fill = guide_colorbar(title = legend.exp.title,
                                              barheight = 4,
                                              direction = "vertical",
                                              frame.colour = "black",    # Set the frame color
                                              frame.linewidth = .5,       # Set the frame line width
                                              frame.linetype = "solid"   # Set the frame line type
  )) +  theme(legend.box = legend.box,
              # Spacing between different legends
              legend.spacing.y = unit(.05, "cm"),
              legend.spacing.x = unit(.25, "cm"))

  if(plot.orientation=='vertical'){
    plot <- plot + coord_flip() + facet_grid(feature.groups ~., scale="free", space="free") +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            strip.text.y = element_text(angle = 0),
            panel.grid = element_line(linewidth = panel.grid.stroke.size),
            legend.position = "right",
            legend.justification = "top")

  }
  if(panel.grid.stroke.size==0){
    plot <- plot + theme(panel.grid = element_blank())
  }

  return(plot)
}

#' Violin + box plots of QC metrics per sample
#'
#' @param seuratObj Seurat object with a `bioNames` meta.data column.
#' @param features QC metrics to plot, one panel each.
#' @param pal Named colours per sample label.
#' @return Invisibly, the list of ggplots.
ViolinPlotQC <- function(seuratObj,
                         features=c('nFeature_RNA','nCount_RNA','percent_mito','percent_ribo'),
                         pal=bionames_palette,
                         title="", ylim=NULL){
  gp <- lapply(features, function(i){
    ggplot(seuratObj@meta.data, aes(x=bioNames, y=get(i), fill=bioNames)) +
      geom_violin(trim = F, scale = 'width') +
      scale_fill_manual("", values = pal) +
      geom_boxplot(fill='white', alpha=0.7, width=0.5, outlier.size = 1) +
      ggtitle(title) + xlab("") + ylab(i) + scale_y_continuous(limits = ylim) +
      theme_classic(14) + theme(axis.text.x = element_text(angle=90), legend.position = 'none')
  })
  #do.call(grid.arrange, c(gp, nrow=2))
}
