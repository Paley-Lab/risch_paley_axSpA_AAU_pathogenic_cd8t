setwd('/storage1/fs1/paleym/Active/Isabel_Risch/MP024/')

library(Seurat)
library(dplyr)
library(ggplot2)

# Create a variable to specify the folder that contains the CellRanger outputs for each of the samples
indir <- "/storage1/fs1/paleym/Active/SeqData/MP164/cellranger_v8.0.1/"

# get the names of all the samples that exist within that folder
sampleNames <- list.dirs(indir, recursive = F, full.names = F)

# create an empty list, with each "slot" named for a sample, so we can store stuff there in a second
storeSamplesHere <- vector("list", length(sampleNames))
names(storeSamplesHere) <- sampleNames

# Now, because I know the way CellRanger structures its output directories, 
# I can loop through the samples one-by-one, read in the CellRanger output for each sample, and then
# store that sample in our list we've created, under the correct name
for (i in sampleNames) {
  getDataFromHere <- paste0(indir, i, "/outs/per_sample_outs/", i, "/count/sample_filtered_feature_bc_matrix")
  fullSample <- Seurat::Read10X(getDataFromHere) 
  ## Create the GEX object, don't do feature filtering here or it messes up the HTO assay
  storeSamplesHere[[i]] <- CreateSeuratObject(fullSample[[1]], project = i, min.cells = 1, min.features = 0)
  ## Add the antibody hashtag data
  storeSamplesHere[[i]][["HTO"]] <- CreateAssayObject(fullSample[[2]], min.cells = 0, min.features = 0)
}

#cleanup
rm(fullSample); rm(getDataFromHere); rm(i)

## merge all the data into one object
tot <- merge(storeSamplesHere[[1]], y=storeSamplesHere[2:length(storeSamplesHere)], add.cell.ids = sampleNames)

####################################
### Filter out low-quality cells ###
####################################

# make sure the default assay is RNA
DefaultAssay(tot) 

#### Do conventional QC on these cells ####
tot[["percent.mt"]] <- PercentageFeatureSet(tot, pattern = "^MT-")

VlnPlot(tot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggplot2::ggsave(filename = "plot/QC/pre_qc_vlnplots_basicStats.pdf", width=7, height=5)

plot1 <- FeatureScatter(tot, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", shuffle = T, pt.size = 0.3)
plot2 <- FeatureScatter(tot, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", shuffle = T, pt.size = 0.3)
plot1 + plot2
ggplot2::ggsave(filename = "plot/QC/pre_qc_scatterplots.pdf", width=15, height=8)

tot <- subset(tot, subset = nFeature_RNA > 200 & percent.mt < 7.5 & nCount_RNA < 15000)

VlnPlot(tot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggplot2::ggsave(filename = "plot/QC/post_qc_vlnplots_basicStats.pdf", width=8, height=5)

table(tot$orig.ident)
#  Day1   Day2 
#  21203  23339 

saveRDS(tot, file = "data/robjects/1_filtered_unclustered.Rds")

##########################
#### Demuxmix package ####
##########################

# export the HTO counts and the UMIs
hto <- as.data.frame(GetAssayData(tot, slot = "counts", assay = "HTO"))
saveRDS(hto, file = "data/robjects/hto_raw.Rds")
meta <- tot@meta.data
saveRDS(meta, file = "data/robjects/metadata_tot.Rds")

### Now run the script "02_demuxmix.R" in an R env that can handle it (I ran it locally) ###

### Now bring that data back into this environment so we can integrate it into the final product
demux <- readRDS("demultiplexing/demuxmix/classes.Rds")

colnames(demux) <- paste0(colnames(demux), "_demuxmix")
tot <- AddMetaData(tot, metadata = demux)

saveRDS(tot, file = "data/robjects/2_filtered_unclustered_demuxmixed.Rds")

#############################################
##### Demux plotting, Demuxmix pipeline #####
#############################################

Idents(tot) <- "Type_demuxmix"
VlnPlot(tot, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
ggsave(filename = "demultiplexing/demuxmix/vlnplot_UMI_by_prediction.pdf", width=8, height=6)


# Remove negative cells from the object
tot.subset <- subset(tot, subset = (Type_demuxmix != "negative") & (Type_demuxmix != "uncertain"))

# Calculate a tSNE embedding of the HTO data
DefaultAssay(tot.subset) <- "HTO"
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
tot.subset <- NormalizeData(tot.subset, assay = "HTO", normalization.method = "CLR", margin = 2)
# Setting margin = 2 means that we're basing our percentiles on within-cell measurements
# as opposed to within-feature measurements
tot.subset <- ScaleData(tot.subset, features = rownames(tot.subset),
                        verbose = FALSE)
tot.subset <- RunPCA(tot.subset, features = rownames(tot.subset), approx = FALSE, npcs=5)
tot.subset <- RunTSNE(tot.subset, dims = 1:5, perplexity = 100, check_duplicates=F)
tot.subset$plot_dmm <- ifelse(tot.subset$Type_demuxmix %in% "singlet", tot.subset$HTO_demuxmix, tot.subset$Type_demuxmix)
DimPlot(tot.subset, group.by = "HTO_demuxmix") + NoLegend()
ggsave(filename = "demultiplexing/demuxmix/HTO_classification_basic_tSNE.pdf", width=10, height=8)
DimPlot(tot.subset, group.by = "plot_dmm", label = T)
ggsave(filename = "demultiplexing/demuxmix/HTO_classification_tSNE.pdf", width=10, height=8)
# DimPlot(tot.subset, group.by = "Type_dmm_YeiH")
# ggsave(filename = "demultiplexing/demuxmix/HTO_classification_basic_tSNE_YeiH.pdf", width=10, height=8)
# DimPlot(tot.subset, group.by = "HTO_dmm_YeiH")
# ggsave(filename = "demultiplexing/demuxmix/HTO_classification_tSNE_YeiH.pdf", width=11, height=8)

################################################
#### Finalize obj based on demuxmix results ####
################################################

# Extract the singlets
tot.singlet <- subset(tot, subset = Type_demuxmix %in% "singlet")
tot.singlet
# An object of class Seurat 
# 27459 features across 40836 samples within 2 assays 
# Active assay: RNA (27446 features, 0 variable features)
# 1 other assay present: HTO

# Change default assay to RNA for downstream clustering
DefaultAssay(tot.singlet) <- "RNA"

# now go forth and cluster your RNA-seq data etc.
saveRDS(tot.singlet, file = "data/robjects/2_singlets_unclustered.Rds")


######################################
### Normalize and cluster the data ###
######################################
library(Seurat)
library(dplyr)
library(ggplot2)

# read in the data
tot <- readRDS(file = "data/robjects/2_singlets_unclustered.Rds")
# add metadata
tot$subject <- sapply(strsplit(tot$HTO_demuxmix, "-"), "[", 1)
tot$tissue <- sapply(strsplit(tot$HTO_demuxmix, "-"), "[", 3)

table(tot$tissue, tot$orig.ident)

# First-pass analysis: remove bystanders
tot <- subset(tot, subset = HTO_demuxmix != "036-001-PBMC-TotalSeqC")

tot
# An object of class Seurat
# 27459 features across 3294 samples within 2 assays
# Active assay: RNA (27446 features, 0 variable features)
# 1 other assay present: HTO



source("../MP011/quietTCRgenes.R")

#normalize RNA data
tot <- NormalizeData(tot, assay = "RNA")
# Find and scale variable features
tot <- FindVariableFeatures(tot, selection.method = "vst") %>% quietTCRgenes() ## use this Trex function to mask the TCR genes for clustering purposes
tot <- ScaleData(tot, features = rownames(tot))
tot <- RunPCA(tot, verbose = FALSE)
ElbowPlot(tot, ndims = 30) 
tot <- FindNeighbors(tot, dims = 1:20)
tot <- FindClusters(tot, resolution = 0.5, verbose = FALSE)
tot <- RunUMAP(tot, dims = 1:20)

#remove some stray contamination
tot <- subset(tot, idents=c("0", "1", "2", "3", "4", "5", "6", "8", "9"))

# Find and scale variable features
tot <- FindVariableFeatures(tot, selection.method = "vst") %>% quietTCRgenes() ## use this Trex function to mask the TCR genes for clustering purposes
tot <- ScaleData(tot, features = rownames(tot))
tot <- RunPCA(tot, verbose = FALSE)
ElbowPlot(tot, ndims = 30)
tot <- FindNeighbors(tot, dims = 1:20)
tot <- FindClusters(tot, resolution = 0.5, verbose = FALSE)
tot <- RunUMAP(tot, dims = 1:20)

DimPlot(tot)

DefaultAssay(tot) <- "HTO"
tot <- NormalizeData(tot, normalization.method = "CLR", margin = 2)

DefaultAssay(tot) <- "RNA"

# 1 and 2 have AAU
# 2, 3, 5 have Crohn's
tot$AAU <- grepl("001|002", tot$subject)
tot$IBD <- ifelse(grepl("001", tot$subject), "Family Hx of IBD", "Crohn's")

saveRDS(tot, file = "data/robjects/3_singlets_clustered_noBystanders.Rds")

####################
### Add TCR data ###
####################
library(scRepertoire)
rm(list=ls())

d <- readRDS("data/robjects/3_singlets_clustered_noBystanders.Rds")

# Create a variable to specify the folder that contains the CellRanger outputs for each of the samples
indir <- "/Volumes/Active/SeqData/MP164/cellranger_v8.0.1/"

# get the names of all the samples that exist within that folder
sampleNames <- list.dirs(indir, recursive = F, full.names = F)

tcrfiles <- lapply(sampleNames, function(i){tmp <- paste0(indir, i, "/outs/per_sample_outs/", i, "/vdj_t/filtered_contig_annotations.csv")})
tcrfiles <- unlist(tcrfiles)

# Now, because I know the way CellRanger structures its output directories, 
# I can loop through the samples one-by-one, read in the CellRanger output for each sample, and then
# store that sample in our list we've created, under the correct name
all <- list()
for (i in tcrfiles) {
  all[[i]] <- read.csv(i)
}
names(all) <- sampleNames

#make barcodes unique
for (i in 1:length(all)) {
  all[[i]]$barcode <- paste0(names(all)[i], "_", all[[i]]$barcode)
  all[[i]]$contig_id <- paste0(names(all)[i], "_", all[[i]]$contig_id)
}

all <- dplyr::bind_rows(all)

contig.list <- createHTOContigList(all,
                                   d,
                                   group.by = "subject")

combined.TCR <- combineTCR(contig.list, 
                           removeNA = FALSE, 
                           removeMulti = F, 
                           filterMulti = T)
names(combined.TCR) <- names(contig.list)
for (i in names(combined.TCR)) combined.TCR[[i]]$subject <- i

## Add the data to Seurat
sce <- combineExpression(combined.TCR, 
                         d, 
                         cloneCall="strict", 
                         group.by = "subject",
                         cloneSize = c(Single = 1, Small = 5, 
                                       Medium = 50, Large = 150),
                         proportion = F)

DimPlot(sce, group.by = "cloneSize")

# see if there are pathogenic cells here by our loosest possible definition
sce$pathogenic <- ifelse(grepl("^.{6,7}[Y|F]ST.{5}", sapply(strsplit(sce$CTaa, "_"), "[", 2)), TRUE, FALSE)
sce$pathogenic <- ifelse(grepl("TRAV21", sce$CTstrict), sce$pathogenic, FALSE)

# look for pathogenic TCRs
cells <- rownames(sce@meta.data[sce$pathogenic,])
DimPlot(sce, cells.highlight = cells)


saveRDS(sce, file = "data/robjects/4_filtered_clustered_tcr.Rds")

################################################################################
### Harmonize these samples to see if we can get a nicer clustering ###
################################################################################
rm(list=ls())
library(harmony)

d <- readRDS("data/robjects/4_filtered_clustered_tcr.Rds")

d <- d %>% 
  RunHarmony("subject", plot_convergence = TRUE)
d <- d %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

# do a rough DE just to locate trash clusters
de <- FindAllMarkers(d)

DimPlot(d, label = TRUE)
DimPlot(d, group.by = "orig.ident")
FeaturePlot(d, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "CD4", "CD3E", "CD8A", "PTPRC", "CD79A", "CD14"), ncol=3)
ggplot2::ggsave(filename = "plot/QC/basic_qc_umap.pdf", width=10, height=9)
FeaturePlot(d, features = c("CD8A", "CD8B", "KLRB1", "LEF1", "GZMA", "MKI67"), ncol=3)
FeaturePlot(d, features = c("CD27", "TCF7", "CCR7", "IRF7", "PDCD1", "HAVCR2"), ncol=3)

new.cluster.ids <- c(`0` = "KLR_KIR_CD8T", 
                     `1` = "Naive_and_Memory_CD8T", 
                     `2` = "KLR_KIR_CD8T",
                     `3` = "GZMK_CD8T", 
                     `4` = "KLR_KIR_CD8T", 
                     `5` = "KLRB1_CD8T", 
                     `6` = "GNLY_CD8T", 
                     `7` = "Naive_and_Memory_CD8T")
d <- RenameIdents(d, new.cluster.ids)
d$celltype1 <- d@active.ident
DimPlot(d, reduction = "umap", label = F) 

d$tissue_merged <- ifelse(grepl("Blood", d$tissue), "Blood", "Gut")

DefaultAssay(d) <- "HTO"
yeih <- d@assays[["HTO"]]@data[c('YeiH-APC-TotalSeqC', 'YeiH-PE-TotalSeqC'),]
cells <- yeih[,intersect(which(yeih[1,] > 1.5), which(yeih[2,] > 1.5))] %>% colnames()

FeatureScatter(d, feature1 = 'YeiH-APC-TotalSeqC', feature2 = 'YeiH-PE-TotalSeqC', cells = cells)

DimPlot(d, cells.highlight = cells)

saveRDS(d, file = "data/robjects/5_final_harmonized.Rds")

