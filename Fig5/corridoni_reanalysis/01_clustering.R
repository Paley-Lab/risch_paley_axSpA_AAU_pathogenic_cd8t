setwd("/Volumes/Active/Isabel_Risch/MP004/")

########### COLONIC T CELL TRANSCRIPTIONAL ATLAS ##########
# https://www.nature.com/articles/s41591-020-1003-4

library(Seurat)
library(dplyr)

loadInGEX <- function(sampleName = NULL) {
  tmp <- Seurat::Read10X(paste0("/Volumes/Active/SeqData/public/GSE148837_corridoni_2020/cellranger_v8.0.1/experiment1/", sampleName, "/outs/per_sample_outs/", sampleName, "/count/sample_filtered_feature_bc_matrix/")) %>% 
    CreateSeuratObject(project = sampleName, min.cells = 1, min.features = 200)
  return(tmp)
}
## Read in the GEX + hashtag data
s21 <- loadInGEX("S21")
s22 <- loadInGEX("S22")
s23 <- loadInGEX("S23")
s24 <- loadInGEX("S24")
s33 <- loadInGEX("S33")
s34 <- loadInGEX("S34")

## merge all the data into one obj
tot <- merge(s21, y=c(s22, s23, s24, s33, s34), add.cell.ids = c("s21", "s22", "s23", "s24", "s33", "s34"))

####################################
### Filter out low-quality cells ###
####################################

DefaultAssay(tot) # RNA

#### Do conventional QC on these cells ####
tot[["percent.mt"]] <- PercentageFeatureSet(tot, pattern = "^MT-")

VlnPlot(tot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggplot2::ggsave(filename = "plot/QC/pre_qc_vlnplots_basicStats.pdf", width=8, height=5)

plot1 <- FeatureScatter(tot, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", shuffle = T, pt.size = 0.3)
plot2 <- FeatureScatter(tot, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", shuffle = T, pt.size = 0.3)
plot1 + plot2
ggplot2::ggsave(filename = "plot/QC/pre_qc_scatterplots.pdf", width=10, height=8)

tot <- subset(tot, subset = nFeature_RNA > 500 & percent.mt < 15 & nCount_RNA < 10000)

VlnPlot(tot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggplot2::ggsave(filename = "plot/QC/post_qc_vlnplots_basicStats.pdf", width=8, height=5)

table(tot$orig.ident)
#  S21  S22  S23  S24  S33  S34 
# 995 1689 1850 1836 1438 1791 

saveRDS(tot, file = "data/experiment1/robjects/1_filtered_unclustered.Rds")

######################################
### Normalize and cluster the data ###
######################################

#normalize RNA data
tot <- NormalizeData(tot, assay = "RNA")
# Find and scale variable features
tot <- FindVariableFeatures(tot, selection.method = "vst") %>% Trex::quietTCRgenes() ## use this Trex function to mask the TCR genes for clustering purposes
tot <- ScaleData(tot, features = rownames(tot), vars.to.regress = "percent.mt")
tot <- RunPCA(tot, verbose = FALSE)
ElbowPlot(tot, ndims = 30) 
tot <- FindNeighbors(tot, dims = 1:30)
tot <- FindClusters(tot, resolution = 0.5, verbose = FALSE)
tot <- RunUMAP(tot, dims = 1:30)

tot <- subset(tot, idents=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                            "12", "13"))

DimPlot(tot, label = TRUE)
DimPlot(tot, group.by = "orig.ident")
FeaturePlot(tot, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "CD4", "CD3E", "CD8A"), ncol=3)
ggplot2::ggsave(filename = "plot/QC/basic_qc_umap.pdf", width=10, height=6)
FeaturePlot(tot, features = c("FOXP3", "SELL", "TCF7", "CD44", "CCR7"), ncol=3)

tot <- FindClusters(tot, resolution = 1.2, verbose = FALSE)

new.cluster.ids <- c(`0` = "Trm", 
                     `1` = "GZMK_effectors_1", 
                     `2` = "Trm", 
                     `3` = "IL26",
                     `4` = "CD8_CD4", 
                     `5` = "MAIT", 
                     `6` = "TYROBPneg_IELs_2", 
                     `7` = "Naive_1",
                     `8` = "TYROBPneg_IELs_1", 
                     `9` = "IL26",
                     `10` = "TYROBPpos_IELs", 
                     `11` = "GZMK_effectors_2", 
                     `12` = "FGFBP2", 
                     `13` = "Naive_2",
                     `14` = "Trm", 
                     `15` = "GZMK_effectors_1", 
                     `16` = "CD8_CD4_FOXP3", 
                     `17` = "Naive_1",
                     `18` = "GZMK_effectors_2")
tot2 <- RenameIdents(tot, new.cluster.ids)
DimPlot(tot2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

tot2$celltype1 <- tot2@active.ident

saveRDS(tot2, file = "data/experiment1/robjects/2_filtered_clustered.Rds")

####################
### Add TCR data ###
####################
library(scRepertoire)
rm(list=ls())

tot <- readRDS("data/experiment1/robjects/2_filtered_clustered.Rds")

all <- list()
for (i in c("S21", "S22", "S23", "S24", "S33", "S34")) {
  all[[i]] <- read.csv(paste0("/Volumes/Active/SeqData/public/GSE148837_corridoni_2020/cellranger_v8.0.1/experiment1/", 
                              i, "/outs/per_sample_outs/", i, "/vdj_t/filtered_contig_annotations.csv"))
}
names(all) <- c("s21", "s22", "s23","s24", "s33", "s34")

combined.TCR <- combineTCR(all, 
                           samples = names(all),
                           removeNA = FALSE, 
                           removeMulti = F, 
                           filterMulti = FALSE)

## Add the data to Seurat
sce <- combineExpression(combined.TCR, 
                         tot, 
                         cloneCall="strict", 
                         group.by = "sample",
                         cloneSize = c(Single = 1, Small = 5, 
                                       Medium = 50, Large = 150),
                         proportion = F)

DimPlot(sce, group.by = "cloneSize")

# see if there are pathogenic cells here by our loosest possible definition
sce$pathogenic <- ifelse(grepl("^.{6,7}[Y|F]ST.{5}", sapply(strsplit(sce$CTaa, "_"), "[", 2)), TRUE, FALSE)
sce$pathogenic <- ifelse(grepl("TRAV21", sce$CTstrict), sce$pathogenic, FALSE)

DimPlot(sce, group.by = "pathogenic") # nope, none there

sce$celltype1 <- sce@active.ident

#### Add condition metadata ####
sce$condition <- ifelse(sce$orig.ident %in% c("S21", "S22", "S23"), "HC", "UC")
sce$condition <- factor(sce$condition, levels = c("UC", "HC"))
DimPlot(sce, group.by = "condition")

saveRDS(sce, file = "data/experiment1/robjects/3_filtered_clustered_withTCR.Rds")
