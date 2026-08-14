# run this on RIS with the following command:
# bsub -Is -M 64GB -R 'rusage[mem=64GB]' -G compute-gfwu -q oncology-interactive -a 'docker(igrisch/seurat:5.2.0_v2)' /bin/bash
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP021")
library(Seurat)
library(ggplot2)
library(dplyr)

d <- readRDS("data/robjects/4datasets_integrated.Rds")

read.gmt <- function (file) 
{
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}

##########################################
### AUCell signature score calculation ### 
##########################################

paths <- read.gmt("../MP003/data/pathogenic_core_signatures.gmt")
paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
paths2 <- paths2[grep("INTERLEUK", names(paths2))]
paths <- append(paths, paths2)
paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.v2023.2.Hs.symbols.gmt")
paths2 <- paths2[grep("IL17", names(paths2))]
paths <- append(paths, paths2)
paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c5.go.v2023.2.Hs.symbols.gmt")
paths2 <- paths2[grep("INTERLEUKIN_17", names(paths2))]
paths <- append(paths, paths2)
paths2 <- read.gmt("data/ramesh_2014_th17_signatures.gmt")
paths <- append(paths, paths2)
paths2 <- read.gmt("../MP023/IR_cd161_genesets.gmt")
paths <- append(paths, paths2)

source("/home/irisch/resources/universal_helpers/seurat_aucell.R")

load("data/counts_and_meta.RData")
seuratAUC(counts, paths = paths, outdir = "AUCell/20260215")

#####################################
### AUCell signature highlighting ### 
#####################################

d <- readRDS("data/robjects/4datasets_integrated.Rds")
outdir <- "AUCell/20260215"

load(paste0(outdir, "/cells_AUC_step2.RData"))

####### Determine the cells with the given gene signatures or active gene sets ######
set.seed(333)
par(mfrow=c(3,3)) 

########## Plot using Seurat ########

allVals <- (cells_AUC@assays@data$AUC)
rownames(allVals) <- stringr::str_replace_all(rownames(allVals), " \\(", "_")
rownames(allVals) <- stringr::str_replace_all(rownames(allVals), "\\)", "")
allVals <- t(allVals)

d.2 <- AddMetaData(d, allVals, colnames(allVals))

name <- "Core_JCI_MP156_c4_vs_c3_noTCR"

name <- colnames(d.2@meta.data)[grep(name, colnames(d.2@meta.data))]
p <- FeaturePlot(d.2, features = name, order=T) + scale_colour_gradientn(colors=colorRamps::matlab.like(8))+ theme(plot.title = element_text(size=8))
ggplot2::ggsave(filename = paste0(outdir,"/plot/umap_", name, ".png"), width=5, height=5)
p <- VlnPlot(d.2, features = name, pt.size = 0) + geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggsave(plot=p, filename = paste0(outdir,"/plot/vlnplot_", name, ".png"), width=8,height=5) 


saveRDS(d.2, file = "data/robjects/4datasets_integrated_withPathogenicSig_update.Rds")


p <- VlnPlot(d.2, features = name, pt.size = 0, group.by = "tissue_consensus") + 
  geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm")) + 
  scale_fill_manual(values = tissuePalette)
ggsave(plot=p, filename = paste0(outdir,"/plot/vlnplot_", name, "_by_tissue.png"), width=4.5,height=5) 

p <- VlnPlot(d.2, features = name, pt.size = 0, group.by = "disease_merged") + 
  geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm")) + 
  scale_fill_manual(values = diseasePalette)
ggsave(plot=p, filename = paste0(outdir,"/plot/vlnplot_", name, "_by_disease.png"), width=4,height=4) 

