# setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP003/")

setwd("/Volumes/Active/Isabel_Risch/MP003")

library(Seurat)
library(VISION)
library(RColorBrewer)
library(ggplot2)

d <- readRDS("data/robjects/9_TCRgenesRemoved_sexRegressed.Rds")

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

paths <- read.gmt("../MP003/data/pathogenic_core_signatures_v2.gmt")
paths2 <- read.gmt("/Volumes/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cgp.v2023.2.Hs.symbols.gmt")
paths2 <- paths2[grep("GOLDRATH", names(paths2))]
paths <- append(paths, paths2)
paths2 <- read.gmt("/Volumes/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt")
paths <- append(paths, paths2)
paths2 <- read.gmt("../MP021/data/ramesh_2014_th17_signatures.gmt")
paths <- append(paths, paths2)
paths2 <- read.gmt("../MP023/IR_cd161_genesets.gmt")
paths <- append(paths, paths2)

source("~/Documents/resources/scripts/universal_helpers/seurat_aucell.R")

seuratAUC(d, paths = paths, outdir = "AUCell/20260504")

#####################################
### AUCell signature highlighting ### 
#####################################

d <- readRDS("data/robjects/9_TCRgenesRemoved_sexRegressed.Rds")
outdir <- "AUCell/20260504"

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

name <- "core_th17"

name <- colnames(d.2@meta.data)[grep(name, colnames(d.2@meta.data), ignore.case = T)]
p <- FeaturePlot(d.2, features = name, order=T) + scale_colour_gradientn(colors=colorRamps::matlab.like(8))+ theme(plot.title = element_text(size=8))
ggplot2::ggsave(filename = paste0(outdir,"/plot/umap_", name, ".pdf"), width=4, height=3, dpi=600)
p <- VlnPlot(d.2, features = name, pt.size = 0) + geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggsave(plot=p, filename = paste0(outdir,"/plot/vlnplot_", name, ".png"), width=8,height=5) 


