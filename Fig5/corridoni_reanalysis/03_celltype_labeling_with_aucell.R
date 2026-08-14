library(Seurat)
library(ggplot2)
library(AUCell)
library(GSEABase)
setwd("/Volumes/Active/Isabel_Risch/MP004")

source("~/Documents/resources/scripts/universal_helpers/helpers.R")
source("~/Documents/resources/scripts/universal_helpers/seurat_aucell.R")

###############################################################
### Make our signatures using the original Supp. Data ### 
###############################################################

deg <- openxlsx::read.xlsx("/Users/irisch/Downloads/41591_2020_1003_MOESM3_ESM.xlsx", 
                           sheet = "10x Cluster Markers")

celltypes_corridoni <- list()
for (i in unique(deg$Cluster)) {
  degs <- deg[(deg$Cluster %in% i) & (deg$avg_logFC > 0.5) & (deg$FDR < 0.01), ]
  celltypes_corridoni[[i]] <- degs$Gene
}


#####################################
### AUCell signature highlighting ### 
#####################################
d <- readRDS("data/experiment1/robjects/3_filtered_clustered_withTCR.Rds")
outdir <- "AUCell/celltype_labeling"

seuratAUC(d, paths = celltypes_corridoni, outdir = outdir)

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

for(i in 20:ncol(d.2@meta.data)) {
  name <- colnames(d.2@meta.data)[i]
  p <- FeaturePlot(d.2, features = name, order=T) + scale_colour_gradientn(colors=colorRamps::matlab.like(8))+ theme(plot.title = element_text(size=8))
  ggplot2::ggsave(filename = paste0(outdir,"/plot/umap_", name, ".png"), width=5, height=5)
  p <- VlnPlot(d.2, features = name, pt.size = 0) + geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
  ggsave(plot=p, filename = paste0(outdir,"/plot/vlnplot_", name, ".png"), width=8,height=5) 
}



