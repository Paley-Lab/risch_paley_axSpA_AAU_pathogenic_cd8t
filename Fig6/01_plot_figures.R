library(Seurat)
library(ggplot2)
library(dplyr)

## Make a little helper function for the AUCell highlighting
aucell_highlight <- function(d, # Seurat object that we're currently using
                             outdir, # directory where the AUCell results live
                             signature_name = "Core_Th17" # name of the signature to highlight
                             ) {
  load(paste0(outdir, "/cells_AUC_step2.RData"))
  set.seed(333)
  allVals <- (cells_AUC@assays@data$AUC)
  rownames(allVals) <- stringr::str_replace_all(rownames(allVals), " \\(", "_")
  rownames(allVals) <- stringr::str_replace_all(rownames(allVals), "\\)", "")
  allVals <- t(allVals)
  d.2 <- AddMetaData(d, allVals, colnames(allVals))
  
  name <- signature_name
  name <- colnames(d.2@meta.data)[grep(name, colnames(d.2@meta.data))]
  p <- FeaturePlot(d.2, features = name, order=T) + scale_colour_gradientn(colors=colorRamps::matlab.like(8))+ theme(plot.title = element_text(size=8))
  ggplot2::ggsave(filename = paste0(outdir,"/plot/umap_", name, ".png"), width=5, height=5)
  p <- VlnPlot(d.2, features = name, pt.size = 0) + geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
  ggsave(plot=p, filename = paste0(outdir,"/plot/vlnplot_", name, ".png"), width=8,height=5) 
  
}

##### ------ YeiH+CD8+ T Cells from HLA-B27+ Patient/HC PBMCs ------ #####
d <- readRDS("/Volumes/Active/Isabel_Risch/MP003/data/robjects/9_TCRgenesRemoved_sexRegressed.Rds")
# FeaturePlots of interesting genes
features <- c( "IL17A", "IL23R", "RORC", "IL22")
for (i in features) {
  FeaturePlot(d, features = i, order = T)
  ggsave(filename = paste0("figures/publication_quality/expression_",i,"_umap.pdf"), 
         width=4.2, height=3.5)
}

### AUCell signature highlighting 
outdir <- "/Volumes/Active/Isabel_Risch/MP003/AUCell/20260504/"
aucell_highlight(d, outdir)

##### ------ 4 public datasets ------ #####
# loading in our data and making sure it's leveled the way we desire
d <- readRDS("/Volumes/Active/Isabel_Risch/MP021/data/robjects/4datasets_integrated.Rds")
t17_genes <- c("RORC", "IL23R", "IL22", "IL17A", "IL26")
p1 <- FeaturePlot(d, features = t17_genes, ncol=2, order = T)
ggsave(plot=p1, filename="plot/umap_type17_genes.pdf", width=7,height=6,units="in")
for (i in t17_genes) {
  p1 <- FeaturePlot(d, features = i, ncol=1, order = T)
  ggsave(plot=p1, filename=paste0("plot/umap_",i,".png"), width=3.5,height=3,units="in", dpi = 600)
}

### AUCell signature highlighting 
outdir <- "/Volumes/Active/Isabel_Risch/MP021/AUCell/20260215/"
aucell_highlight(d, outdir)

##### ------ HLA-B27+ IBD gut/blood dataset ------ #####
d <- readRDS("/Volumes/Active/Isabel_Risch/MP024/data/robjects/5_final_harmonized.Rds")
d$celltype1 <- factor(as.character(d@active.ident), 
                      levels=unique(c("GZMK_CD8T", 
                                      "Naive_and_Memory_CD8T" ,
                                      "KLRB1_CD8T",
                                      "KLR_KIR_CD8T", 
                                      "GNLY_CD8T")))
# UMAPs of type 17 features and IL-26
t17 <- c("IL22", "RORC", "IL17A", "IL23R", "IL26")
for (i in t17) {
  p <- FeaturePlot(d, features = i, order = T)
  ggsave(plot=p, filename = paste0("plot/umap_",i,".pdf"), width=3.2, height=2.9)
}

### AUCell signature highlighting 
outdir <- "/Volumes/Active/Isabel_Risch/MP024/AUCell/20260719/"
aucell_highlight(d, outdir)