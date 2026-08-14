.libPaths(c("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/R_libraries/4.0.3", 
            "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP003/")

library(Seurat)

######################################################
### Reclustering the data without TCR genes or sex ### July 19, 2024
######################################################

d <- readRDS("data/robjects/singletsNoBystanders_clustered_withTCR.Rds")
d <- NormalizeData(d) %>% 
  FindVariableFeatures() %>% 
  Trex::quietTCRgenes() %>% ## use this Trex function to mask the TCR genes for clustering purposes
  ScaleData(features = rownames(d), vars.to.regress = "sex") %>% 
  RunPCA(verbose = FALSE)

ElbowPlot(d) #honestly looks like we only need 10 PCs at most
d <- FindNeighbors(d, dims = 1:10)
d <- FindClusters(d, resolution = 0.2, verbose = FALSE)
d <- RunUMAP(d, dims = 1:10)
DimPlot(d, label = TRUE)
DimPlot(d, group.by = "orig.ident")
DimPlot(d, group.by = "HTO_dmm")
DimPlot(d, group.by = "condition")
DimPlot(d, group.by = "pathogenic")

saveRDS(d, file = "data/robjects/9_TCRgenesRemoved_sexRegressed.Rds")

trav21_cells <- d@meta.data[grep("TRAV21", d$CTgene),]

DimPlot(d, cells.highlight = rownames(trav21_cells))



