# .libPaths(c("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/R_libraries/4.1.0/", .libPaths()))

prefix <- ifelse(dir.exists("/Volumes/"), "/Volumes/", "/storage1/fs1/paleym/")
setwd(paste0(prefix, "/Active/Isabel_Risch/MP024/"))
library(Seurat)
library(dplyr)
library(ggplot2)
# library(qusage)

d <- readRDS("data/robjects/5_final_harmonized.Rds")
outdir <- "AUCell/20260719/"

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
paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
paths2 <- paths2[grep("INTERLEUK", names(paths2))]
paths <- append(paths, paths2)
paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.v2023.2.Hs.symbols.gmt")
paths2 <- paths2[grep("IL17", names(paths2))]
paths <- append(paths, paths2)
paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c5.go.v2023.2.Hs.symbols.gmt")
paths2 <- paths2[grep("INTERLEUKIN_17", names(paths2))]
paths <- append(paths, paths2)
paths2 <- read.gmt("../MP021/data/ramesh_2014_th17_signatures.gmt")
paths <- append(paths, paths2)
paths2 <- read.gmt("../MP023/IR_cd161_genesets.gmt")
paths <- append(paths, paths2)

source("/home/irisch/resources/universal_helpers/seurat_aucell.R")

seuratAUC(d, paths = paths, outdir = outdir)

#####################################
### AUCell signature highlighting ### 
#####################################

d <- readRDS("data/robjects/5_final_harmonized.Rds")
# outdir <- "AUCell/20251211"

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

name <- "Core_Th17_genes_11g"

name <- colnames(d.2@meta.data)[grep(name, colnames(d.2@meta.data), ignore.case = T)]
p <- FeaturePlot(d.2, features = name, ) + scale_colour_gradientn(colors=colorRamps::matlab.like(8))+ theme(plot.title = element_text(size=8))
ggplot2::ggsave(filename = paste0(outdir,"/plot/umap_", name, ".png"), width=4, height=3)
p <- VlnPlot(d.2, features = name, pt.size = 0) + geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggsave(plot=p, filename = paste0(outdir,"/plot/vlnplot_", name, ".png"), width=8,height=5) 


agg <- d@meta.data %>% 
  group_by(subject, celltype1, tissue_merged) %>%
  summarize(Mean = mean(Core_JCI_MP156_c4_vs_c3_noTCR_20g, na.rm=TRUE))

for (j in levels(d@active.ident)) {
  print(j)
  # blood <- d@meta.data[(d$tissue_merged %in% "Blood") & (d$celltypeFinal %in% j), name]
  # gut <- d@meta.data[(d$tissue_merged %in% "Gut") & (d$celltypeFinal %in% j), name]
  
  current <- agg[agg$celltype1 %in% j,] %>% as.data.frame()
  print(ggplot(data = current, aes(x=tissue_merged,y=Mean)) + 
          geom_boxplot() + geom_point() + ggtitle(j) + 
          geom_line(aes(group = subject, colour = "lightblue")))
  
  tmp <- t.test(x = current[current$tissue_merged %in% "Blood", "Mean"], 
                y=current[current$tissue_merged %in% "Gut", "Mean"], paired=T)
  print(tmp)
}
  

saveRDS(d.2, file = "data/robjects/5_final_harmonized.Rds")
