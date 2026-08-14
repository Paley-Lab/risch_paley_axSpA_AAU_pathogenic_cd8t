library(Seurat)
library(ggplot2)
library(AUCell)
library(GSEABase)
setwd("/Volumes/Active/Isabel_Risch/MP004")

d <- readRDS("data/experiment1/robjects/3_filtered_clustered_withTCR.Rds")

source("~/resources/universal_helpers/helpers.R")
source("/home/irisch/resources/universal_helpers/seurat_aucell.R")
read.gmt <- readgmt

###############################################
### highlight our pathogenic gene signature ###
###############################################

paths <- read.gmt("../MP003/data/pathogenic_core_signatures_v2.gmt") 
paths <- paths[2:3]
# paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
# paths2 <- paths2[grep("INTERLEUK", names(paths2))]
# paths <- append(paths, paths2)
# paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.v2023.2.Hs.symbols.gmt")
# paths2 <- paths2[grep("IL17", names(paths2))]
# paths <- append(paths, paths2)
# paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c5.go.v2023.2.Hs.symbols.gmt")
# paths2 <- paths2[grep("INTERLEUKIN_17", names(paths2))]
# paths <- append(paths, paths2)
paths2 <- read.gmt("../MP021/data/ramesh_2014_th17_signatures.gmt")
paths <- append(paths, paths2)
paths2 <- read.gmt("../MP023/IR_cd161_genesets.gmt")
paths <- append(paths, paths2)

seuratAUC(d, paths = paths, outdir = "AUCell/20260501/")

#####################################
### AUCell signature highlighting ### 
#####################################

d <- readRDS("data/experiment1/robjects/3_filtered_clustered_withTCR.Rds")
outdir <- "AUCell/20260501"

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

# name <- "Core_JCI_MP156_pathogenic_removeTcrGenes_20g"
name <- "Core_Th17"

name <- colnames(d.2@meta.data)[grep(name, colnames(d.2@meta.data))]
p <- FeaturePlot(d.2, features = name, order=T) + scale_colour_gradientn(colors=colorRamps::matlab.like(8))+ theme(plot.title = element_text(size=8))
ggplot2::ggsave(filename = paste0(outdir,"/plot/umap_", name, ".png"), width=5, height=5)
p <- VlnPlot(d.2, features = name, pt.size = 0) + geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggsave(plot=p, filename = paste0(outdir,"/plot/vlnplot_", name, ".png"), width=8,height=5) 

cowplot::plot_grid(p, p1, FeaturePlot(d, features = "IL26"), FeaturePlot(d, features = "IL23R"), DimPlot(d),DimPlot(d, group.by = "condition"), ncol = 3)
ggsave(filename = "tmp.png", width=16,height=10) 

# saveRDS(d.2, file = "data/robjects/4datasets_integrated_withPathogenicSig_update.Rds")


###############################################
### highlight our pathogenic gene signature ###
###############################################

paths <- read.gmt("../MP003/data/pathogenic_core_signatures_v2.gmt") 
paths <- paths[2:3]
# paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
# paths2 <- paths2[grep("INTERLEUK", names(paths2))]
# paths <- append(paths, paths2)
# paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.v2023.2.Hs.symbols.gmt")
# paths2 <- paths2[grep("IL17", names(paths2))]
# paths <- append(paths, paths2)
# paths2 <- read.gmt("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c5.go.v2023.2.Hs.symbols.gmt")
# paths2 <- paths2[grep("INTERLEUKIN_17", names(paths2))]
# paths <- append(paths, paths2)
paths2 <- read.gmt("../MP021/data/ramesh_2014_th17_signatures.gmt")
paths <- append(paths, paths2)
paths2 <- read.gmt("../MP023/IR_cd161_genesets.gmt")
paths <- append(paths, paths2)

seuratAUC(d, paths = paths, outdir = "AUCell/20260501/")


# ####################################
# ########## Fisher's Exact ##########
# ####################################
# 
# # read in the genesets
# gs2gene <- list("custom" = qusage::read.gmt("../MP003/data/pathogenic_core_signatures.gmt")
# )
# 
# 
# ##### Cluster 4 vs Clusters 2 and 3 ####
# 
# # read in the DE results
# de <- FindMarkers(d, ident.2 = c("Trm", "IL26", "Gzmk_effector_1", "Naive_1", 
#                                  "Naive_2", "Gzmk_effector_2", "CD4_CD8_DP", "TYROBPneg_IELs", 
#                                  "FGFBP2_effector", "TYROBPpos_IELs", "Cycling"), ident.1 = "IL26_CCR6", logfc.threshold = 0)
# de_up <- de[de$avg_log2FC > 1 & de$p_val_adj < 0.05,]
# de_up <- rownames(de_up)
# 
# # Use the Jiyang Yu Lab's package NetBID2 
# res_up <- NetBID2::funcEnrich.Fisher(input_list=de_up, bg_list=rownames(de), Pv_thre=0.5,Pv_adj = 'BH', gs2gene = gs2gene, 
#                                      use_gs = "custom", min_gs_size = 5, max_gs_size = 50)
# NetBID2::draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=10,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,Pv_thre = 1,
#                                  # pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
#                                  cluster_gs=TRUE,cluster_gene = TRUE,h=0.95, Pv_col = "Adj_P")
# 
# 
# 
