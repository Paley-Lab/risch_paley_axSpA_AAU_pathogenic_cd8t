# .libPaths(c("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/R_libraries/4.0.3", 
#             "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
# 
# setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP003/")

setwd("/Volumes/Active/Isabel_Risch/MP003")

library(Seurat)
library(dplyr)

d <- readRDS("data/robjects/9_TCRgenesRemoved_sexRegressed.Rds")

##############################################
###### Getting markers for each cluster ######
##############################################

de <- FindAllMarkers(d, test.use = "MAST")

write.csv(de, file="data/de/de_all_clusters_seuratMAST.csv")

markers <- list()

for(i in unique(de$cluster)) {
  tmp <- de[(de$cluster %in% i) & de$p_val_adj < 0.01  & de$avg_log2FC > 0.5,]
  tmp <- tmp[order(tmp$avg_log2FC, decreasing = T),] %>% head(5)
  markers[[i]] <- tmp
}

markers <- bind_rows(markers)

DotPlot(d, features = unique(markers$gene), cols = c("lightgrey", "red")) + 
  theme(axis.text.x=element_text(angle=45, hjust=1), plot.margin = margin(b = 0.5,t = 1.5,r = 0.5, l = 0.5, "cm"))+ 
  scale_y_discrete(labels=c("1", "2", "3", "4", "5")) + ylab("Cluster")
ggsave(filename = "figures/publication_quality/dotplot_cluster_markers.pdf", width=7.75,height=3.25)

############################################################
###### DE for pathogenic vs. non-pathogenic effectors ######
############################################################

# Cluster 3 (pathogenic) vs. 2 (non-pathogenic effector)
de <- FindMarkers(d, ident.1 = "3", ident.2 = c("2"), logfc.threshold = 0, test.use = "MAST")
write.table(de, file="data/de/de_c3_vs_c2.txt", sep = "\t", row.names = T)
de <- read.delim("data/de/de_c3_vs_c2.txt")
de$gene <- rownames(de)
top10 <- head(de, 10) %>% rownames()
source("~/Documents/resources/scripts/universal_helpers/netbid_volcano_edit.R")
draw.volcanoPlot.edit(de, label_col = "gene", logFC_col = "avg_log2FC", Pv_col = "p_val_adj", 
                      logFC_thre = 0.5, Pv_thre = 0.05, 
                      main = "Cluster 3 (Pathogenic) vs. Cluster 2 (Non-Pathogenic Effector)", show_label = T,
                      pdf_file = "figures/de_volcanoplots/de_volcanoplot_cluster3_vs_cluster2.pdf", highlight_genes=top10)
draw.volcanoPlot.edit(de, label_col = "gene", logFC_col = "avg_log2FC", Pv_col = "p_val_adj", 
                      logFC_thre = 0.5, Pv_thre = 0.05, 
                      main = "Cluster 3 (Pathogenic) vs. Cluster 2 (Non-Pathogenic Effector)", show_label = F,
                      pdf_file = "figures/de_volcanoplots/de_volcanoplot_cluster3_vs_cluster2_nolabels.pdf", 
                      highlight_genes=top10)

#### compare the pathogenic genes to JCI Insights paper, Fig. 3C
jci <- c("TRAV21", "TRBV9", "LGALS3", "PRR5", "CEBPD", 
         "AQP3", "LTB", "CCR6", "S100A6", "LST1", "KLRB1", "CTSH", "S100A4", 
         "CHN1", "IL7R", "TNFRSF25", "SLC4A10", "LTK", "ANXA5", "NR1D1", 
         "CD28", "S100A10", "VIM", "JAML", "RUNX2", "ITGAE", "HPGD", "WDR86-AS1", 
         "CAPG")

# intersect pathogenic genes with the JCI paper
de <- read.delim("data/de/de_c3_vs_c2.txt")
de_up <- rownames(de[de$p_val_adj < 0.05 & de$avg_log2FC > 0.5 ,])
core_pathogenic <- intersect(jci, de_up)

# remove TCR genes from the core pathogenic signature
core_pathogenic_tcrRemoved <- core_pathogenic[!grepl("^TRA|^TRB",core_pathogenic )]

core <- list("Core_JCI_MP156_pathogenic"=core_pathogenic, 
             "Core_JCI_MP156_pathogenic_removeTcrGenes"=core_pathogenic_tcrRemoved)
source("~/Documents/resources/scripts/universal_helpers/helpers.R")
gmtlist2file(core, filename = "data/pathogenic_core_signatures_v2.gmt")
