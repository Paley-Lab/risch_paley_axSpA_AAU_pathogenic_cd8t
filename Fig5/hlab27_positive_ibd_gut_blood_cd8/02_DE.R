setwd("/Volumes/Active/Isabel_Risch/MP024/")
source("~/Documents/resources/scripts/universal_helpers/netbid_volcano_edit.R")

library(Seurat)
library(scRepertoire)
library(dplyr)
library(ggplot2)
library(NetBID2)

d <- readRDS("data/robjects/5_final_harmonized.Rds")

########## Each cluster vs. all others ###########
  
# Markers for all the various clusters
de <- FindAllMarkers(d, logfc.threshold = 0, test.use = "MAST")
write.table(de, file="data/de/de_allClusters.txt", sep = "\t", row.names = T)
# de <- read.delim("data/de/de_allClusters.txt")

for (i in unique(de$cluster)) {
  
  de2 <- de[de$cluster %in% i,]
  de2[de2$p_val_adj==0,"p_val_adj"] <- min(de2[de2$p_val_adj!=0,"p_val_adj"])*0.1
  top10 <- head(de2, 10) %>% rownames()
  
  padj <- 0.01
  lfc <- 0.75
  
  draw.volcanoPlot.edit(de2, label_col = "gene", logFC_col = "avg_log2FC", Pv_col = "p_val_adj", 
                        logFC_thre = lfc, Pv_thre = padj, 
                        main = sprintf("%s vs. Others",i), show_label = F,
                        pdf_file = sprintf("plot/de_volcanoplots/de_volcanoplot_%s_vs_others_nolabels.pdf",i))
  draw.volcanoPlot.edit(de2, label_col = "gene", logFC_col = "avg_log2FC", Pv_col = "p_val_adj", 
                        logFC_thre = lfc, Pv_thre = padj, 
                        main = sprintf("%s vs. Others",i), show_label = T,
                        pdf_file = sprintf("plot/de_volcanoplots/de_volcanoplot_%s_vs_others_labeled.pdf",i))
}



########## Gut vs. Blood ###########

# markers for gut vs. blood
de <- FindMarkers(d, logfc.threshold = 0, test.use = "MAST", group.by="tissue_merged", ident.1="Gut", ident.2="Blood", )
de$gene <- rownames(de)
write.table(de, file="data/de/de_gut_vs_blood_all.txt", sep = "\t", row.names = T)
# de <- read.delim("data/de/de_gut_vs_blood_all.txt")

de[de$p_val_adj==0,"p_val_adj"] <- min(de[de$p_val_adj!=0,"p_val_adj"])*0.1
top10 <- head(de, 10) %>% rownames()

padj <- 0.01
lfc <- 1.25

draw.volcanoPlot.edit(de, label_col = "gene", logFC_col = "avg_log2FC", Pv_col = "p_val_adj", 
                      logFC_thre = lfc, Pv_thre = padj, 
                      main = "Gut vs. Blood, All Cells", show_label = T,
                      pdf_file = "plot/de_volcanoplots/de_volcanoplot_gut_vs_blood_top10.pdf",
                      highlight_genes=top10)
draw.volcanoPlot.edit(de, label_col = "gene", logFC_col = "avg_log2FC", Pv_col = "p_val_adj", 
                      logFC_thre = lfc, Pv_thre = padj, 
                      main = "Gut vs. Blood, All Cells", show_label = F,
                      pdf_file = "plot/de_volcanoplots/de_volcanoplot_gut_vs_blood_nolabels.pdf")
draw.volcanoPlot.edit(de, label_col = "gene", logFC_col = "avg_log2FC", Pv_col = "p_val_adj", 
                      logFC_thre = lfc, Pv_thre = padj, 
                      main = "Gut vs. Blood, All Cells", show_label = T,
                      pdf_file = "plot/de_volcanoplots/de_volcanoplot_gut_vs_blood_labeled.pdf")

