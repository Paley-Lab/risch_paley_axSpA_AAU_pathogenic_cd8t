if(dir.exists("/Volumes/")) {
  setwd("/Volumes/Active/Isabel_Risch/MP024")
} else {
  setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP024")
}
library(Seurat)
library(ggplot2)
library(dplyr)

d <- readRDS("data/robjects/5_final_harmonized.Rds")
d$celltype1 <- factor(as.character(d@active.ident), 
                      levels=unique(c("GZMK_CD8T", 
                               "Naive_and_Memory_CD8T" ,
                               "KLRB1_CD8T",
                               "KLR_KIR_CD8T", 
                               "GNLY_CD8T")))

d$pathogenic_old <- d$pathogenic
newPathogenicSeq <- read.csv("../pathogenic_motifs", header = F)
d$pathogenic <- ifelse(grepl(newPathogenicSeq[1,2], sapply(strsplit(d$CTaa, "_"), "[", 2)), TRUE, FALSE)
d$pathogenic <- ifelse(grepl("TRAV21", d$CTstrict), d$pathogenic, FALSE)
d$pathogenic <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", d$CTstrict), d$pathogenic, FALSE)

#### ---------- Dotplots of cluster-by-cluster DE genes ---------- ####
degs <- read.delim("data/de/de_allClusters.txt")
degs <- degs[(degs$p_val_adj < 0.01) & degs$avg_log2FC > 0,]
degs <- degs[order(degs$avg_log2FC, decreasing = T),]
per_cluster_markers <- character(0)
for (i in unique(degs$cluster)) {
  tmp <- degs[degs$cluster==i,"gene"] %>% head(5)
  per_cluster_markers <- c(per_cluster_markers, tmp)
}
p1 <- DotPlot(d, group.by = "celltype1",
              features = unique(c(per_cluster_markers, 
                           "GZMA", "GZMB", "PRF1", "PDCD1", "GNLY", "NKG7")),
              cols = c("lightgrey", "red")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot=p1, filename="plot/cluster_markers.pdf", width=10,height=2.25,units="in")


# UMAPs of type 17 features and IL-26
t17 <- c("IL22", "RORC", "IL17A", "IL23R", "IL26")
for (i in t17) {
  p <- FeaturePlot(d, features = i, order = T)
  ggsave(plot=p, filename = paste0("plot/umap_",i,".pdf"), width=3.2, height=2.9)
}

# UMAP by tissue (non-granular)
DimPlot(d, group.by = "tissue_merged", cols = c("Blood"="darkred", "Gut" = "#5790C1"))
ggsave(filename = "plot/umap_tissue.pdf", width=4, height=2.9)

# UMAP by tissue (granular)
DimPlot(d, group.by = "tissue")
ggsave(filename = "plot/umap_tissue_granular.pdf", width=5.75, height=2.9)

# UMAP by cluster numbers
DimPlot(d, group.by = "seurat_clusters")
ggsave(filename = "plot/umap_clusters.pdf", width=3.5, height=2.9)

# UMAP by my labeled celltypes
DimPlot(d, cols=c("GZMK_CD8T" = "#F9766D", "Naive_and_Memory_CD8T" = "#C49B00", "KLRB1_CD8T" = "#A68AFF",
                  "KLR_KIR_CD8T" = "#00C095", "GNLY_CD8T" = "#00B6ED"))
ggsave(filename = "plot/umap_celltypes.pdf", width=5.75, height=2.9)

# UMAP by pathogenic
cells <- rownames(d@meta.data[d$pathogenic,])
DimPlot(d, cells.highlight = cells, split.by = "IBD")
ggsave(filename = "plot/umap_pathogenicTCR_by_dx.pdf", width=6, height=2.9)

# UMAP by cluster numbers
DimPlot(d, cells.highlight = cells)
ggsave(filename = "plot/umap_pathogenicTCR.pdf", width=4.5, height=2.9)

# UMAP by pathogenic, split by tissue and IBD status
cells <- rownames(d@meta.data[d$pathogenic,])
gut <- subset(d, tissue_merged %in% "Gut")
blood <- subset(d, tissue_merged %in% "Blood")
p1 <- DimPlot(gut, cells.highlight = cells, split.by = "IBD") + NoLegend()
p2 <- DimPlot(blood, cells.highlight = cells, split.by = "IBD") + NoLegend()
p3 <- cowplot::plot_grid(p1, p2, nrow = 2)
ggsave(plot=p3, filename = "plot/umap_pathogenicTCR_by_dx_and_tissue.pdf", width=6, height=6)


# Barplot of celltypes by tissue
tmp <- table(d@active.ident, d$tissue_merged)
tmp <- reshape2::melt(tmp)
tmp$Var1 <- as.character(tmp$Var1)
colnames(tmp) <- c("Celltype", "Tissue", "value")
tmp$Celltype <- factor(tmp$Celltype, 
                       levels=c("GZMK_CD8T", "Naive_and_Memory_CD8T", 
                                "KLR_KIR_CD8T", "KLRB1_CD8T", "GNLY_CD8T"))
tmp <- tmp %>% group_by(Celltype) %>% 
  mutate(`Cells in Celltype` = sum(value)) %>%
  mutate(`Proportion` = value/`Cells in Celltype`)

ggplot(tmp, aes(fill=Tissue, y=value, x=Celltype)) + 
  geom_bar(position="fill", stat="identity") + guides(x=guide_axis(angle=90)) + theme_minimal() +
scale_fill_manual(values=c("Blood"="darkred", "Gut" = "#5790C1")) 
ggsave(filename = "plot/tissue_barplot_byCluster.pdf", width=3, height=3)
write.csv(tmp, file = "plot/tissue_barplot_byCluster.csv")

# Pathogenic signature violinplot in gut vs. blood in each 
d$celltype1 <- factor(as.character(d@active.ident), levels = c("GZMK_CD8T", "Naive_and_Memory_CD8T", "KLR_KIR_CD8T", "KLRB1_CD8T", "GNLY_CD8T"))
p <- VlnPlot(d, features = "Core_JCI_MP156_c4_vs_c3_noTCR_20g", group.by = "celltype1",
             split.by = "tissue_merged", pt.size = 0, cols = c("Blood"="darkred", "Gut" = "#5790C1")) +
  geom_boxplot(position = position_dodge(width=0.9))
ggsave(plot=p, filename = "plot/vlnplot_pathogenicSignature_split_by_tissue.pdf", width=5, height=4)

# average pathogenic score per patient per tissue, only in the KLRB1+ cluster
tmp <- subset(d, idents = "KLRB1_CD8T")
meta <- tmp@meta.data
means <- meta %>%
  group_by(subject, tissue_merged) %>%
  dplyr::summarize(Mean_pathogenic_score = mean(Core_JCI_MP156_c4_vs_c3_noTCR_20g, na.rm=TRUE))
write.csv(means, file = "data/mean_pathogenic_score_per_patient_KLRB1cluster.csv", row.names = F)

df <- read.csv("data/mean_pathogenic_score_per_patient_KLRB1cluster.csv")

df$tissue_merged <- factor(df$tissue_merged, levels = c("Blood", "Gut"))
df$condition <- ifelse(df$subject == "IBD001", "AAU", "IBD")
df$condition <- factor(df$condition, levels = c("IBD", "AAU"))

# ---- plot ---------------------------------------------------------------
ggplot(df, aes(x = tissue_merged, y = Mean_pathogenic_score)) +
  geom_line(aes(group = subject), color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = tissue_merged), size = 3) +
  scale_color_manual(values = c("Blood" = "#8b0000", "Gut" = "#5790c1")) +
  facet_wrap(~condition, scales = "free_x") +
  labs(x = NULL, y = "Mean pathogenic score", color = "Tissue") +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("plot/mean_pathogenic_score_per_patient_KLRB1cluster.png", width = 3, height = 2.5, dpi = 600)
