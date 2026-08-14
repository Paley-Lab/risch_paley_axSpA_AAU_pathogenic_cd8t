# run this on RIS with the following command:
# bsub -Is -M 64GB -R 'rusage[mem=64GB]' -G compute-gfwu -q oncology-interactive -a 'docker(igrisch/seurat:5.2.0_v2)' /bin/bash
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP021")
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

# loading in our data and making sure it's leveled the way we desire
d <- readRDS("data/robjects/4datasets_integrated.Rds")
d$B27 <- factor(as.character(d$B27), levels = c("B27-", "B27+"))
d$tissue_consensus <- factor(as.character(d$tissue_consensus), levels=c("PBMC", "Aqueous","Synovial_Fluid"))
d$disease_merged <- factor(as.character(d$disease_merged), levels=c("HC", "AS", "AU", "PsA"))
d$subject <- factor(as.character(d$subject), levels=c("HC1564", "HC1660", "HC1788", "HC546", "73-4", "76", "81", "74", "SpA-08", "SpA-07", "SpA-01", "SpA-03", "SpA-04", "UV019",
                                                      "UV027", "UV122", "UV180", "PSA1", "PSA2", "PSA3", "PSA4"))

newPathogenicSeq <- read.csv("../pathogenic_motifs", header = F)
d$pathogenic_updated <- ifelse(grepl(newPathogenicSeq[1,2], sapply(strsplit(d$CTaa, "_"), "[", 2)), TRUE, FALSE)
d$pathogenic_updated <- ifelse(grepl("TRAV21", d$CTstrict), d$pathogenic_updated, FALSE)
d$pathogenic_updated <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", d$CTstrict), d$pathogenic_updated, FALSE)


source("../color_palettes.R")

#### ---------- 12/10/25: Plots for the Psoriasis Foundation grant ---------- ####

# 3 rows. First row, pathogenic cells on the UMAP, split by disease
tmp <- subset(d, subset= disease_merged != "HC")
tmp$tissue_consensus <- factor(tmp$tissue_consensus, levels = c("PBMC", "Aqueous", "Synovial_Fluid"))
cells_highlight <- rownames(tmp@meta.data[tmp$pathogenic,])
p1 <- DimPlot(tmp, cells.highlight = cells_highlight, split.by = "disease_merged") + NoAxes()
p2 <- DimPlot(tmp, cells.highlight = cells_highlight, split.by = "tissue_consensus") + NoAxes()
p3 <- FeaturePlot(tmp, features = c("KLRB1", "RORC", "IL23R"), ncol=3, order = T) & NoAxes()

final <- plot_grid(p1, p2, nrow = 2, align = "hv", axis = "bt")
ggsave(final, filename="tmp.pdf", width=8,height=5.8,units="in")
ggsave(p3, filename="tmp2.pdf", width=9,height=3,units="in")


# 1. Create the violin plot without the default dots
p <- VlnPlot(d, split.by = "tissue_consensus", group.by = "disease_merged", features = "IL26", pt.size = 0)
# 2. Overlay a clean jitter layer and color dots by identity
p1 <- p + geom_jitter(mapping = aes(colour = tissuePalette), data = p$data, alpha = 0.5, size = 1)
ggsave(plot=p1, filename="plot/vlnplot_IL26.pdf", width=6,height=4,units="in")

#### ---------- Dotplots of cluster-by-cluster DE genes ---------- ####
per_cluster_markers <- readRDS("data/per_cluster_markers.Rds")
p1 <- DotPlot(d, features = c(per_cluster_markers, "GZMA", "GZMB", "PRF1", "PDCD1", "GNLY", "NKG7"),
              cols = c("lightgrey", "red")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot=p1, filename="plot/cluster_markers.pdf", width=13,height=2.75,units="in")

#### ---------- Various UMAP plots ---------- ####

# Plotting a UMAP of the datasets
p1 <- DimPlot(d,reduction = "umap.harmony",combine = FALSE, label.size = 2, shuffle = T, group.by = "dataset")
ggsave(plot=p1[[1]], filename="plot/umap_dataset.png", width=5,height=4,units="in")


# Plotting the clusters, first all together, then split by dataset/disease/etc
p1 <- DimPlot(d,reduction = "umap.harmony",combine = FALSE, label.size = 2, shuffle = T)
ggsave(plot=p1[[1]], filename="plot/clusters.png", width=6,height=4,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony", split.by = "dataset", combine = FALSE, label.size = 2, shuffle = T)
ggsave(plot=p1[[1]], filename="plot/clusters_by_dataset.pdf", width=22,height=5,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony", split.by = "disease_merged", combine = FALSE, label.size = 2, shuffle = T)
ggsave(plot=p1[[1]], filename="plot/clusters_by_disease.pdf", width=17,height=5,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony",combine = FALSE, label.size = 2, shuffle = T, group.by = "seurat_clusters")
ggsave(plot=p1[[1]], filename="plot/clusters_original_seurat.png", width=4,height=4,units="in", dpi = 600)


# UMAP of pathogenic cells, first all together, then split by B27 status, tissue, dataset, disease
cells <- rownames(d@meta.data[d$pathogenic_updated,])
p1 <- DimPlot(d,reduction = "umap.harmony",cells.highlight = cells, label.size = 2, order = cells)
ggsave(plot=p1, filename="plot/umap_pathogenic.pdf", width=5,height=4,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony",cells.highlight = cells, label.size = 2, order = cells, split.by = "B27") + NoLegend()
ggsave(plot=p1, filename="plot/umap_pathogenic_by_b27.pdf", width=4.75,height=3,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony",cells.highlight = cells, split.by = "tissue_consensus", combine = FALSE, label.size = 2)
ggsave(plot=p1[[1]], filename="plot/pathogenic_by_tissue.pdf", width=10,height=4,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony",cells.highlight = cells, split.by = "dataset", combine = FALSE, label.size = 2)
ggsave(plot=p1[[1]], filename="plot/pathogenic_by_dataset.pdf", width=12,height=4,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony",cells.highlight = cells, split.by = "disease_merged", combine = FALSE, label.size = 2)
ggsave(plot=p1[[1]], filename="plot/pathogenic_by_disease.png", width=14,height=4,units="in")


# UMAP of MAIT cells, selected by TCR
cells <- rownames(d@meta.data[d$mait2,])
p1 <- DimPlot(d,reduction = "umap.harmony",cells.highlight = cells, label.size = 2, order = cells, cols.highlight = "forestgreen")
ggsave(plot=p1, filename="plot/umap_mait.pdf", width=5,height=4,units="in")

# UMAPs of the various tissues represented, split by disease and dataset etc.
p1 <- DimPlot(d,reduction = "umap.harmony",group.by = "tissue_consensus", shuffle = T, cols = tissuePalette, pt.size = 0.1, order = c("Aqueous"))
ggsave(plot=p1[[1]], filename="plot/umap_tissue.pdf", width=6,height=5,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony", split.by = "disease_merged", combine = FALSE, label.size = 2, shuffle = T, group.by = "tissue_consensus", cols= tissuePalette)
ggsave(plot=p1[[1]], filename="plot/tissue_by_disease.pdf", width=17,height=5,units="in")
p1 <- DimPlot(d,reduction = "umap.harmony",group.by = "tissue_consensus", split.by = "dataset", combine = FALSE, label.size = 2, shuffle = T, cols = tissuePalette)
ggsave(plot=p1[[1]], filename="plot/tissue_by_dataset.pdf", width=15,height=4,units="in")

# UMAP of tissue, split by B27 status.
p1 <- DimPlot(d,reduction = "umap.harmony", split.by = "B27", combine = FALSE, label.size = 2, shuffle = T, group.by = "tissue_consensus", cols= tissuePalette)
ggsave(plot=p1[[1]], filename="plot/tissue_by_b27.pdf", width=9,height=5,units="in")
ggsave(plot=p1[[1]], filename="plot/tissue_by_b27.png", width=9,height=5,units="in", dpi=600)


# UMAP of which disease the patients have
p1 <- DimPlot(d,reduction = "umap.harmony",group.by = "disease_merged", split.by = "dataset", 
              combine = FALSE, label.size = 2, shuffle = T, cols = diseasePalette)
ggsave(plot=p1[[1]], filename="plot/disease_by_dataset.pdf", width=15,height=4,units="in")


## Marker plots, various 
p1 <- FeaturePlot(d, features = c("KLRB1", "RORC", "IL23R"), ncol=3, order = T)
ggsave(plot=p1, filename="plot/umap_rorc_1l23r.pdf", width=14,height=4,units="in")

t17_genes <- c("RORC", "IL23R", "IL22", "IL17A", "IL26")
p1 <- VlnPlot(d, features = t17_genes, ncol=4, group.by = "dataset")
ggsave(plot=p1, filename="plot/vlnplot_type17_genes.pdf", width=14,height=4,units="in")
p1 <- FeaturePlot(d, features = t17_genes, ncol=2, order = T)
ggsave(plot=p1, filename="plot/umap_type17_genes.pdf", width=7,height=6,units="in")
for (i in t17_genes) {
  p1 <- FeaturePlot(d, features = i, ncol=1, order = T)
  ggsave(plot=p1, filename=paste0("plot/umap_",i,".png"), width=3.5,height=3,units="in", dpi = 600)
}


p <- VlnPlot(d, features = "Core_JCI_MP156_c4_vs_c3_noTCR_20g", pt.size = 0.2) + geom_boxplot(position=position_dodge(1)) + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggplot2::ggsave(filename = paste0("plot/aucell_vlnplot_byCluster_", "Core_JCI_MP156_c4_vs_c3_noTCR", ".png"), width=7, height=5)

p <- VlnPlot(d, features = "Core_JCI_MP156_c4_vs_c3_noTCR", pt.size = 0.2, group.by = "pathogenic") + geom_boxplot(position=position_dodge(1)) + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggplot2::ggsave(filename = paste0("plot/aucell_vlnplot_byPathogenic_", "Core_JCI_MP156_c4_vs_c3_noTCR", ".png"), width=4, height=5)

p <- VlnPlot(d, features = "Core_JCI_MP156_c4_vs_c3_noTCR_20g", pt.size = 0.2, split.by = "disease_merged") + geom_boxplot(position=position_dodge(1)) + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggplot2::ggsave(filename = paste0("plot/aucell_vlnplot_byDisease_", "Core_JCI_MP156_c4_vs_c3_noTCR", ".png"), width=10, height=5)

p <- VlnPlot(d, features = "Core_JCI_MP156_c4_vs_c3_noTCR_20g", pt.size = 0.2, group.by = "disease_merged") + geom_boxplot(position=position_dodge(1)) + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggplot2::ggsave(filename = paste0("plot/aucell_vlnplot_byDisease_", "Core_JCI_MP156_c4_vs_c3_noTCR", ".png"), width=5, height=5)

p <- VlnPlot(d, features = "Core_JCI_MP156_c4_vs_c3_noTCR_20g", pt.size = 0.2, group.by = "tissue_consensus") + geom_boxplot(position=position_dodge(1)) + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggplot2::ggsave(filename = paste0("plot/aucell_vlnplot_byTissue_", "Core_JCI_MP156_c4_vs_c3_noTCR", ".png"), width=5, height=5)

#### ---------- Various barplots ---------- ####

## Barplots of proportions of MAITs, pathogenic
source("../helper_functions.R")
d <- annotateMaits(d)
d$celltype_17 <- ifelse(d$mait2, "MAIT", ifelse(d$pathogenic, "pathogenic", "other"))
df <- as.data.frame(table(d$celltype1, d$celltype_17))
colnames(df) <- c("Cluster", "Celltype", "Frequency")
df$Celltype <- factor(df$Celltype, levels=c("other", "MAIT", "pathogenic"))
p<-ggplot(df, aes(x=Cluster, y=Frequency, fill=Celltype)) +
  geom_bar(stat="identity", position = "fill")+theme_classic()+
  # labs(title="Title") +
  theme(plot.title = element_text(hjust = 0.5, size=16)) +
  scale_fill_manual(values=c("lightgrey", "forestgreen", "red"))
# geom_text(aes(label=round(`-Log10(p-value)`, digits=1)), vjust=1.6, color="white", size=6)
ggsave(plot=p, filename = "plot/barplot_mait_pathogenic_by_annotated_cluster.pdf", device="pdf", width=4, height=2)

df <- as.data.frame(table(d$subject, d$celltype_17))
colnames(df) <- c("Subject", "Celltype", "Frequency")
df$Celltype <- factor(df$Celltype, levels=c("other", "MAIT", "pathogenic"))
p<-ggplot(df, aes(x=Subject, y=Frequency, fill=Celltype)) +
  geom_bar(stat="identity", position = "fill")+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=16), axis.text.x = element_text(hjust = 1,vjust=0.5, angle = 90)) +
  scale_fill_manual(values=c("lightgrey", "forestgreen", "red"))
# geom_text(aes(label=round(`-Log10(p-value)`, digits=1)), vjust=1.6, color="white", size=6)
ggsave(plot=p, filename = "plot/barplot_mait_pathogenic_by_subject.pdf", device="pdf", width=4, height=2)


## Barplots of location
source("../helper_functions.R")
df <- as.data.frame(table(d$celltype1, d$tissue_consensus))
colnames(df) <- c("Cluster", "Tissue", "Frequency")
df$Tissue <- factor(as.character(df$Tissue), levels=c("PBMC", "Aqueous","Synovial_Fluid"))
p<-ggplot(df, aes(x=Cluster, y=Frequency, fill=Tissue)) +
  geom_bar(stat="identity", position = "fill")+theme_classic()+
  # labs(title="Title") +
  theme(plot.title = element_text(hjust = 0.5, size=16)) +
  scale_fill_manual(values=tissuePalette)
# geom_text(aes(label=round(`-Log10(p-value)`, digits=1)), vjust=1.6, color="white", size=6)
ggsave(plot=p, filename = "plot/barplot_tissue_by_cluster.pdf", device="pdf", width=4, height=2)

# ---------- Plot barplots of cluster by disease/tissue ------------- #
d$DiseaseTissue <- paste0(d$disease_merged, " ", d$tissue_consensus)
plot <- d@meta.data %>% group_by(celltype1, DiseaseTissue) %>% 
  summarise(NumberOfCells = n()) %>% 
  group_by(DiseaseTissue) %>%
  mutate(TotalCellsInCondition = sum(NumberOfCells))
plot$Proportion <- plot$NumberOfCells/plot$TotalCellsInCondition
write.csv(plot, "plot/barplot_cluster_by_disease_and_tissue_inputData.csv", row.names = F)

df <- read.csv("plot/barplot_cluster_by_disease_and_tissue_inputData.csv", 
               stringsAsFactors = FALSE)
dt_levels <- c("HC PBMC", "AS PBMC", "AU PBMC",
               "AS Synovial_Fluid", "AU Aqueous", "PsA Synovial_Fluid")
df$DiseaseTissue <- factor(df$DiseaseTissue, levels = dt_levels)
df$celltype1 <- factor(df$celltype1, levels = c("GZMK_effector_CD8T", 
                                                      "naive_CD8T",
                                                      "MAIT",
                                                      "NKgene_CD8T", 
                                                      "GNLY_effector_CD8T", 
                                                      "pathogenic_CD8T", 
                                                      "cycling"))

p <- ggplot(df, aes(x = DiseaseTissue, y = Proportion, fill = celltype1)) +
  geom_bar(stat="identity", position = "fill") +
  labs(x = NULL, y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = p, filename = "plot/barplot_cluster_by_disease.pdf", device = "pdf", width = 6, height = 4)

#### ---------- Barplot showing how many samples have pathogenic cells ---------- ####

meta <- d@meta.data
tmp <- meta %>% group_by(B27, pathogenic_updated, subject) %>% summarize(n())
duplicates <- tmp$subject[duplicated(tmp$subject)]
remove <- (tmp$subject %in% duplicates) & !tmp$pathogenic_updated
tmp <- tmp[!remove,]
tmp <- tmp[,1:3]
tmp$B27_by_condition <- paste0(tmp$B27, "_", ifelse(grepl("HC", tmp$subject), "HC", "Disease"))
tmp$B27_by_condition <- factor(tmp$B27_by_condition, levels=c("B27+_HC", "B27-_Disease", "B27+_Disease"))
df <- tmp %>% group_by(B27_by_condition, pathogenic_updated) %>% summarize(number = n())
write.csv(df, file = "../MP021/data/samples_with_pathogeniccells_info.csv", row.names = F)


# df$Tissue <- factor(as.character(df$Tissue), levels=c("PBMC", "Aqueous","Synovial_Fluid"))
p<-ggplot(df, aes(x=B27_by_condition, y=number, fill=pathogenic_updated)) +
  geom_bar(stat="identity")+theme_classic()+scale_y_continuous(breaks = c(0,2,4,6,8,10,12),
                                                               limits = c(0,12)) +
  # labs(title="Title") +
  theme(plot.title = element_text(hjust = 0.5, size=16)) +
  scale_fill_manual(values=c("lightgrey","darkred")) + guides(y = guide_axis(minor.ticks = TRUE))
# geom_text(aes(label=round(`-Log10(p-value)`, digits=1)), vjust=1.6, color="white", size=6)
ggsave(plot=p, filename = "plot/barplot_pathogenic_by_b27.pdf", device="pdf", width=3, height=2)





