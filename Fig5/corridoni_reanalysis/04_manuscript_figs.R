setwd("/Volumes/Active/Isabel_Risch/MP004")
library(Seurat)
library(ggplot2)
library(dplyr)

d <- readRDS("data/experiment1/robjects/3_filtered_clustered_withTCR.Rds")


d$celltype_merged <- ifelse(grepl("aive", as.character(d@active.ident)), "Naive", as.character(d@active.ident))
d$celltype_merged <- ifelse(grepl("TYROBPneg_IELs", d$celltype_merged), "TYROBPneg_IELs", d$celltype_merged)
d$celltype_merged <- ifelse(grepl("GZMK_effectors", d$celltype_merged), "GZMK_effectors", d$celltype_merged)
d$celltype_merged <- factor(d$celltype_merged, levels = c("MAIT", "CD8_CD4", "IL26", "GZMK_effectors", "TYROBPneg_IELs", 
                                                              "TYROBPpos_IELs", "FGFBP2", "CD8_CD4_FOXP3", "Trm","Naive"))

DimPlot(d, group.by = "celltype_merged")
ggsave(filename = "plot/clustering.pdf", width = 6,height = 3.5)

DimPlot(d, group.by = "condition", cols= c("HC" = "#7FB37A", "UC" = "#F1A83D"))
ggsave(filename = "plot/manuscriptFigs/umap_condition.pdf", width = 4.75,height = 3.5)

########################################
######## AUCell, Experiment 1 ##########
########################################
library(AUCell)
load("AUCell/20260501/cells_AUC_step2.RData")
load("AUCell/20260501/cells_rankings.RData")

####### Determine the cells with the given gene signatures or active gene sets ######
set.seed(333)
par(mfrow=c(3,3)) 
# cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

########## Plot using Seurat ########

allVals <- (cells_AUC@assays@data$AUC)
rownames(allVals) <- stringr::str_replace_all(rownames(allVals), " \\(", "_")
rownames(allVals) <- stringr::str_replace_all(rownames(allVals), "\\)", "")
allVals <- t(allVals)


d.2 <- AddMetaData(d, allVals, colnames(allVals))

name <- "Core_JCI_MP156_pathogenic_removeTcrGenes_20g"

name <- colnames(d.2@meta.data)[grep(name, colnames(d.2@meta.data))]
p <- FeaturePlot(d.2, features = name, pt.size = 0.7, order=T) + scale_colour_gradientn(colors=colorRamps::matlab.like(8))+ theme(plot.title = element_text(size=8))
ggplot2::ggsave(filename = paste0("plot/manuscriptFigs/umap_", name, ".png"), width=6, height=5, dpi=600)
p <- VlnPlot(d.2, features = name, pt.size = 0, group.by="celltype_merged") + geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggsave(plot=p, filename = paste0("plot/manuscriptFigs/vlnplot_", name, ".pdf"), width=8,height=4) 
p <- VlnPlot(d.2, features = name, pt.size = 0, group.by = "condition") + geom_boxplot() + theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
ggsave(plot=p, filename = paste0("plot/manuscriptFigs/vlnplot_byCondition_", name, ".png"), width=8,height=4) 

########## Barplot of celltypes by disease ######## 

# export a table 
celltypes_by_condition <- d@meta.data %>% group_by(orig.ident, celltype_merged, condition) %>%
  summarise(
    count = n()
  )
write.csv(celltypes_by_condition, file = "data/celltypes_by_condition.csv", row.names = F)


# --- Load data ---
df <- read.csv("data/celltypes_by_condition.csv")
df$celltype_merged <- factor(df$celltype_merged, levels = c("MAIT", "CD8_CD4", "IL26", "GZMK_effectors", "TYROBPneg_IELs", 
                                     "TYROBPpos_IELs", "FGFBP2", "CD8_CD4_FOXP3", "Trm","Naive"))


# --- Summarize: mean and SEM per celltype x condition ---
summary_df <- df %>%
  group_by(celltype_merged, condition) %>%
  summarise(
    mean_count = mean(count),
    sem = sd(count) / sqrt(n()),
    .groups = "drop"
  )

# --- Plot ---
ggplot(summary_df, aes(x = celltype_merged, y = mean_count, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = mean_count - sem, ymax = mean_count + sem),
    position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.4
  ) +
  geom_point(
    data = df, aes(x = celltype_merged, y = count, fill = condition),
    position = position_dodge(width = 0.8), size = 1.5, shape = 21, color = "black"
  ) +
  scale_fill_manual(values = c("HC" = "#7FB37A", "UC" = "#F1A83D")) +
  labs(x = NULL, y = "Average cell number", fill = "Condition") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "top"
  )


# --- Compute per-sample frequency (% of total cells in that sample) ---
freq_df <- df %>%
  group_by(orig.ident) %>%
  mutate(freq = 100 * count / sum(count)) %>%
  ungroup()

# --- Summarize: mean and SEM frequency per celltype x condition ---
summary_freq <- freq_df %>%
  group_by(celltype_merged, condition) %>%
  summarise(
    mean_freq = mean(freq),
    sem = sd(freq) / sqrt(n()),
    .groups = "drop"
  )

# --- Plot ---
ggplot(summary_freq, aes(x = celltype_merged, y = mean_freq, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = mean_freq - sem, ymax = mean_freq + sem),
    position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.4
  ) +
  geom_point(
    data = freq_df, aes(x = celltype_merged, y = freq, fill = condition),
    position = position_dodge(width = 0.8), size = 1.5, shape = 21, color = "black"
  ) +
  scale_fill_manual(values = c("HC" = "#7FB37A", "UC" = "#F1A83D")) +
  labs(x = NULL, y = "% of sample", fill = "Condition") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "top"
  )
ggsave(filename = "plot/barplot_celltype_frequency.pdf", width = 3.5, height=4)


