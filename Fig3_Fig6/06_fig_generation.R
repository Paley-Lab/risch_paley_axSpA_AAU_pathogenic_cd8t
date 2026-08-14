#### Figure generation for various things

setwd("/Volumes/Active/Isabel_Risch/MP003")

library(Seurat)
library(scRepertoire)
library(dplyr)
library(ggplot2)
library(NetBID2)

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


d <- readRDS("data/robjects/9_TCRgenesRemoved_sexRegressed.Rds")
newPathogenicSeq <- read.csv("../pathogenic_motifs", header = F)
d$pathogenic_updated <- ifelse(grepl(newPathogenicSeq[1,2], sapply(strsplit(d$CTaa, "_"), "[", 2)), TRUE, FALSE)
d$pathogenic_updated <- ifelse(grepl("TRAV21", d$CTstrict), d$pathogenic_updated, FALSE)
d$pathogenic_updated <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", d$CTstrict), d$pathogenic_updated, FALSE)


source("../helper_functions.R")

# - UMAP annotated / labelled (Memory (1&2), Cytotoxic (3), KLRB1+ (4), Transitional? (2), Other (5))
# new.cluster.ids <- c("Memory_like","Memory_like","Effector_Transitional", "Effector_Cytotoxic", "Effector_KLRB1", "SIK3_NFKB1_hi")
# names(new.cluster.ids) <- levels(d)
# d <- RenameIdents(d, new.cluster.ids)
# d$celltype <- d@active.ident
DimPlot(d, reduction = "umap", label = F, pt.size = 0.5) 
ggsave(filename = "figures/publication_quality/celltypes_umap.pdf", width=4, height=2.9)
ggsave(filename = "figures/publication_quality/celltypes_umap.png", width=4, height=2.9, dpi = 600)


# ## CHANGE THESE COLORS
# DimPlot(d, reduction = "umap", label = F, pt.size = 0.5, cols = c("#C49B00", "palegreen3", "#F9766D", "#A68AFF", "salmon4")) 
# ggsave(filename = "figures/publication_quality/celltypes_umap_updatedColors.pdf", width=4, height=2.9)


# c("GZMK_CD8T" = "#F9766D", "Naive_CD8T" = "#C49B00", "MAIT" = "#53B400", "KLRB1_CD8T" = "#A68AFF",
#   "KIR_CD8T" = "#00C095", "GNLY_CD8T" = "#00B6ED", "Cycling" = "#FC61D9")


# - UMAP of pathogenic cells highlighted red (from patients) and blue (from healthy) on grey cell background
# --- do the pathogenic TCRs in the healthy individual cluster elsewhere from those in the patients?
d$tmp <- ifelse(d$pathogenic_updated, ifelse(d$condition %in% "HC", "healthy_pathogenic", "patient_pathogenic"), "non-pathogenic")
d$tmp <- factor(d$tmp, levels=c("non-pathogenic", "healthy_pathogenic", "patient_pathogenic"))
p <- DimPlot(d, group.by = "tmp", cols = c("lightgrey", "#eb3b93", "#189d75"), order = c("healthy_pathogenic", "patient_pathogenic", "non-pathogenic"), pt.size = 0.5) + ggtitle("Pathogenic Cells")
ggsave(plot=p, filename = "figures/publication_quality/pathogenic_umap.pdf", width=5.3, height=3)
ggsave(plot=p, filename = "figures/publication_quality/pathogenic_umap.png", width=5.3, height=3, dpi=600, device="png")


# - UMAP of genotype
DimPlot(d, group.by = "condition", cols = c("hotpink", "#189d75")) + ggtitle("Disease Status")
ggsave(filename = "figures/publication_quality/condition_umap.pdf", width=5, height=3)
ggsave(filename = "figures/publication_quality/condition_umap.png", width=5, height=2.9, dpi=600, device="png")


# FeaturePlots of interesting genes
features <- c("GNLY", "PRF1", "GZMA", "GZMB", "CCR6", "KLRB1", "GZMK", "LEF1", #"IL17A",
              "IL23R", "RORC", "IL22")
for (i in features) {
  FeaturePlot(d, features = i, order = T)
  ggsave(filename = paste0("figures/publication_quality/expression_",i,"_umap.pdf"), width=4.2, height=3.5)
}

# FeaturePlots of interesting genes (Dotplot version)
DotPlot(d, features = c("TCF7", "SELL", "LEF1","CCR7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMK", "IFNG"), 
        cols = c("lightgrey", "red")) + 
  theme(axis.text.x=element_text(angle=45, hjust=1), plot.margin = margin(2,1,1,2, "cm")) + 
  scale_y_discrete(labels=c("1", "2", "3", "4", "5")) + ylab("Cluster")
ggsave(filename = "figures/publication_quality/dotplot_markers_supplementary.pdf", width=6,height=3)

# UMAP of clone size in each cluster
DimPlot(d, group.by = "cloneSize", cols = c('#FDE725FF', '#3CBB75FF','#2d708EFF','#440154FF'  ))
ggsave(filename = "figures/publication_quality/umap_cloneSize.pdf", width=6, height=3)

# --------- Bar plots of clone size --------- #
my_sum <- d@meta.data %>%
  group_by(cloneSize, seurat_clusters) %>%
  summarise( 
    n=n()
  ) %>% group_by(seurat_clusters) %>%
  mutate(
    total_in_cluster=sum(n)
  )
my_sum$percent_per_cluster <- (my_sum$n/my_sum$total_in_cluster)*100
my_sum$cloneSize <- as.character(my_sum$cloneSize)
my_sum[is.na(my_sum$cloneSize), "cloneSize"] <- "NA"
my_sum$cloneSize <- factor(my_sum$cloneSize, levels=c("Large (50 < X <= 150)", "Medium (5 < X <= 50)", "Small (1 < X <= 5)", "Single (0 < X <= 1)",
                                                         "NA"
                                                          ))

# adjust cluster numbering to match 1-indexing in paper
my_sum$seurat_clusters <- as.numeric(my_sum$seurat_clusters) %>%
  factor(levels=c("1","2","3","4","5"))
ggplot(my_sum, aes(fill=cloneSize, y=percent_per_cluster, x=seurat_clusters)) + 
  geom_bar(position="fill", stat="identity") + guides(x=guide_axis(angle=90)) + theme_minimal() +
  scale_fill_manual(values = c('#FDE725FF', '#3CBB75FF','#2d708EFF','#440154FF','lightgrey'))
ggsave(filename = "figures/barplot_proportion_cloneSize_per_cluster.pdf", width=5, height=3)
ggplot(my_sum, aes(fill=cloneSize, y=n, x=seurat_clusters)) + 
  geom_bar( stat="identity") + guides(x=guide_axis(angle=90)) + theme_minimal() +
  scale_fill_manual(values = c('#FDE725FF', '#3CBB75FF','#2d708EFF','#440154FF','lightgrey'))
ggsave(filename = "figures/barplot_cellnumber_cloneSize_per_cluster.pdf", width=5, height=3)

my_sum$percent_per_cluster <- signif(my_sum$percent_per_cluster, digits = 3)
write.csv(my_sum, file="figures/publication_quality/suppFig3d_barplot_cloneSize_per_cluster.csv", row.names = F)


# -------- Fraction of pathogenic clones per individual -------- #

d$tmp <- ifelse(d$pathogenic_updated, "pathogenic", "non-pathogenic")
tmp <- table(d$tmp, d$HTO_dmm)
tmp <- reshape2::melt(tmp)
colnames(tmp) <- c("Pathogenic", "Sample", "value")
tmp$Sample <- factor(as.character(tmp$Sample), levels = c("HC284-002-TotalSeqC", "HC352-001-2-TotalSeqC", 
                                                          "HC357-001-TotalSeqC", "HC366-001-TotalSeqC", "HC387-001-TotalSeqC",
                                                          "003-001-TotalSeqC", "026-001-TotalSeqC", "030-001-TotalSeqC", 
                                                          "049-001-TotalSeqC"))
tmp$total <- 0
for (i in unique(tmp$Sample)) {
  sumTib <- tmp[tmp$Sample %in% i,]
  sumTib <- sum(sumTib$value)
  tmp[tmp$Sample %in% i,"total"] <- sumTib
}
tmp$pct <- (tmp$value/tmp$total)*100
tmp <- tmp[!grepl("non", tmp$Pathogenic),]
tmp$condition <- ifelse(grepl("HC", tmp$Sample), "Healthy", "Patient")
write.table(tmp, file="figures/publication_quality/pathogenic_percent_bySample.txt", sep="\t", row.names=F, col.names=T)

p<-ggplot(tmp, aes(x=condition, y=pct)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center') + scale_fill_manual(values=c("hotpink", "#189d75"))
p 

# -------- Number and percent of sample in each cluster, barplots -------- #

cluster_breakdown <- d@meta.data %>% 
  group_by(HTO_dmm, seurat_clusters, condition) %>%
  summarise(number_of_cells = n()) %>% 
  group_by(HTO_dmm) %>%
  mutate(total_cells_in_sample = sum(number_of_cells))
cluster_breakdown$percentage_of_sample_in_cluster <- cluster_breakdown$number_of_cells /
  cluster_breakdown$total_cells_in_sample
cluster_breakdown$percentage_of_sample_in_cluster <- signif(cluster_breakdown$percentage_of_sample_in_cluster, digits = 3) 
cluster_breakdown$seurat_clusters <- as.numeric(cluster_breakdown$seurat_clusters)

write.csv(cluster_breakdown, file = "figures/publication_quality/suppFig3a_sample_breakdown_by_cluster.csv", 
          row.names = F)

## =============================================================
## PRISM-style barplots: cell number & % of sample per cluster
## Healthy Control (HC) vs Disease
## =============================================================


## -------------------- 1. Load data ---------------------------
df <- read.csv("figures/publication_quality/suppFig3a_sample_breakdown_by_cluster.csv", stringsAsFactors = FALSE)

df$seurat_clusters <- factor(df$seurat_clusters)
df$condition <- factor(df$condition, levels = c("HC", "Disease"))
df$percentage_of_sample_in_cluster <- df$percentage_of_sample_in_cluster * 100  # convert to %

## -------------------- 2. Summary stats (mean +/- SEM) ---------
summarise_metric <- function(data, metric) {
  data %>%
    group_by(seurat_clusters, condition) %>%
    summarise(
      mean_val = mean(.data[[metric]]),
      sem_val  = sd(.data[[metric]]) / sqrt(n()),
      sd_val  = sd(.data[[metric]]),
      n        = n(),
      .groups = "drop"
    )
}

summary_cellnum <- summarise_metric(df, "number_of_cells")
summary_pct     <- summarise_metric(df, "percentage_of_sample_in_cluster")

## -------------------- 3. PRISM-style theme --------------------
theme_prism <- function() {
  theme_classic(base_size = 14) +
    theme(
      axis.line       = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.8, colour = "black"),
      axis.text       = element_text(colour = "black", face = "bold"),
      axis.title      = element_text(colour = "black", face = "bold"),
      legend.title    = element_blank(),
      legend.position = "top",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
}

fill_colors <- c("HC" = "#2CA02C", "Disease" = "#D620A0")  # green / magenta

## -------------------- 4. Plot 1: Cell number per cluster -------
p_cellnum <- ggplot(summary_cellnum, aes(x = seurat_clusters, y = mean_val, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, colour = "black", linewidth = 0.6) +
  geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.6) +
  geom_point(data = df, aes(x = seurat_clusters, y = number_of_cells, fill = condition),
             position = position_dodge(width = 0.8), shape = 21, size = 2,
             colour = "black", alpha = 0.7) +
  scale_fill_manual(values = fill_colors) +
  coord_cartesian(ylim = c(0, 500)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cluster", y = "Cell number per cluster", title = "Cell Number by Cluster") +
  theme_prism()

## -------------------- 5. Plot 2: Percentage of sample per cluster ----
p_pct <- ggplot(summary_pct, aes(x = seurat_clusters, y = mean_val, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, colour = "black", linewidth = 0.6) +
  geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.6) +
  geom_point(data = df, aes(x = seurat_clusters, y = percentage_of_sample_in_cluster, fill = condition),
             position = position_dodge(width = 0.8), shape = 21, size = 2,
             colour = "black", alpha = 0.7) +
  scale_fill_manual(values = fill_colors) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cluster", y = "% of sample in cluster", title = "Percentage of Sample by Cluster") +
  theme_prism()

## -------------------- 6. Save plots ----------------------------
ggsave("cluster_cell_number_barplot.pdf", p_cellnum, width = 6, height = 5)
ggsave("cluster_percentage_barplot.pdf", p_pct, width = 6, height = 5)

## =============================================================
## 7. Mann-Whitney U tests: Disease vs HC, per cluster
## =============================================================
run_mwu_per_cluster <- function(data, metric) {
  clusters <- levels(data$seurat_clusters)
  results <- lapply(clusters, function(cl) {
    sub <- data[data$seurat_clusters == cl, ]
    hc_vals      <- sub[[metric]][sub$condition == "HC"]
    disease_vals <- sub[[metric]][sub$condition == "Disease"]
    
    test <- tryCatch(
      wilcox.test(disease_vals, hc_vals, exact = FALSE),
      error = function(e) NULL
    )
    
    data.frame(
      cluster  = cl,
      metric   = metric,
      n_HC     = length(hc_vals),
      n_Disease = length(disease_vals),
      W_statistic = if (!is.null(test)) unname(test$statistic) else NA,
      p_value  = if (!is.null(test)) test$p.value else NA
    )
  })
  do.call(rbind, results)
}

mwu_cellnum <- run_mwu_per_cluster(df, "number_of_cells")
mwu_pct     <- run_mwu_per_cluster(df, "percentage_of_sample_in_cluster")

## Combine and apply multiple-testing correction (BH/FDR) within each metric
mwu_results <- rbind(mwu_cellnum, mwu_pct)
mwu_results$p_adj_BH <- ave(mwu_results$p_value, mwu_results$metric,
                            FUN = function(p) p.adjust(p, method = "BH"))

print(mwu_results)
write.csv(mwu_results, "data/mannwhitney_disease_vs_HC_per_cluster.csv", row.names = FALSE)

######################################################################################
### Fisher's Exact Test to compare path. signature to type 17, CD161int, and MAIT  ### March 9, 2026
######################################################################################

# get a more representative list of background genes so as not to artificially inflate the P-values
bg <- rownames(d)

# pull our pathogenic signature
pathogenic <- read.gmt("../MP003/data/pathogenic_core_signatures_v2.gmt")

# test against our curated Type 17 and CD161 signatures, and Canonical Pathways signatures 
paths <- read.gmt("../MP021/data/ramesh_2014_th17_signatures.gmt")
paths2 <- read.gmt("../MP023/IR_cd161_genesets_noRedundancies.gmt")
paths <- append(paths, paths2)
paths <- append(paths, list("background"= bg) )
# paths2 <- read.gmt("../MP023/paper_cd161_genesets.gmt")
# paths <- append(paths, paths2)
paths2 <- read.gmt("/Volumes/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c2.cp.v2023.2.Hs.symbols.gmt")
paths <- append(paths, paths2)
paths2 <- read.gmt("../MP021/data/4publicdataset_signatures.gmt")
paths <- append(paths, paths2)
fisher <- NetBID2::funcEnrich.Fisher(input_list = pathogenic[["Core_JCI_MP156_pathogenic_removeTcrGenes"]],
                                     bg_list = bg, 
                                     Pv_thre = 1, gs2gene = paths, Pv_adj = "BH")
NetBID2::draw.funcEnrich.cluster(fisher, Pv_thre = 0.05, Pv_col = "Adj_P")
write.table(fisher, file = 'FisherExactTest/results_type17_signatures_curated.tsv', sep = "\t")

# Test Gene Ontology signatures
paths<- read.gmt("/Volumes/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/c5.go.v2023.2.Hs.symbols.gmt")
paths <- append(paths, list("background"= bg) )
fisher <- NetBID2::funcEnrich.Fisher(input_list = pathogenic[["Core_JCI_MP156_pathogenic_removeTcrGenes"]], bg_list = bg, 
                                     Pv_thre = 0.05, gs2gene = paths, Pv_adj = "BH")
# NetBID2::draw.funcEnrich.cluster(fisher)
# NO RESULTS FOUND 

# Test Hallmark signatures
paths<- read.gmt("/Volumes/mgriffit/Active/griffithlab/adhoc/irisch/msigdb/genesets/h.all.v2023.2.Hs.symbols.gmt")
paths <- append(paths, list("background"= bg) )
fisher <- NetBID2::funcEnrich.Fisher(input_list = pathogenic[["Core_JCI_MP156_pathogenic_removeTcrGenes"]], bg_list = bg, 
                                     Pv_thre = 0.05, gs2gene = paths, Pv_adj = "BH")
# NetBID2::draw.funcEnrich.cluster(fisher)
# NO RESULTS FOUND 







