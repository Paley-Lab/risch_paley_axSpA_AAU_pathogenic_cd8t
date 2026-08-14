# run this on RIS with the following command:
# bsub -Is -M 128GB -R 'rusage[mem=128GB]' -G compute-mgriffit -q oncology-interactive -a 'docker(igrisch/seurat:5.2.0_v2)' /bin/bash
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP021")
library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(openxlsx)

outdir <- "data/de/DESeq2/"
smallestGroupSize = 3
minPctExpressed = 0.1
load("/storage1/fs1/paleym/Active/Isabel_Risch/MP021/data/counts_and_meta.RData")
filtered <- CreateSeuratObject(counts = counts, meta.data = meta, min.cells = minPctExpressed*dim(counts)[2])

# create a pseudobulked matrix, splitting by patient/tissue and adding disease for ease of visualiz.
filtered$cluster_patient_tissue <- paste0(filtered$celltype1, "_", filtered$subject,"_", filtered$tissue_consensus)
counts <- Seurat::AggregateExpression(filtered, assays = "RNA", group.by = "cluster_patient_tissue", slot = "counts")[[1]]

# create a simple metadata dataframe
meta <- meta[,c("celltype1", "subject", "dataset", "disease_merged", "tissue_consensus")]
colnames(meta) <- c("cluster", "subject", "dataset", "dx", "tissue")
meta <- unique(meta)
rownames(meta) <- paste0(meta$cluster, "_", meta$subject,"_", meta$tissue)
meta <- meta[colnames(counts),]

# make HC the reference level for diagnosis
meta$dx <- as.factor(meta$dx)
meta$dx <- relevel(meta$dx, "HC")

# make PBMC the reference level for tissue
meta$tissue <- as.factor(meta$tissue)
meta$tissue <- relevel(meta$tissue, "PBMC")


######## Set up our object

dds <- DESeqDataSetFromMatrix(counts,
                              colData = meta,
                              design = as.formula("~0 + subject + tissue + cluster"))

if (!dir.exists(outdir)) {
  dir.create(paste0(outdir, "/QC"), recursive = T)
}

# perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)

## QUALITY CONTROL, SAMPLE LEVEL ##

# Transform counts for data visualization
# rld <- rlog(dds) 

vst <- DESeq2::varianceStabilizingTransformation(dds, fitType = "local") # try vst instead

if (!dir.exists(paste0(outdir, "/QC/"))) dir.create(paste0(outdir, "/QC/"), recursive=T)

# Plot PCA
# png(paste0(outdir, "/QC/rld_disease_pca.png"), width = 8, height=5, res = 300, units = "in")
# print(DESeq2::plotPCA(rld, ntop = 500, intgroup = "dx"))
# dev.off()
png(paste0(outdir, "/QC/vst_disease_pca.png"), width = 8, height=5, res = 300, units = "in")
print(DESeq2::plotPCA(vst, ntop = 500, intgroup = "dx"))
dev.off()

# png(paste0(outdir, "/QC/rld_tissue_pca.png"), width = 8, height=5, res = 300, units = "in")
# print(DESeq2::plotPCA(rld, ntop = 500, intgroup = "tissue"))
# dev.off()
png(paste0(outdir, "/QC/vst_tissue_pca.png"), width = 8, height=5, res = 300, units = "in")
print(DESeq2::plotPCA(vst, ntop = 500, intgroup = "tissue"))
dev.off()

# png(paste0(outdir, "/QC/rld_dataset_pca.png"), width = 8, height=5, res = 300, units = "in")
# print(DESeq2::plotPCA(rld, ntop = 500, intgroup = "dataset"))
# dev.off()
png(paste0(outdir, "/QC/vst_dataset_pca.png"), width = 8, height=5, res = 300, units = "in")
print(DESeq2::plotPCA(vst, ntop = 500, intgroup = "dataset"))
dev.off()

# png(paste0(outdir, "/QC/rld_patient_pca.png"), width = 8, height=5, res = 300, units = "in")
# print(DESeq2::plotPCA(rld, ntop = 500, intgroup = "patient"))
# dev.off()
png(paste0(outdir, "/QC/vst_patient_pca.png"), width = 8, height=5, res = 300, units = "in")
print(DESeq2::plotPCA(vst, ntop = 500, intgroup = "subject"))
dev.off()

# png(paste0(outdir, "/QC/rld_batch_pca.png"), width = 8, height=5, res = 300, units = "in")
# print(DESeq2::plotPCA(rld, ntop = 500, intgroup = "submission"))
# dev.off()
png(paste0(outdir, "/QC/vst_cluster_pca.png"), width = 8, height=5, res = 300, units = "in")
print(DESeq2::plotPCA(vst, ntop = 500, intgroup = "cluster"))
dev.off()

# png(paste0(outdir, "/QC/rld_cellCount_pca.png"), width = 8, height=5, res = 300, units = "in")
# print(DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count"))
# dev.off()
# png(paste0(outdir, "/QC/vst_cellCount_pca.png"), width = 8, height=5, res = 300, units = "in")
# print(DESeq2::plotPCA(vst, ntop = 500, intgroup = "cell_count"))
# dev.off()

# Extract the rlog matrix from the object and compute pairwise correlation values
# rld_mat <- assay(rld)
# rld_cor <- cor(rld_mat)

vst_mat <- assay(vst)
vst_cor <- cor(vst_mat)

# Plot heatmap
pdf(paste0(outdir, "/QC/vst_corr_heatmap.pdf"), width = 12,height=15)
pheatmap(vst_cor, annotation = meta[, c("dataset", "dx", "cluster"), drop=F])
dev.off()
# pdf(paste0(outdir, "/QC/rld_corr_heatmap.pdf"), width = 8,height=7)
# pheatmap(rld_cor, annotation = meta[, c("patient"), drop=F])
# dev.off()

# Run DESeq2 differential expression analysis
dds <- DESeq(dds, fitType = "local")

# Plot dispersion estimates
pdf(paste0(outdir, "/QC/dispersionplot.pdf"), width = 8,height=5)
plotDispEsts(dds) #acceptable i guess? idk
dev.off()

saveRDS(dds, file = paste0(outdir,"/final_dds.Rds"))

res <- results(dds, contrast = c("cluster", "pathogenic_CD8T", "MAIT"))
View(as.data.frame(res))
write.table(as.data.frame(res), file = paste0(outdir, "/pathogenic_vs_MAIT_DESeq2_res_total.txt"), sep = "\t", row.names = T, col.names = T)


#############################################
#### Look at each cell type individually ####
#############################################

dds <- readRDS("/Volumes/Active/Isabel_Risch/MP021/data/de/DESeq2/pathogeniccluster_vs_MAIT_subject_tissue/final_dds.Rds")
rn <- resultsNames(dds)

# grab the 6 explicit cluster coefficients
cluster_coefs <- grep("^cluster", rn, value = TRUE)
cluster_coefs

cluster_labels <- sub("^cluster", "", cluster_coefs)     # strip prefix for readable names
all_clusters <- c(cluster_labels, "GZMK_effector_CD8T")  # add back the implicit reference

all_celltype_markers <- list()

for (ct in all_clusters) {
  outdir <- paste0("data/de/DESeq2/", ct, "_vs_others/")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  contrast_vec <- setNames(numeric(length(rn)), rn)   # subject/tissue terms stay 0
  
  if (ct == "GZMK_effector_CD8T") {
    # reference level: its own coefficient is implicitly 0
    # so this cluster's mean = 0 - mean(other 6 cluster coefficients)
    contrast_vec[cluster_coefs] <- -1/6
  } else {
    my_coef <- paste0("cluster", ct)
    other_coefs <- setdiff(cluster_coefs, my_coef)
    contrast_vec[my_coef]    <- 1
    contrast_vec[other_coefs] <- -1/6
    # GZMK's implicit 0 is the 6th "other" cluster — already correctly
    # folded into the denominator without needing an explicit term
  }
  
  res <- results(dds, contrast = contrast_vec)
  
  res2 <- as.data.frame(res)
  res2$gene <- rownames(res2)
  res2 <- res2[order(res2$log2FoldChange, decreasing = TRUE), ]
  res2$celltype <- ct
  
  write.table(res2, file = paste0(outdir, "/", ct, "_vs_others_DESeq2_res_total.txt"),
              sep = "\t", row.names = TRUE, col.names = TRUE)
  all_celltype_markers[[ct]] <- res2
}

per_cluster_markers <- character(0)
for (i in 1:length(all_celltype_markers)) {
  signif_markers <- all_celltype_markers[[i]][ all_celltype_markers[[i]]$padj < 0.01,]
  per_cluster_markers <- c(per_cluster_markers, signif_markers[1:5,"gene"])
}
per_cluster_markers <- unique(per_cluster_markers)
saveRDS(per_cluster_markers, "data/per_cluster_markers.Rds")
saveRDS(all_celltype_markers, "data/DESeq2_all_celltype_markers.Rds")

# export to an excel file for easy reporting
wb <- createWorkbook()
for (ct in names(all_celltype_markers)) {
  addWorksheet(wb, sheetName = ct)
  writeData(wb, sheet = ct, x = all_celltype_markers[[ct]], rowNames = TRUE)
}
saveWorkbook(wb, "data/de/DESeq2/all_celltype_markers.xlsx", overwrite = TRUE)

#############################################
#### MAIT vs. others plotting (volcanoplots) ####
#############################################
source("~/Documents/resources/scripts/universal_helpers/helpers.R")
source("~/Documents/resources/scripts/universal_helpers/netbid_volcano_edit.R")

outdir <- "data/de/DESeq2/MAIT_vs_others/"
res2 <- read.delim(paste0(outdir, "/MAIT_vs_others_DESeq2_res_total.txt"))
mait_marker_genes <- res2[(res2$log2FoldChange > 1) & (res2$padj < 0.01),"gene"]
gmtlist2file(gmtlist = list("mait_vs_others_4publicdataset" = mait_marker_genes), 
             filename = "data/4publicdataset_signatures.gmt")

draw.volcanoPlot.edit(dat = res2, label_col = "gene", logFC_col = "log2FoldChange", Pv_col = "padj", 
                      logFC_thre = 1, Pv_thre = 0.01, show_plot = TRUE, xlab = "log2 Fold Change", 
                      ylab = "P-value", show_label = T, label_cex = 0.5, legend_cex = 0.8, 
                      label_type = "distribute", main = "", pdf_file = paste0(outdir,"/MAIT_vs_others_labeled.pdf"), highlight_genes=NULL,
                      pipeline_script="/Users/irisch/Documents/software/NetBID-master/R/pipeline_functions.R") 

# res2 <- res2[order(res2$padj, decreasing = F),]
# highlight_genes <- c(head(res2[res2$log2FoldChange > 1,"gene"], 10),
#                      head(res2[res2$log2FoldChange < -1,"gene"], 10))

res2 <- res2[order(res2$log2FoldChange, decreasing = T),]
highlight_genes <- c(head(res2[res2$padj < 0.01,"gene"], 10),
                     tail(res2[res2$padj < 0.01,"gene"], 10))

draw.volcanoPlot.scaleFixed(dat = res2, label_col = "gene", logFC_col = "log2FoldChange", Pv_col = "padj", 
                            logFC_thre = 1, Pv_thre = 0.01, show_plot = TRUE, xlab = "log2 Fold Change", 
                            ylab = "P-value", show_label = T, label_cex = 0.5, legend_cex = 0.8, 
                            label_type = "distribute", main = "", 
                            pdf_file = paste0(outdir,"/MAIT_vs_others_labeledTop.pdf"), 
                            highlight_genes=highlight_genes,
                            pipeline_script="/Users/irisch/Documents/software/NetBID-master/R/pipeline_functions.R") 


