# run this on RIS with the following command:
# bsub -Is -M 128GB -R 'rusage[mem=128GB]' -G compute-gfwu -q general-interactive -a 'docker(igrisch/seurat:5.2.0_v2)' /bin/bash
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP021")
library(Seurat)
library(ggplot2)
library(dplyr)

#################################################
#### Remove the YeiH sort data and recluster ####
#################################################

mat <- read.delim("data/scanpy/raw_counts_all.tsv")
meta <- read.delim("data/scanpy/meta.tsv")
d <- CreateSeuratObject(counts=mat, meta.data = meta)
d <- subset(d, dataset != "paley_1")


## Add annotation based on disease status
d$disease_merged <- ifelse(!is.na(d$disease_merged), d$disease_merged, 
                           ifelse(!is.na(d$Anatomy), d$Anatomy, "PsA"))
d <- subset(d, disease_merged != "PanU")

## Add annotation based on B27 status
d$B27 <- ifelse(grepl("ositive",d$HLA_B27), "B27+",
                ifelse(grepl("ositive",d$`HLA.B27`), "B27+",
                       ifelse(grepl("egative",d$`HLA.B27`), "B27-",
                           ifelse(grepl("PSA3|PSA4", d$subject), "B27-", # in the Polvoleri (PsA) data, PSA3 and PSA4 were B27-
                                  ifelse(d$subject %in% "76", "B27-","B27+") # in the Tang data, only patient 76 was B27-
                                  )
                       )
                       ))
# remove redundant/confusing columns 
d$HLA_B27 <- NULL
d$`HLA.B27` <- NULL
d$Anatomy <- NULL
d$Etiology <- NULL
d$disease <- NULL
d$Tissue_1 <- NULL
d$Tissue_2 <- NULL
d$tissue <- NULL
d$`Sample.Name` <- NULL
d$mait <- NULL
d$mait2 <- NULL
d$rough_celltype <- NULL
d$celltype2 <- NULL
d$Sample <- NULL
d$Subject <- NULL

# split by your variable of interest
d[["RNA"]] <- split(d[["RNA"]], f = d$subject)
# run standard workflow
d <- NormalizeData(d)
d <- FindVariableFeatures(d)
d <- ScaleData(d)
d <- RunPCA(d)

# After preprocessing, we integrate layers with added parameters specific to Harmony:
d <- IntegrateLayers(object = d, method = HarmonyIntegration, orig.reduction = "pca",
                     new.reduction = 'harmony', verbose = FALSE)

d <- FindNeighbors(d, reduction = "harmony", dims = 1:15)
d <- RunUMAP(d, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
d <- FindClusters(d, resolution = 0.3, cluster.name = "harmony_clusters")

p1 <- DimPlot(d,reduction = "umap.harmony",group.by = "dataset", combine = FALSE, label.size = 2, shuffle = T)
ggsave(plot=p1[[1]], filename="plot/umap_dataset.pdf", width=6,height=4,units="in")

p1 <- DimPlot(d,reduction = "umap.harmony",group.by = "harmony_clusters", combine = FALSE, label.size = 2, shuffle = T)
ggsave(plot=p1[[1]], filename="plot/umap_clusters.pdf", width=6,height=4,units="in")

p1 <- FeaturePlot(d, features = c("PDCD1", "HAVCR2", "KLRB1", "KIR3DL1", "KLRC3", "FCER1G", "KIR2DL4", "KLRC2"), ncol=4)
ggsave(plot=p1, filename="tmp2.pdf", width=18,height=8,units="in")

new.cluster.ids <- c(`0` = "GZMK_effector_CD8T",
                     `1` = "GZMK_effector_CD8T",
                     `2` = "naive_CD8T",
                     `3` = "MAIT",
                     `4` = "NKgene_CD8T",
                     `5` = "GZMK_effector_CD8T",
                     `6` = "GNLY_effector_CD8T",
                     `7` = "GZMK_effector_CD8T",
                     `8` = "pathogenic_CD8T",
                     `9` = "cycling")
names(new.cluster.ids) <- levels(d)
d <- RenameIdents(d, new.cluster.ids)

saveRDS(d, file = "data/robjects/4datasets_integrated.Rds")





