# .libPaths(c("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/R_libraries/4.0.3", 
#             "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
# 
# setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP003/")

setwd("/Volumes/Active/Isabel_Risch/MP003")

library(Seurat)
library(scRepertoire)
library(dplyr)

d <- readRDS("data/robjects/6_singletsNoBystanders_clustered.Rds")

#########################################
########## Read in TCR data #############
#########################################

S1 <- read.csv("/Volumes/Active/SeqData/MP156/Paley_vdj_multi_SR004645_10X/cellranger_v8.0.1/YeiH_Sort_1/outs/outs/per_sample_outs/YeiH_Sort_1/vdj_t/filtered_contig_annotations.csv")
S2 <- read.csv("/Volumes/Active/SeqData/MP156/Paley_vdj_multi_SR004645_10X/cellranger_v8.0.1/YeiH_Sort_2/outs/outs/per_sample_outs/YeiH_Sort_2/vdj_t/filtered_contig_annotations.csv")
S3 <- read.csv("/Volumes/Active/SeqData/MP156/Paley_vdj_multi_SR004645_10X/cellranger_v8.0.1/YeiH_Sort_3/outs/outs/per_sample_outs/YeiH_Sort_3/vdj_t/filtered_contig_annotations.csv")

#make barcodes unique
S1$barcode <- paste0("Sort1_", S1$barcode)
S2$barcode <- paste0("Sort2_", S2$barcode)
S3$barcode <- paste0("Sort3_", S3$barcode)

S1$contig_id <- paste0("Sort1_", S1$contig_id)
S2$contig_id <- paste0("Sort2_", S2$contig_id)
S3$contig_id <- paste0("Sort3_", S3$contig_id)

all <- rbind(S1, S2) %>% rbind(S3)

contig.list <- createHTOContigList(all,
                                   d,
                                   group.by = "HTO_dmm")

combined.TCR <- combineTCR(contig.list, 
                           samples = names(contig.list),
                           removeNA = FALSE, 
                           removeMulti = F, 
                           filterMulti = FALSE)

# for visualization, put the patients and the controls next to each other
combined.TCR <- combined.TCR[sort(names(combined.TCR))]

# clonalScatter(combined.TCR, 
#               cloneCall ="gene", 
#               x.axis = "P18B", 
#               y.axis = "P18L",
#               dot.size = "total",
#               graph = "proportion")

saveRDS(combined.TCR, file = "data/robjects/combined_tcrs.Rds")

#############################################
######## Add the TCR data to Seurat #########
#############################################

combine.tmp <- lapply(combined.TCR, function(x) {
  x[[1]] <- sapply(strsplit(x[[1]], "_Sort"), "[", 2)
  x[[1]] <- paste0("Sort", x[[1]])
  return(x)
  })
for (i in names(combine.tmp)) combine.tmp[[i]]$HTO_dmm <- i

sce <- combineExpression(combine.tmp, 
                         d, 
                         cloneCall="strict", 
                         group.by = "HTO_dmm",
                         cloneSize = c(Single = 1, Small = 5, 
                                       Medium = 50, Large = 150),
                         proportion = F)

DimPlot(sce, group.by = "cloneSize")

# too-strict pathogenic definition below, not what was used in manuscript
sce$pathogenic <- ifelse(grepl("^.{5,6}[L|T][Y|F]ST.{5}", sapply(strsplit(sce$CTaa, "_"), "[", 2)), TRUE, FALSE)
sce$pathogenic <- ifelse(grepl("TRAV21", sce$CTstrict), sce$pathogenic, FALSE)

DimPlot(sce, group.by = "pathogenic")

## Barplots to determine if all the pathogenic cells are from one or two patients
tmp <- table(d$pathogenic, d$HTO_dmm)
tmp <- reshape2::melt(tmp)

ggplot(tmp, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="stack", stat="identity")

ggplot(tmp, aes(fill=Var2, x=Var1, y=value)) + 
  geom_bar(position="fill", stat="identity")

# saveRDS(sce, file = "data/robjects/7_singletsNoBystanders_clustered_withTCR.Rds")
