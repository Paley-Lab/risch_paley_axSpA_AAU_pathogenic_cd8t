# run this on RIS with the following command:
# bsub -Is -M 128GB -R 'rusage[mem=128GB]' -G compute-gfwu -q oncology-interactive -a 'docker(igrisch/seurat:5.2.0_v2)' /bin/bash
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP021")
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

d <- readRDS("data/robjects/4datasets_integrated.Rds")

###########################################################################
########## Adapting Michael's code to check for "hidden" pathogenic TCRs ###########
###########################################################################

datasetLocations <- c(`deschler`="/storage1/fs1/paleym/Active/SeqData/public/SRP403458_deschler_JA_2022/cellranger_v8.0.1/",
                      `tang`="/storage1/fs1/paleym/Active/SeqData/public/GSE288581_tang/cellranger_v8.0.1/",
                      `polvoleri`="/storage1/fs1/paleym/Active/SeqData/public/GSE216914_polvoleri_2023/cellranger_v8.0.1/",
                      `paley_2` = "/storage1/fs1/paleym/Active/SeqData/Eye_Blood_Pairs/AB_cellranger_v8.0.1/cellranger_v8.0.1/")

#### Deschler dataset ####

des <- subset(d, dataset %in% "deschler")
# there's weird mapping between the SRA names and the sample names so we have to do funny stuff with the metadata
meta <- read.csv("data/deschler_metadata.csv")
seqs <- list()
for (i in unique(des$subject)) { # for each sample...
    sra <- meta[(meta$Sample.Name %in% i) & (meta$type =="tcr"),1] # pull the corresponding TCR SRA number for each sample name
    sra_gex <- meta[(meta$Sample.Name %in% i) & (meta$type =="gex"),1] # pull the corresponding GEX SRA number for each sample name
    seqs[[sra_gex]] <- read.csv(paste0("/storage1/fs1/paleym/Active/SeqData/public/SRP403458_deschler_JA_2022/cellranger_v8.0.1/",sra,"/outs/filtered_contig_annotations.csv")) # read in the TCR seq data
    seqs[[sra_gex]]$cell <- paste0(sra_gex, "_", seqs[[sra_gex]]$barcode)
    seqs[[sra_gex]]$cell <- stringr::str_replace_all(seqs[[sra_gex]]$cell, "-", "\\.")
}
seqs <- dplyr::bind_rows(seqs) # collapse everything into one dataframe for search purposes

seqs_1 <- seqs

#### Tang dataset ####

tmp <- subset(d, dataset %in% "tang")
seqs <- list()
files <- list.dirs(path="/storage1/fs1/paleym/Active/SeqData/public/GSE288581_tang/cellranger_v8.0.1/", recursive = F, full.names = F)
files <- files[files %in% unique(tmp$subject)] 
names(files) <- files
files <- as.list(files)
seqs <- lapply(files, function(x){read.csv(paste0("/storage1/fs1/paleym/Active/SeqData/public/GSE288581_tang/cellranger_v8.0.1/",x,"/outs/per_sample_outs/",x,"/vdj_t/filtered_contig_annotations.csv"))})

for (i in 1:length(seqs)) { # for each sample...
  seqs[[i]]$cell <- paste0(names(seqs)[i], "_", seqs[[i]]$barcode)
  seqs[[i]]$cell <- stringr::str_replace_all(seqs[[i]]$cell, "-", "\\.")
}
seqs <- dplyr::bind_rows(seqs) # collapse everything into one dataframe for search purposes
seqs_2 <- seqs
#### Polvoleri dataset ####

tmp <- subset(d, dataset %in% "polvoleri")
seqs <- list()
files <- list.dirs(path="/storage1/fs1/paleym/Active/SeqData/public/GSE216914_polvoleri_2023/cellranger_v8.0.1/", recursive = F, full.names = F)
files <- files[files %in% unique(tmp$subject)] 
names(files) <- files
files <- as.list(files)
seqs <- lapply(files, function(x){read.csv(paste0("/storage1/fs1/paleym/Active/SeqData/public/GSE216914_polvoleri_2023/cellranger_v8.0.1/",x,"/outs/per_sample_outs/",x,"/vdj_t/filtered_contig_annotations.csv"))})

for (i in 1:length(seqs)) { # for each sample...
  seqs[[i]]$cell <- paste0(tolower(names(seqs)[i]), "_", seqs[[i]]$barcode)
  seqs[[i]]$cell <- stringr::str_replace_all(seqs[[i]]$cell, "-", "\\.")
}
seqs <- dplyr::bind_rows(seqs) # collapse everything into one dataframe for search purposes
seqs_3 <- seqs

#### Paley uveitis dataset ####

tmp <- subset(d, dataset %in% "paley_2")
meta <- read.csv("data/paley_2_metadata.csv")
meta <- meta[grep(subject, meta$Subject),]

files <- paste0(meta$TCR_CRoutput, "/filtered_contig_annotations.csv") %>% as.list()
names(files) <- meta$Sample
seqs <- lapply(files, read.csv)

for (i in 1:length(seqs)) { # for each sample...
  seqs[[i]]$cell <- paste0(names(seqs)[i], "_", seqs[[i]]$barcode)
  seqs[[i]]$cell <- stringr::str_replace_all(seqs[[i]]$cell, "-", "\\.")
}
seqs <- dplyr::bind_rows(seqs) # collapse everything into one dataframe for search purposes
seqs_4 <- seqs


#### Universal across all datasets ####

seqs <- rbind(seqs_1, seqs_2) %>% rbind(seqs_3) %>% rbind(seqs_4)

# If the cell name starts with a number, add an "X" to the beginning to make it match the Seurat obj
addX <- grep("^[0-9]", seqs$cell)
seqs[addX,"cell"] <- paste0("X",seqs[addX,"cell"])

# Michael's code. Within each cell, looks for a pathogenic TCR and will return TRUE if any pathogenic TCRs are found.
airr <- seqs %>% 
  group_by(cell) %>%  
  mutate(   
    pathogenic.2 = any( str_detect( v_gene , "TRAV21" ) ) & any( str_detect( cdr3[chain == "TRB"] , "^.{6,7}[Y|F]ST.{5}" ) )  
  ) %>%  
  ungroup()

airr$clonotype_id <- paste0(substr(airr$cell,1,nchar(airr$cell)-19), "_", airr$raw_clonotype_id)

# Find clonotypes with both TRUE and FALSE entries for the pathogenic motif
clonotype_mixed <- airr %>%
  group_by(clonotype_id) %>%
  summarize(
    has_true = any(pathogenic.2 == TRUE),
    has_false = any(pathogenic.2 == FALSE)
  ) %>%
  filter(has_true & has_false)
clonotype_mixed


add <- airr[,c("cell", "pathogenic.2")] %>% unique()
add <- as.data.frame(add)
rownames(add) <- add$cell

d <- AddMetaData(d, metadata = add)



###########################################################################
########## Getting a dataframe of TCR sequences to check ###########
###########################################################################

# pull the metadata
meta <- d@meta.data

# # get a unified disease column
# meta$disease_merged <- ifelse(!is.na(meta$disease_merged), meta$disease_merged, 
#                               ifelse(!is.na(meta$Anatomy), meta$Anatomy, 
#                                      "PsA"))
write.table(meta, file = "data/4datasets_integrated_withSig_meta.tsv", sep="\t", row.names = T, col.names = T)

# get the mean pathogenic score for each clone
means <- meta %>%
  group_by(CTstrict, subject) %>%
  dplyr::summarize(Mean_pathogenic_score = mean(Core_JCI_MP156_c4_vs_c3_noTCR_20g, na.rm=TRUE))
# add the mean pathogenic signature score data to the metadata
meta <- merge(meta, means, by = c("CTstrict", "subject"))
#subset to the pathogenic cluster
meta_sub <- meta[meta$celltype1=="pathogenic_CD8T",]
#subset to only those cells with clonality > 1
meta_sub <- meta_sub[meta_sub$clonalFrequency > 1, ]
#subset to only those cells that use TRAV21
meta_sub <- meta_sub[grep("TRAV21", meta_sub$CTgene),]
#remove cells with no beta chain
meta_sub <- meta_sub[!grepl("_NA", meta_sub$CTaa),]
#remove cells with multiple chains
meta_sub <- meta_sub[!grepl("\\;", meta_sub$CTgene),]
# subset to only the columns we're interested in and unique
meta_sub_unique <- meta_sub[,c("tissue_consensus", "subject", "pathogenic", "CTgene", 
                               "CTnt", "CTaa", "CTstrict", "clonalFrequency", "cloneSize", 
                               "disease_merged", "Mean_pathogenic_score", "dataset")] %>% unique()

# get the total number of cells in each sample-by-tissue
tot_cells <- meta[,c("subject", "tissue_consensus", "CTstrict")]
tot_cells <- tot_cells[complete.cases(tot_cells),] %>% group_by(subject, tissue_consensus) %>%
  dplyr::summarize(totalCells=n())

# add a column to meta_sub_unique that says what proportion of the cells with clonal information this constitutes
meta_sub_unique <- merge(meta_sub_unique, tot_cells, by=c("subject", "tissue_consensus"))
meta_sub_unique$cloneProportion <- meta_sub_unique$clonalFrequency/meta_sub_unique$totalCells

# identify the clones that are found across multiple tissues
collapse <- meta_sub_unique[duplicated(meta_sub_unique$CTstrict),"CTstrict"]
meta_sub_unique$multitissue <- meta_sub_unique$CTstrict %in% collapse

# prioritize clones that are more expanded in the tissues than in the PBMCs, where possible
meta_sub_unique$keep <- TRUE
for (i in 1:nrow(meta_sub_unique)) {
  if(!meta_sub_unique[i,"multitissue"]) next
  if(!meta_sub_unique[i,"keep"]) next
  tcr <- meta_sub_unique[i,"CTstrict"]
  pbmc <- meta_sub_unique[meta_sub_unique$CTstrict == tcr & meta_sub_unique$tissue_consensus=="PBMC","cloneProportion"]
  other <- meta_sub_unique[meta_sub_unique$CTstrict == tcr & meta_sub_unique$tissue_consensus!="PBMC","cloneProportion"]
  if (pbmc > other) {
    meta_sub_unique[meta_sub_unique$CTstrict == tcr,"keep"] <- FALSE
  }
}

meta_sub_unique <- meta_sub_unique[meta_sub_unique$keep,]
meta_sub_unique$keep <- NULL

# # OK, looks like all of these are from just one dataset soooo let's do this the ugly way
meta_sub_unique$tissue_consensus <- as.character(meta_sub_unique$tissue_consensus)
meta_sub_unique[meta_sub_unique$multitissue,"tissue_consensus"] <- "PBMC,Synovial_Fluid"
meta_sub_unique <- unique(meta_sub_unique)


meta_sub_unique$TCRa_AA <- sapply(strsplit(meta_sub_unique$CTaa, "_"), "[", 1)
meta_sub_unique$TCRb_AA <- sapply(strsplit(meta_sub_unique$CTaa, "_"), "[", 2)
meta_sub_unique <- tidyr::separate(meta_sub_unique, col="CTgene", sep="\\.|_", into=c("TRAV", "TRAJ", NA, "TRBV", NA, "TRBJ", "TRBC"), remove = F)

# # check if multi-tissue clones are more expanded in the inflamed tissue than the blood
# multitissue <- meta_sub_unique[meta_sub_unique$multitissue,]

#remove that tissue proportion  column to prevent duplicates
meta_sub_unique$cloneProportion <- NULL
meta_sub_unique$totalCells <- NULL
meta_sub_unique <- unique(meta_sub_unique)

# add quartile info
summary_info <- summary(meta_sub_unique$Mean_pathogenic_score)
meta_sub_unique$quartile <- ifelse(meta_sub_unique$Mean_pathogenic_score >= summary_info[5], "Q4",
                                   ifelse(meta_sub_unique$Mean_pathogenic_score >= summary_info[3], "Q3",
                                          ifelse(meta_sub_unique$Mean_pathogenic_score >= summary_info[2], "Q2",
                                          "Q1")))

#sort by descending path. score
meta_sub_unique <- meta_sub_unique[order(meta_sub_unique$Mean_pathogenic_score, decreasing = T),]

# give each clone a number
for (i in unique(meta_sub_unique$subject)) {
  meta_sub_unique[meta_sub_unique$subject %in% i,"subject"] <- paste0(meta_sub_unique[meta_sub_unique$subject %in% i,"subject"], ".", 1:length(meta_sub_unique[meta_sub_unique$subject %in% i,"subject"]))
}

# print out to csv
write.table(meta_sub_unique, file = "data/4datasets_integrated_TCRs_of_interest.csv", row.names = F, sep = ",")


####### Run barplots locally for RColorBrewer #####

meta_sub_unique <- openxlsx::read.xlsx("data/4datasets_integrated_TCRs_of_interest_handAnnotated.xlsx")

library(RColorBrewer)
f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
(cols <- f("Paired")) %>% dput()

## Barplots of TCR usage by quartile
df <- as.data.frame(table(meta_sub_unique$TRBV, meta_sub_unique$quartile))
colnames(df) <- c("TRBV", "Quartile", "Frequency")
p<-ggplot(df, aes(x=Quartile, y=Frequency, fill=TRBV)) +
  geom_bar(stat="identity", position = "fill")+theme_classic()+
  # labs(title="Title") +
  theme(plot.title = element_text(hjust = 0.5, size=16)) + scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99","#B15928", 
                                                                                        "#FDBF6F", "#FF7F00",  "#6A3D9A", "#FFFF99",  "#E31A1C",
                                                                                        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
                                                                                        "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#FFED6F"))
# geom_text(aes(label=round(`-Log10(p-value)`, digits=1)), vjust=1.6, color="white", size=6)
p
ggsave(plot=p, filename = "plot/barplot_trbv_usage_by_quartile.pdf", device="pdf", width=4, height=2)

df <- as.data.frame(table(meta_sub_unique$pathogenic, meta_sub_unique$quartile))
colnames(df) <- c("Pathogenic", "Quartile", "Frequency")
df$pathogenic <- factor(df$Pathogenic, levels=c("FALSE", "maybe", "TRUE"))
p<-ggplot(df, aes(x=Quartile, y=Frequency, fill=Pathogenic)) +
  geom_bar(stat="identity", position = "fill")+theme_classic()+
  # labs(title="Title") +
  theme(plot.title = element_text(hjust = 0.5, size=16))  +
  scale_fill_manual(values=c("lightgrey","hotpink", "red"))
# geom_text(aes(label=round(`-Log10(p-value)`, digits=1)), vjust=1.6, color="white", size=6)
p
ggsave(plot=p, filename = "plot/barplot_pathogenicity_by_quartile.pdf", device="pdf", width=4, height=2)


####### automatically pull the full TCR sequences that we want to validate #####
rm(list=ls())
tcrs <- openxlsx::read.xlsx("data/4datasets_integrated_TCRs_of_interest_handAnnotated.xlsx")
tcrs <- tcrs[grep("yes", tcrs$Validate_new, ignore.case = T),]

# datasetLocations <- c(`deschler`="/Volumes/Active/SeqData/public/SRP403458_deschler_JA_2022/cellranger_v8.0.1/",
#                       `tang`="/Volumes/Active/SeqData/public/GSE288581_tang/cellranger_v8.0.1/",
#                       `polvoleri`="/Volumes/Active/SeqData/public/GSE216914_polvoleri_2023/cellranger_v8.0.1/",
#                       `paley_2` = "/Volumes/Active/SeqData/Eye_Blood_Pairs/AB_cellranger_v8.0.1/cellranger_v8.0.1/")

for (i in 1:nrow(tcrs)) {
  
  subject <- sapply(strsplit(tcrs[i,"subject"], "\\."), "[", 1)
  
  if (tcrs[i,"dataset"] == "deschler") {
    meta <- read.csv("data/deschler_metadata.csv")
    meta <- meta[(meta$Sample.Name %in% subject) & (meta$type =="tcr"),]
    seqs <- read.delim(paste0("/Volumes/Active/SeqData/public/SRP403458_deschler_JA_2022/cellranger_v8.0.1/",meta[1,1],"/outs/airr_rearrangement.tsv"))
    
  } else if (tcrs[i,"dataset"] == "tang") {
    files <- list.dirs(path="/Volumes/Active/SeqData/public/GSE288581_tang/cellranger_v8.0.1/", recursive = F, full.names = F)
    files <- files[grep(subject,files)] %>% as.list()
    seqs <- lapply(files, function(x){read.delim(paste0("/Volumes/Active/SeqData/public/GSE288581_tang/cellranger_v8.0.1/",x,"/outs/per_sample_outs/",x,"/vdj_t/airr_rearrangement.tsv"))})
    seqs <- dplyr::bind_rows(seqs)
    
  } else if (tcrs[i,"dataset"] == "polvoleri") {
    seqs <- read.delim(paste0("/Volumes/Active/SeqData/public/GSE216914_polvoleri_2023/cellranger_v8.0.1/",subject,"/outs/per_sample_outs/",subject,"/vdj_t/airr_rearrangement.tsv"))
    
  } else if (tcrs[i,"dataset"] == "paley_2") {
    meta <- read.csv("data/paley_2_metadata.csv")
    meta <- meta[grep(subject, meta$Subject),]
    files <- paste0(stringr::str_replace_all(meta$TCR_CRoutput, "storage1/fs1/paleym", "Volumes"), "/airr_rearrangement.tsv") %>% as.list()
    seqs <- lapply(files, read.delim)
    seqs <- dplyr::bind_rows(seqs)
    
  }
  
  # filter down to those cells with one and ONLY ONE alpha chain, and with one and ONLY ONE beta chain
  paired <- seqs[duplicated(seqs$cell_id), "cell_id"]
  
  seqs <- seqs[seqs$cell_id %in% paired,]
  
  seqs$num_unique_seqs <- NA
  
  for (j in unique(seqs$cell_id)) {
    seqs[seqs$cell_id %in% j, "num_unique_seqs"] <- length(unique(seqs[seqs$cell_id %in% j, "sequence_aa"]))
  }
  
  # seqs2 <- seqs[seqs$num_unique_seqs == 2, ]
  
  tcra <- sapply(strsplit(tcrs[i, "CTnt"], "_"),"[", 1)
  tcrb <- sapply(strsplit(tcrs[i, "CTnt"], "_"),"[", 2)
  
  full_seq_alpha <- seqs[grepl(tcra, seqs$junction) & grepl("TRA", seqs$v_call),]
  full_seq_beta <- seqs[grepl(tcrb, seqs$junction) & grepl("TRB", seqs$v_call),]
  
  full_seq_alpha <- full_seq_alpha[grep("TRA", full_seq_alpha$c_call),]
  full_seq_beta <- full_seq_beta[grep("TRB", full_seq_beta$c_call),]
  
  full_seq_alpha <- full_seq_alpha[full_seq_alpha$cell_id %in% full_seq_beta$cell_id,]
  full_seq_beta <- full_seq_beta[full_seq_beta$cell_id %in% full_seq_alpha$cell_id,]
  
  if ((length(unique(full_seq_alpha$sequence_aa)) != 1) | (length(unique(full_seq_beta$sequence_aa)) != 1)) next
  
  tcrs[i,"full_tcra"] <- unique(full_seq_alpha$sequence_aa)
  tcrs[i,"full_tcrb"] <- unique(full_seq_beta$sequence_aa)
  
}

openxlsx::write.xlsx(tcrs, "data/4datasets_integrated_TCRs_of_interest_handAnnotated_withFullTCR.xlsx")












