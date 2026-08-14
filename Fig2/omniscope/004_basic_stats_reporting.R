.libPaths(c("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/R_libraries/4.4.0", 
            "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP002/")

# library(Seurat)
library(openxlsx)

# ----- Read in the beta-only and paired data ------ # 
beta <- readRDS("data/beta_only_chains.Rds")
alpha <- readRDS("data/alpha_only_chains.Rds")
paired_reshaped <- readRDS("data/paired_object.Rds")


# ----- Get average number of unique sequences per sample (TCRb, TCRa, paired) ------ # 

# On average, CAIRR-seq returned [X] unique TCRb sequences 
mean_unique_beta <- unique(beta[,c("sample", "contig_id")]) %>% 
  group_by(sample) %>% 
  summarise(numberUniqueBetas = n())

# and [X] unique TCRa sequences per sample, 
mean_unique_alpha <- unique(alpha[,c("sample", "contig_id")]) %>% 
  group_by(sample) %>% 
  summarise(numberUniqueAlphas = n())

# as well as [X] unique paired TCRs per sample
mean_unique_paired <- unique(paired_reshaped[,c("sample", "contig_id_alpha", "contig_id_beta")]) %>% 
  group_by(sample) %>% 
  summarise(numberUniquePaired = n())

# Merge all tables and save
total_unique_sequences_table <- merge(mean_unique_alpha, mean_unique_beta) %>% merge(mean_unique_paired)
write.csv(total_unique_sequences_table, file = "data/CAIRRseq_unique_sequence_counts_per_sample.csv", row.names = F)



