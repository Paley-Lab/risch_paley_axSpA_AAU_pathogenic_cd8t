
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP028")
all_samples_beta <- list()
all_samples_alpha <- list()
input_dir <- "data/inputs/SRR011138/beta_chain_only_data/"
sample_list <- list.files(input_dir)
sample_list <- sapply(strsplit(sample_list, "_"), "[", 1) %>% unique()
for(i in sample_list) {
  all_samples_beta[[i]] <- read.delim(gzfile(paste0(input_dir, i, "_pseudobulk_TRB.tsv.gz")))
  all_samples_beta[[i]]$sample <- i
  all_samples_alpha[[i]] <- read.delim(gzfile(paste0("data/inputs/SRR011138/alpha_chain_only_data/", i, "_pseudobulk_TRA.tsv.gz")))
  all_samples_alpha[[i]]$sample <- i
}
input_dir <- "data/inputs/SRR009265/samples/"
sample_list <- list.files(input_dir)
for(i in sample_list) {
  all_samples_beta[[paste0(i, "_2")]] <- read.delim(gzfile(paste0(input_dir, i,"/", i, "_pseudobulk_TRB.tsv.gz")))
  all_samples_beta[[paste0(i, "_2")]]$sample <- paste0(i, "_2")
  all_samples_alpha[[paste0(i, "_2")]] <- read.delim(gzfile(paste0(input_dir, i,"/", i, "_pseudobulk_TRA.tsv.gz")))
  all_samples_alpha[[paste0(i, "_2")]]$sample <- paste0(i, "_2")
}

## collapse everything into one dataframe
all_samples_beta <- dplyr::bind_rows(all_samples_beta)
all_samples_alpha <- dplyr::bind_rows(all_samples_alpha)

write.csv(all_samples_alpha, "data/outputs/all_alpha_data.csv")
# write.csv(all_samples_beta, "data/outputs/all_beta_data.csv") #redundant with the csv saved in "02_beta_only_analysis", figure out deduplication


# ----- Get average number of unique sequences per sample (TCRb, TCRa, paired) ------ # 

all_samples_alpha <- read.csv("data/inputs/all_alpha_data.csv")
all_samples_beta <- read.csv("data/inputs/all_samples_BETA_ONLY.csv")
all_samples_paired <- read.csv("data/inputs/all_samples_paired.csv")

# On average, CAIRR-seq returned [X] unique TCRb sequences 
mean_unique_beta <- unique(all_samples_beta[,c("sample", "targetSequences", "v", "j", "aaSeqCDR3")]) %>% 
  group_by(sample) %>% 
  summarise(numberUniqueBetas = n())

# and [X] unique TCRa sequences per sample, 
mean_unique_alpha <- unique(all_samples_alpha[,c("sample", "targetSequences", "v", "j", "aaSeqCDR3")]) %>% 
  group_by(sample) %>% 
  summarise(numberUniqueAlphas = n())

# as well as [X] unique paired TCRs per sample
mean_unique_paired <- unique(all_samples_paired[,c("sample", "alpha_nuc", "beta_nuc", "cdr3a", "va", "da", "ja", "ca", "cdr3b", "vb", "db", "jb", "cb")]) %>% 
  group_by(sample) %>% 
  summarise(numberUniquePaired = n())

# Merge all tables and save
total_unique_sequences_table <- merge(mean_unique_alpha, mean_unique_beta) %>% merge(mean_unique_paired)
write.csv(total_unique_sequences_table, file = "data/TIRTLseq_unique_sequence_counts_per_sample.csv", row.names = F)
