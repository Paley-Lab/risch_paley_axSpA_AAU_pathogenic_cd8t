library(ggplot2)
library(tidyr)
library(dplyr)
library(Biostrings)

# ── 0. Translate count matrix and collapse redundant sequences ──────────────────────────────────────────────────
count_matrix <- read.csv("outputs/count_matrix.csv", row.names = 1)

## collapse by sample
count_matrix <- t(count_matrix) %>% as.data.frame()
count_matrix$sample <- sapply(strsplit(rownames(count_matrix), "_lane"), "[", 1)
count_matrix_sampleCollapsed <- aggregate(
  . ~ count_matrix[["sample"]],
  data = count_matrix[, !colnames(count_matrix) %in% "sample", drop = FALSE],
  FUN = sum
)

# Set as row index and drop the column
rownames(count_matrix_sampleCollapsed) <- count_matrix_sampleCollapsed[[1]]
count_matrix_sampleCollapsed[[1]] <- NULL
count_matrix_sampleCollapsed <- t(count_matrix_sampleCollapsed)
write.csv(count_matrix_sampleCollapsed, file = "outputs/count_matrix_lanesCollapsed.csv", row.names = T, quote = F)


count_matrix <- count_matrix_sampleCollapsed %>% as.data.frame()

# -------- Translate the sequences in the matrix to amino acids --------- #
count_matrix$aa_seq <- Biostrings::translate(DNAStringSet(rownames(count_matrix)), no.init.codon = TRUE) %>%
  as.character() %>% paste0("_")
count_matrix <- count_matrix[order(count_matrix$aa_seq),]
count_matrix$aa_seq <- make.unique(count_matrix$aa_seq, sep = "")
count_matrix[grep("_$", count_matrix$aa_seq), "aa_seq"] <- paste0(count_matrix[grep("_$", count_matrix$aa_seq), "aa_seq"], "0")
rownames(count_matrix) <- count_matrix$aa_seq

# subset count matrix to only the columns with actual count_matrix in them
count_matrix <- count_matrix[,1:(ncol(count_matrix)-1)]

write.csv(count_matrix, file = "outputs/count_matrix_translated.csv", row.names = T, quote = F)


rm(list = ls())

df <- read.csv("outputs/count_matrix_translated.csv", row.names = 1, check.names = FALSE)
df$aa <- sapply(strsplit(rownames(df), "_"), "[", 1)

## collapse by amino acid sequence
collapsed <- aggregate(
  . ~ df[["aa"]],
  data = df[, !colnames(df) %in% "aa", drop = FALSE],
  FUN = sum
)

# Set as row index and drop the column
rownames(collapsed) <- collapsed[[1]]
collapsed[[1]] <- NULL

# save
write.csv(collapsed, "outputs/collapsed_count_matrix.csv", row.names = TRUE, quote=F)

