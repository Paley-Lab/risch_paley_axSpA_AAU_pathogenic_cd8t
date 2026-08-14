# Plot the correlation of number of wells with beta sequence (paired) to number of reads (unpaired)
# to demonstrate that it can be used as a surrogate
library(ggplot2)

all_samples_beta <- read.csv(file = "data/inputs/all_samples_BETA_ONLY.csv")
paired <- read.csv(file = "data/inputs/all_samples_paired.csv")

betaOnlyCounts <- all_samples_beta[,c(1,2,12)]
colnames(betaOnlyCounts) <- c("sequence_BetaOnly", "counts_BetaOnly", "subject_BetaOnly")
pairedBetaCounts <- paired[,c(6,4,1,5)]
colnames(pairedBetaCounts) <- c("sequence_Paired", "counts_Paired", "subject_Paired", "alpha_Paired")

all <- merge(betaOnlyCounts, pairedBetaCounts, by.x = c("sequence_BetaOnly", "subject_BetaOnly"), by.y=  c("sequence_Paired", "subject_Paired"))
all <- unique(all)

write.csv(all, file = "data/outputs/correlate_beta_counts_paired_unpaired.csv", row.names = F)

p <- ggplot(all, aes(x = log10(counts_BetaOnly), y = log10(counts_Paired), color = subject_BetaOnly)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "Beta-Only Counts vs Paired Counts",
    x = "log10(counts_BetaOnly)",
    y = "log10(counts_Paired)",
    color = "Subject"
  ) +
  # scale_x_log10() + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("data/outputs/scatterplot_correlate_beta_counts_paired_unpaired.png", plot = p, width = 10, height = 7, dpi = 300)

# Are the CDR3b sequences that are being counted as 'pathogenic' different between the 
# patients and the healthy controls? (Make a SeqLogo?)


