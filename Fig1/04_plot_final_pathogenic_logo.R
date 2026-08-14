setwd("/Volumes/Active/Isabel_Risch/MP033")

library(ggseqlogo)
library(ggplot2)

# ---- build a custom position-frequency matrix for positions 6-10 ----------
AA_ORDER <- c('A','R','N','D','C','Q','E','G','H','I',
              'L','K','M','F','P','S','T','W','Y','V')

positions <- as.character(1:15)

pfm <- matrix(0, nrow = length(AA_ORDER), ncol = length(positions),
              dimnames = list(AA_ORDER, positions))

# position 1: C only
pfm["C", "1"] <- 1

# position 2: A only
pfm["A", "2"] <- 1

# position 3: S only
pfm["S", "3"] <- 1

# position 4: S only
pfm["S", "4"] <- 1

# position 5: any amino acid
pfm[c('A','R','N','D','C','Q','E','G','H','I',
      'L','K','M','F','P','S','T','W','Y','V'), "5"] <- 20

# position 6: A, G, D, S, E, Q, R in equal proportions
pfm[c("A","G","D","S","E","Q","R"), "6"] <- 1

# position 7: F, R, L, Y, H, V, I, T, S, M in equal proportions
pfm[c("F","R","L","Y","H","V","I","T","S","M"), "7"] <- 1

# position 8: Y, F in equal proportions
pfm[c("Y","F"), "8"] <- 1

# position 9: S only
pfm["S", "9"] <- 1

# position 10: T only
pfm["T", "10"] <- 1

# position 11: D only
pfm["D", "11"] <- 1

# position 12: T only
pfm["T", "12"] <- 1

# position 13: Q only
pfm["Q", "13"] <- 1

# position 14: Y only
pfm["Y", "14"] <- 1

# position 15: F only
pfm["F", "15"] <- 1

# ---- plot as proportions (not bits) ----------------------------------------
p <- ggseqlogo(pfm, method = "probability") +
  labs(x = "Position", y = "Proportion") +
  ggtitle("Custom motif: positions 6-10")

ggsave(plot=p, filename = "outputs/final_pathogenic_search_logo.pdf", width = 6.5, height = 2)
