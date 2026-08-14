## =======================================
## Log-ratio of pre-screen replicates 
## =======================================
library(dplyr)
library(tidyverse)
setwd("/Volumes/Active/Isabel_Risch/MP033/")

MIN_COUNT <- 1

df <- read.csv("../MP027/outputs/count_matrix_translated.csv", row.names = 1)

#### CPM normalize 
cpm_norm <- function(df) {
  # Each sample (column) divided by its own total counts, scaled to per-million
  lib_sizes <- colSums(df)
  cpm <- sweep(df, 2, lib_sizes, FUN = "/") * 1e6
  return(cpm)
}

cpm_df <- cpm_norm(df)

#### Add a pseudocount to avoid undefined values
enrich_df <- cpm_df + 1

log_ratio_df <- enrich_df %>%
  filter(Pre.YST_Pre.screen_Rep1 >= MIN_COUNT, Pre.YST_Pre.screen_Rep2 >= MIN_COUNT) %>%
  transmute(log_ratio = log2(Pre.YST_Pre.screen_Rep1 / Pre.YST_Pre.screen_Rep2))

mu    <- 0 #mean(log_ratio_df$log_ratio)
sigma <- sd(log_ratio_df$log_ratio)

ggplot(log_ratio_df, aes(log_ratio)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80,
                 fill = "#F8D686", color = "white", alpha = 0.9) +
  stat_function(fun = dnorm, args = list(mean = mu, sd = sigma),
                color = "firebrick", linewidth = 1) +
  coord_cartesian(xlim = c(-2, 2)) +
  labs(title = paste0("log2(startingLibraryPilot1/startingLibraryPilot2) (n=", nrow(log_ratio_df), ")"),
       x = "log(Rep1/Rep2)", y = "density") +
  theme_bw(base_size = 12)
ggsave("outputs/logNormal_plot.pdf", width=4,height=3,units="in")
