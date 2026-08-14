library(dplyr)
setwd("/Volumes/Active/Isabel_Risch/MP033/")

rep1 <- read.csv("/Volumes/Active/Isabel_Risch/MP035/outputs/count_matrix_translated.csv", row.names = 1)
rep2 <- read.csv("/Volumes/Active/Isabel_Risch/MP031/outputs/count_matrix_translated.csv", row.names = 1)
colnames(rep2) <- stringr::str_replace_all(colnames(rep2), "rep1", "rep2")

# ------- What % of the library is degenerate? -------- #
allCodons <- rownames(rep1)
allCodons2<- sapply(strsplit(allCodons, "_"), "[", 1)
# allCodons2 <- allCodons2[-grep("\\*", allCodons2)]
allCodons2 <- table(allCodons2) %>% as.data.frame()
table(allCodons2$Freq!=1)
# 7063/(7063+2197)
# 0.762743

# ------- In each replicate: -------- #
#### CPM normalize 
cpm_norm <- function(df) {
  # Each sample (column) divided by its own total counts, scaled to per-million
  lib_sizes <- colSums(df)
  cpm <- sweep(df, 2, lib_sizes, FUN = "/") * 1e6
  return(cpm)
}

cpm_rep1 <- cpm_norm(rep1)
cpm_rep2 <- cpm_norm(rep2)

#### Add a pseudocount to avoid undefined values
enrich_rep1 <- cpm_rep1 + 1
enrich_rep2 <- cpm_rep2 + 1

#### Now take the ratios of selected to unselected
enrich_rep1 <- enrich_rep1[,c("Pre_YST_rep3_YeiH", "Pre_YST_rep3_GPER1", 
                              "Pre_YST_rep3_PRPF3", "Pre_YST_rep3_PSG5",  "Pre_YST_rep3_RNASEH2b", "Pre_YST_rep3_no_peptide"
                              )] / enrich_rep1[,c("Pre_YST_rep3_startingLibrary")]
enrich_rep2 <- enrich_rep2[,c("Pre_YST_library_rep2_yeih","Pre_YST_library_rep2_GPER1", "Pre_YST_library_rep2_PRPF3", 
                           "Pre_YST_library_rep2_PSG5", "Pre_YST_library_rep2_RNASEH2b","Pre_YST_library_rep2_no_peptide"
                           )] / enrich_rep2[,c("Pre_YST_library_rep2_preSort")]


#### Null hypothesis: the ratio is 1

##### Inspect the distributions. Enrich uses a 2-sided Poisson exact test, 
##### which means when we see our density distribution of the pre-enrich samples, it should be a Poisson
##### distribution centered around 1

# ── Reshape to long format ─────────────────────────────────────────────────
long_counts <- cbind(enrich_rep1,enrich_rep2) |>
  pivot_longer(everything(), names_to = "sample", values_to = "count") |>
  mutate(
    group = case_when(
      grepl("rep3", sample)  ~ "rep3",
      grepl( "rep2", sample)  ~ "Rep2",
      TRUE ~ "Unknown"
    )
  )

# ── Plot density curves ──────────────────────────────
ggplot(long_counts, aes(x = log2(count), colour = sample, fill = sample)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    name   = expression(log[2](enrichmentRatio))
  ) +
  scale_y_continuous(name = "Density") +
  labs(
    title   = "Enrichment Ratio Distribution",
    colour  = "Sample",
    fill    = "Sample"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )


#### Okay, so we clearly see when we LOG normalize, there's a deselected and a selected population, and based
#### on the distribution we see with the pre-selection counts, our null hypothesis is that no selection occurred 
#### either way and the log-transformed values are just a nice normal distribution.
#### We can do a T-test on this if we log transform the counts

allEnrichScores <- cbind(enrich_rep1, enrich_rep2)

logEnrichScores <- log2(allEnrichScores)

### For each set of degenerate codons in each condition, do a T test to see if the values differ significantly from 1 
logEnrichScores$aa <- sapply(strsplit(rownames(logEnrichScores), "_"), "[", 1)

allSamplesExceptNoPeptide_plots <- list()
allSamplesExceptNoPeptide_res <- lapply(1:6, function(x) data.frame(`triplet` = NA, `p.val` = NA, avg_log2FC = NA, `n_obs` = NA))
names(allSamplesExceptNoPeptide_res) <- c("YeiH", "RNASEH2B", "GPER1", "PRPF3", "PSG5", "allAntigens")
for (triplet in unique(logEnrichScores$aa)) {
  tmp <- logEnrichScores[logEnrichScores$aa %in% triplet, !grepl("no_peptide", colnames(logEnrichScores))]
  melted <- suppressMessages(reshape2::melt(tmp))
  melted$group <- ifelse(grepl("yeih", melted$variable, ignore.case = T), "YeiH", 
                         ifelse(grepl("rnaseh2b", melted$variable, ignore.case = T), "RNASEH2B", 
                                ifelse(grepl("gper1", melted$variable, ignore.case = T), "GPER1", 
                                       ifelse(grepl("PRPF3", melted$variable, ignore.case = T), "PRPF3", 
                                              ifelse(grepl("psg5", melted$variable, ignore.case = T), "PSG5", 
                                              "other")))))
  melted$replicate <- ifelse(grepl("rep3", melted$variable, ignore.case = T), "Rep3", "Rep2")
  
  # ----- T test to determine if values are different from 0 ------- #
  test_allgroups <- t.test(melted$value, mu = 0)
  avg_log2FC_allgroups <- mean(melted$value)
  allSamplesExceptNoPeptide_res[["allAntigens"]] <- rbind(allSamplesExceptNoPeptide_res[["allAntigens"]], 
                                                          data.frame(`triplet` = triplet, 
                                                                     `p.val` = test_allgroups$p.value, 
                                                                     avg_log2FC = avg_log2FC_allgroups, 
                                                                     `n_obs` = nrow(melted)))
  
  for (i in unique(melted$group)) {
    test <- t.test(melted[melted$group %in% i,"value"], mu = 0)
    avg_log2FC <- mean(melted[melted$group %in% i,"value"])
    allSamplesExceptNoPeptide_res[[i]] <- rbind(allSamplesExceptNoPeptide_res[[i]], 
                                           data.frame(`triplet` = triplet, 
                                                      `p.val` = test$p.value, 
                                                      avg_log2FC = avg_log2FC, 
                                                      `n_obs` = nrow(melted[melted$group %in% i,])))
  }
  
  
  
  
  # ----- Dotplots of all the values ------- #
  dotplot <- ggplot(melted, aes(x=group, y = value, fill=replicate)) + 
    geom_dotplot(binaxis='y', stackdir='center') + ggtitle(triplet) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_continuous(limits = c(-10, 10))
  
  allSamplesExceptNoPeptide_plots[[triplet]] <- dotplot
  
}

# adjust p-values and add in the information of how many data points were tested in each case
for (i in 1:length(allSamplesExceptNoPeptide_res)) {
  allSamplesExceptNoPeptide_res[[i]] <- allSamplesExceptNoPeptide_res[[i]][complete.cases(allSamplesExceptNoPeptide_res[[i]]),]
  allSamplesExceptNoPeptide_res[[i]]$padj <- p.adjust(allSamplesExceptNoPeptide_res[[i]]$p.val, method = "BH")
}


# ---- Save outputs ------------------------------------------
saveRDS(allSamplesExceptNoPeptide_res, file = "outputs/differential_selection_results/allSignificanceResults.Rds")
saveRDS(allSamplesExceptNoPeptide_plots, file = "outputs/differential_selection_results/allEnrichRatioDotPlots.Rds")

for (i in 1:length(allSamplesExceptNoPeptide_res)) {
  write.csv(allSamplesExceptNoPeptide_res[[i]], paste0("outputs/differential_selection_results/",names(allSamplesExceptNoPeptide_res)[i],"_enrichRatio_Results.csv"), row.names = F)
}

write.csv(allEnrichScores, "outputs/enrichment_ratios.csv", row.names = TRUE)
write.csv(logEnrichScores, "outputs/log2_enrichment_ratios.csv", row.names = TRUE)









