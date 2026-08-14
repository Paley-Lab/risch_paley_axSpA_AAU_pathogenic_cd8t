# ==========================================================================
# Merged from 08_plot_TCR_barplots.R and 12_plot_random_TCR_barplots.R
#
# Generates TCR clone percentile tables and barplots ranked by BOTH the
# pathogenic-signature score (Mean_pathogenic_score, from 08) and the
# random-gene-set score (Mean_random_score, from 12), plus the
# TRBV-usage-by-patient plot from 08.
#
# The shared setup / data-loading code and the shared `makeTcrTable()`
# logic (which differed only in which score column was used for sorting
# and percentiles) were consolidated so they run once. The streamlined
# plotting helpers from 12_plot_random_TCR_barplots.R (make_decile,
# summarize_counts, prep_and_order, plot_decile_bar) are reused for both
# the pathogenic-score plots and the random-score plots.
# ==========================================================================

# run this on RIS with the following command:
# bsub -Is -M 64GB -R 'rusage[mem=64GB]' -G compute-mgriffit -q general-interactive -a 'docker(igrisch/seurat:5.2.0_v2)' /bin/bash
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP021")
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(patchwork)

# loading in our data and making sure it's leveled the way we desire
d <- readRDS("data/robjects/4datasets_integrated.Rds")
d$B27 <- factor(as.character(d$B27), levels = c("B27-", "B27+"))
d$tissue_consensus <- factor(as.character(d$tissue_consensus), levels=c("PBMC", "Aqueous","Synovial_Fluid"))
# d$disease_merged <- ifelse(!is.na(d$disease_merged), d$disease_merged, 
#                               ifelse(!is.na(d$Anatomy), d$Anatomy, 
#                                      "PsA"))
d$disease_merged <- factor(as.character(d$disease_merged), levels=c("HC", "PanU", "AS", "AU", "PsA"))
d$subject <- factor(as.character(d$subject), levels=c("HC1564", "HC1660", "HC1788", "HC546", "73-4", "76", "81", "74", "SpA-08", "SpA-07", "SpA-01", "SpA-03", "SpA-04", "UV019",
                                                      "UV027", "UV122", "UV180", "PSA1", "PSA2", "PSA3", "PSA4"))


# save old pathogenic sequence and update to the new one (Jul 12, 2026)
d$pathogenic_old <- d$pathogenic
newPathogenicSeq <- read.csv("../pathogenic_motifs", header = F)
d$pathogenic <- ifelse(grepl(newPathogenicSeq[1,2], sapply(strsplit(d$CTaa, "_"), "[", 2)), TRUE, FALSE)
d$pathogenic <- ifelse(grepl("TRAV21", d$CTstrict), d$pathogenic, FALSE)
d$pathogenic <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", d$CTstrict), d$pathogenic, FALSE)


# source our palettes to keep everything consistent
source("../color_palettes.R")

###########################################################################
########## Getting a dataframe of TCR sequence frequencies, NOT split by tissue ###########
###########################################################################

## FUNCTION to get the exact table we want, with each clone's average pathogenic
## AND average random-signature score, plus a quartile/percentile ranking based
## on whichever score `score_col` points to ("Mean_pathogenic_score" or
## "Mean_random_score").
makeTcrTable <- function(meta, score_col = "Mean_pathogenic_score") {
  # get the mean pathogenic + random signature score for each clone
  means <- meta %>%
    group_by(CTstrict, subject) %>%
    dplyr::summarize(Mean_pathogenic_score = mean(Core_JCI_MP156_c4_vs_c3_noTCR_20g, na.rm=TRUE),
                     Mean_random_score = mean(Random_20g, na.rm=TRUE))
  # add the mean signature score data to the metadata
  meta <- merge(meta, means, by = c("CTstrict", "subject"))
  meta_sub <- meta
  # #subset to the non-MAIT cells
  # meta_sub <- meta_sub[meta_sub$celltype1!="MAIT",]
  #subset to only those cells with clonality > 1
  meta_sub <- meta_sub[meta_sub$clonalFrequency > 1, ]
  #remove cells with no beta chain
  meta_sub <- meta_sub[!grepl("_NA", meta_sub$CTaa),]
  #remove cells with no alpha chain
  meta_sub <- meta_sub[!grepl("^NA_", meta_sub$CTaa),]
  #remove cells with multiple chains
  meta_sub <- meta_sub[!grepl("\\;", meta_sub$CTgene),]
  # subset to only the columns we're interested in and unique
  meta_sub_unique <- meta_sub[,c("subject", "pathogenic", "CTgene", 
                                 "CTnt", "CTaa", "CTstrict", "clonalFrequency", "cloneSize", 
                                 "disease_merged", "Mean_pathogenic_score", "Mean_random_score", "dataset")] %>% unique()
  
  #sort by descending score_col
  meta_sub_unique <- meta_sub_unique[order(meta_sub_unique[[score_col]], decreasing = T),]
  meta_sub_unique <- meta_sub_unique[complete.cases(meta_sub_unique),]
  
  # add quartile/percentile info, based on score_col
  summary_info <- summary(meta_sub_unique[[score_col]])
  meta_sub_unique$quartile <- ifelse(meta_sub_unique[[score_col]] >= summary_info[5], "Q4",
                                     ifelse(meta_sub_unique[[score_col]] >= summary_info[3], "Q3",
                                            ifelse(meta_sub_unique[[score_col]] >= summary_info[2], "Q2",
                                                   "Q1")))
  meta_sub_unique$percentile <- cut(meta_sub_unique[[score_col]], 
                                    breaks = quantile(meta_sub_unique[[score_col]], 0:20 / 20), 
                                    labels=1:20/20, include.lowest=TRUE)
  
  #sort by descending score_col
  meta_sub_unique$TCRa_AA <- sapply(strsplit(meta_sub_unique$CTaa, "_"), "[", 1)
  meta_sub_unique$TCRb_AA <- sapply(strsplit(meta_sub_unique$CTaa, "_"), "[", 2)
  meta_sub_unique <- tidyr::separate(meta_sub_unique, col="CTgene", sep="\\.|_", into=c("TRAV", "TRAJ", NA, "TRBV", NA, "TRBJ", "TRBC"), remove = F)
  
  return(meta_sub_unique)
}

meta_full <- d@meta.data
write.table(meta_full, file = "data/4datasets_integrated_withSig_meta.tsv", sep="\t", row.names = T, col.names = T)
meta_full <- read.delim("data/4datasets_integrated_withSig_meta.tsv")

# non-MAIT cells (base population for all four percentile tables below)
tmp <- subset(d, subset = celltype1 != "MAIT")
meta_nonMAIT <- tmp@meta.data
# non-MAIT, TRAV21+ cells only
meta_nonMAIT_TRAV21pos <- meta_nonMAIT[grep("TRAV21", meta_nonMAIT$CTgene),]

## -- Pathogenic-score-ranked tables (from 08_plot_TCR_barplots.R) --
meta_sub_unique_pathogenic_nonMAIT <- makeTcrTable(meta_nonMAIT, score_col = "Mean_pathogenic_score")
write.table(meta_sub_unique_pathogenic_nonMAIT, file = "data/TCR_proportions_and_percentiles_nonMAIT.txt", sep="\t", row.names = T, col.names = T)

meta_sub_unique_pathogenic_TRAV21pos <- makeTcrTable(meta_nonMAIT_TRAV21pos, score_col = "Mean_pathogenic_score")
write.table(meta_sub_unique_pathogenic_TRAV21pos, file = "data/TCR_proportions_and_percentiles_nonMAIT_TRAV21pos.txt", sep="\t", row.names = T, col.names = T)

## -- Random-score-ranked tables (from 12_plot_random_TCR_barplots.R) --
meta_sub_unique_random_nonMAIT <- makeTcrTable(meta_nonMAIT, score_col = "Mean_random_score")
write.table(meta_sub_unique_random_nonMAIT, file = "data/TCR_proportions_and_percentiles_nonMAIT_randomScore.txt", sep="\t", row.names = T, col.names = T)

meta_sub_unique_random_TRAV21pos <- makeTcrTable(meta_nonMAIT_TRAV21pos, score_col = "Mean_random_score")
write.table(meta_sub_unique_random_TRAV21pos, file = "data/TCR_proportions_and_percentiles_nonMAIT_TRAV21pos_randomScore.txt", sep="\t", row.names = T, col.names = T)


###########################################################################
########## Barplots of TRBV / CDR3b usage by decile (streamlined helpers) #########
###########################################################################

## get a whole bunch of colors to use for the TRBVs.
## (Originally generated via brewer.pal() + sample(); hardcoded here for
## reproducibility -- the random draw in the source scripts was dead code,
## since it was always overwritten by this fixed vector.)
palette <- c("#C7E9C0", "#252525", "#7BCCC4", "#006D2C", "#084081", "#74C476",
             "#C6DBEF", "#2171B5", "#8C96C6", "#0868AC",
             "#DEEBF7", "#CCEBC5", "#A1D99B", "#525252", "#E0F3DB", "#BFD3E6",
             "#08519C", "#9EBCDA", "#238B45", "#810F7C",
             "#F7FCF5", "#6BAED6", "#F0F0F0", "#D9D9D9", "#41AB5D", "#737373",
             "#8C6BB1", "#E5F5E0", "#4292C6", "#A8DDB5", "#4EB3D3", "#00441B",
             "#88419D", "#969696", "#F7FCFD", "#9ECAE1", "#2B8CBE",
             "#4D004B", "#08306B", "#BDBDBD")

## ---------------------------------------------------------------------
## Helper: collapse the 5%-resolution `percentile` column into 10 decile
## bins (width = 0.1 = 10%), labeled "0-10%", "10-20%", ..., "90-100%".
## Because the source values are exact multiples of 0.05, pairing them up
## with ceiling(x / 0.1) * 0.1 is exact -- no floating point edge cases.
## ---------------------------------------------------------------------
make_decile <- function(percentile_vec) {
  decile_num <- ceiling(round(percentile_vec, 10) / 0.1) * 0.1
  decile_num <- round(decile_num, 1)  # guard against floating point drift
  decile_lab <- paste0(round((decile_num - 0.1) * 100), "-", round(decile_num * 100), "%")
  levels_ordered <- paste0(seq(0, 90, 10), "-", seq(10, 100, 10), "%")
  factor(decile_lab, levels = levels_ordered)
}

## ---- Helper: unweighted count of `group_col` per decile ----
summarize_counts <- function(data, group_col) {
  tab <- as.data.frame(table(data[[group_col]], data$Decile))
  colnames(tab) <- c("Group", "Decile", "Frequency")
  tab
}

## ---- Helper: order by descending frequency + apply sequential releveling ----
## (each successive level in `relevel_to` is moved to the front, so the
## LAST element ends up as the first factor level -- matches the original
## chained relevel() %>% relevel() calls)
prep_and_order <- function(df, relevel_to = NULL) {
  df <- df[order(df$Frequency, decreasing = TRUE), ]
  df$Group <- factor(df$Group, levels = unique(df$Group))
  for (lvl in relevel_to) df$Group <- relevel(df$Group, lvl)
  df
}

## ---- Helper: stacked-proportion decile barplot ----
plot_decile_bar <- function(df, fill_values, scale_y_reverse = TRUE) {
  p <- ggplot(df, aes(x = Decile, y = Frequency, fill = Group)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(values = fill_values)
  if (scale_y_reverse) p <- p + scale_y_reverse()
  return(p)
}

###########################################################################
## PATHOGENIC-SCORE plots (from 08_plot_TCR_barplots.R)
###########################################################################
meta_pathogenic_TRAV21pos <- read.delim("data/TCR_proportions_and_percentiles_nonMAIT_TRAV21pos.txt")
meta_pathogenic_TRAV21pos$Decile <- make_decile(meta_pathogenic_TRAV21pos$percentile)

meta_pathogenic_nonMAIT <- read.delim("data/TCR_proportions_and_percentiles_nonMAIT.txt")
meta_pathogenic_nonMAIT$Decile <- make_decile(meta_pathogenic_nonMAIT$percentile)

## 1. TRBV usage, TRAV21+ cells -- multicolor
df1 <- summarize_counts(meta_pathogenic_TRAV21pos, "TRBV") %>% prep_and_order(relevel_to = "TRBV5-5")
# write.csv(df1, "plot/supplementary/barplot_trbv_usage_by_percentile_nonMAIT_TRAV21pos_multicolor_summary_stats.csv", row.names = FALSE)
p1 <- plot_decile_bar(df1, c("#FB9A99", "hotpink", palette))
ggsave(plot = p1, filename = "plot/supplementary/barplot_trbv_usage_by_percentile_nonMAIT_TRAV21pos_multicolor.pdf",
       device = "pdf", width = 5, height = 2)

## 2. TRBV usage, TRAV21+ cells -- non-multicolor (easier to see the point in the main figures)
df2 <- summarize_counts(meta_pathogenic_TRAV21pos, "TRBV") %>% prep_and_order(relevel_to = c("TRBV9", "TRBV5-5"))
write.csv(df2, "plot/trbv_usage_barplots/barplot_trbv_usage_by_percentile_nonMAIT_TRAV21pos_summary_stats.csv", row.names = FALSE)
p2 <- plot_decile_bar(df2, c("#FB9A99", "hotpink", rep("lightgrey", 40)))
ggsave(plot = p2, filename = "plot/trbv_usage_barplots/barplot_trbv_usage_by_percentile_nonMAIT_TRAV21pos.pdf",
       device = "pdf", width = 5, height = 2)

## 3. TRBV usage, non-MAIT cells -- non-multicolor (proportion)
df3 <- summarize_counts(meta_pathogenic_nonMAIT, "TRBV") %>% prep_and_order(relevel_to = c("TRBV9", "TRBV5-5"))
write.csv(df3, "plot/trbv_usage_barplots/barplot_trbv_usage_by_percentile_nonMAIT_summary_stats.csv", row.names = FALSE)
p3 <- plot_decile_bar(df3, c("#FB9A99", "hotpink", rep("lightgrey", 43)))
ggsave(plot = p3, filename = "plot/trbv_usage_barplots/barplot_trbv_usage_by_percentile_nonMAIT.pdf",
       device = "pdf", width = 5, height = 2)

## 4. Pathogenic CDR3b motif by percentile -- TRAV21+ cells
df4 <- summarize_counts(meta_pathogenic_TRAV21pos, "pathogenic") %>% prep_and_order()
write.csv(df4, "plot/cdr3b_barplots/barplot_CDR3b_by_percentile_nonMAIT_TRAV21pos.csv", row.names = FALSE)
p4 <- plot_decile_bar(df4, c("lightgrey", "maroon"), scale_y_reverse = FALSE)
ggsave(plot = p4, filename = "plot/cdr3b_barplots/barplot_CDR3b_by_percentile_nonMAIT_TRAV21pos.pdf",
       device = "pdf", width = 3.3, height = 2)

## 5. Pathogenic CDR3b motif by percentile -- non-MAIT cells
df5 <- summarize_counts(meta_pathogenic_nonMAIT, "pathogenic") %>% prep_and_order()
write.csv(df5, "plot/cdr3b_barplots/barplot_CDR3b_by_percentile_nonMAIT.csv", row.names = FALSE)
p5 <- plot_decile_bar(df5, c("lightgrey", "maroon"), scale_y_reverse = FALSE)
ggsave(plot = p5, filename = "plot/cdr3b_barplots/barplot_CDR3b_by_percentile_nonMAIT.pdf",
       device = "pdf", width = 3.3, height = 2)


###########################################################################
## RANDOM-SCORE plots (from 12_plot_random_TCR_barplots.R)
###########################################################################
meta_random_TRAV21pos <- read.delim("data/TCR_proportions_and_percentiles_nonMAIT_TRAV21pos_randomScore.txt")
meta_random_TRAV21pos$Decile <- make_decile(meta_random_TRAV21pos$percentile)

meta_random_nonMAIT <- read.delim("data/TCR_proportions_and_percentiles_nonMAIT_randomScore.txt")
meta_random_nonMAIT$Decile <- make_decile(meta_random_nonMAIT$percentile)

## 1. TRBV usage, TRAV21+ cells -- multicolor (excluded from combined fig)
rdf1 <- summarize_counts(meta_random_TRAV21pos, "TRBV") %>% prep_and_order(relevel_to = "TRBV5-5")
write.csv(rdf1, "plot/supplementary/barplot_RANDOMSCORE_trbv_usage_by_percentile_nonMAIT_TRAV21pos_multicolor_summary_stats.csv", row.names = FALSE)
rp1 <- plot_decile_bar(rdf1, c("#FB9A99", "hotpink", palette))
ggsave(plot = rp1, filename = "plot/supplementary/barplot_RANDOMSCORE_trbv_usage_by_percentile_nonMAIT_TRAV21pos_multicolor.pdf",
       device = "pdf", width = 5, height = 2)

## 2. TRBV usage, TRAV21+ cells -- non-multicolor
rdf2 <- summarize_counts(meta_random_TRAV21pos, "TRBV") %>% prep_and_order(relevel_to = c("TRBV9", "TRBV5-5"))
write.csv(rdf2, "plot/trbv_usage_barplots/barplot_RANDOMSCORE_trbv_usage_by_percentile_nonMAIT_TRAV21pos_summary_stats.csv", row.names = FALSE)
rp2 <- plot_decile_bar(rdf2, c("#FB9A99", "hotpink", rep("lightgrey", 40)))
ggsave(plot = rp2, filename = "plot/trbv_usage_barplots/barplot_RANDOMSCORE_trbv_usage_by_percentile_nonMAIT_TRAV21pos.pdf",
       device = "pdf", width = 5, height = 2)

## 3. TRBV usage, non-MAIT cells -- non-multicolor
rdf3 <- summarize_counts(meta_random_nonMAIT, "TRBV") %>% prep_and_order(relevel_to = c("TRBV9", "TRBV5-5"))
write.csv(rdf3, "plot/trbv_usage_barplots/barplot_RANDOMSCORE_trbv_usage_by_percentile_nonMAIT_summary_stats.csv", row.names = FALSE)
rp3 <- plot_decile_bar(rdf3, c("#FB9A99", "hotpink", rep("lightgrey", 43)))
ggsave(plot = rp3, filename = "plot/trbv_usage_barplots/barplot_RANDOMSCORE_trbv_usage_by_percentile_nonMAIT.pdf",
       device = "pdf", width = 5, height = 2)

## 4. Pathogenic CDR3b motif by percentile -- TRAV21+ cells
rdf4 <- summarize_counts(meta_random_TRAV21pos, "pathogenic") %>% prep_and_order()
write.csv(rdf4, "plot/barplot_RANDOMSCORE_CDR3b_by_percentile_nonMAIT_TRAV21pos.csv", row.names = FALSE)
rp4 <- plot_decile_bar(rdf4, c("lightgrey", "maroon"), scale_y_reverse = FALSE)
ggsave(plot = rp4, filename = "plot/barplot_RANDOMSCORE_CDR3b_by_percentile_nonMAIT_TRAV21pos.pdf",
       device = "pdf", width = 3.3, height = 2)

## 5. Pathogenic CDR3b motif by percentile -- non-MAIT cells
rdf5 <- summarize_counts(meta_random_nonMAIT, "pathogenic") %>% prep_and_order()
write.csv(rdf5, "plot/barplot_RANDOMSCORE_CDR3b_by_percentile_nonMAIT.csv", row.names = FALSE)
rp5 <- plot_decile_bar(rdf5, c("lightgrey", "maroon"), scale_y_reverse = FALSE)
ggsave(plot = rp5, filename = "plot/barplot_RANDOMSCORE_CDR3b_by_percentile_nonMAIT.pdf",
       device = "pdf", width = 3.3, height = 2)

## 6. Combine the 4 non-multicolor random-score barplots (rp2-rp5) into a
##    single row, no legends, saved as one PDF
combined_random <- (rp5 + rp4 + rp3 + rp2) +
  patchwork::plot_layout(nrow = 1) &
  theme(legend.position = "none")

ggsave(plot = combined_random, filename = "plot/cdr3b_barplots/combined_RANDOMSCORE_barplots_non_multicolor_single_row.pdf",
       device = "pdf", width = 8, height = 2)


##############################################################################################
########## Plot vbeta usage by person to demonstrate the diversity of vbetas used ###########
##############################################################################################
# # NOTE: in the original 08_plot_TCR_barplots.R, this section re-used the
# # `meta` variable that was left over from building the TRAV21+ percentile
# # table above (via `meta <- meta[grep("TRAV21", meta$CTgene),]`), so it
# # actually plots TRAV21+, non-MAIT cells rather than literally "all cells"
# # as the output filename suggests. That original behavior is preserved
# # here explicitly via `meta_nonMAIT_TRAV21pos` -- flagging this in case it
# # was unintentional in the source script and you'd rather plot on
# # `meta_nonMAIT` (all non-MAIT cells) instead.
# vbetaUse <- meta_nonMAIT_TRAV21pos %>%
#   group_by(CTgene, subject) %>%
#   dplyr::summarize(total = n())
# vbetaUse <- vbetaUse[complete.cases(vbetaUse),]
# vbetaUse$CTgene <- sapply(strsplit(vbetaUse$CTgene, "_"), "[", 2) # remove the TCRa sequences. We only care about vbeta here and it's possible to have 2 alpha chains, which throws off the following steps
# vbetaUse <- vbetaUse[!grepl(";",vbetaUse$CTgene),] # remove any cells with two TCRb chains, this isn't really a thing
# vbetaUse <- vbetaUse[!grepl("NA", vbetaUse$CTgene),] # remove cells without a beta chain
# vbetaUse <- tidyr::separate(vbetaUse, col="CTgene", sep="\\.", into=c("TRBV", NA, "TRBJ", "TRBC"), remove = F)
# 
# # get vbeta usage if we weight each CLONE equally
# flatVbeta <- vbetaUse %>%
#   group_by(TRBV, subject) %>%
#   dplyr::summarize(total = n())
# flatVbeta$TRBV <- factor(flatVbeta$TRBV, levels=unique(flatVbeta$TRBV)) %>% 
#   relevel(c("TRBV9")) %>% 
#   relevel(c("TRBV5-5"))
# p_vbeta <- ggplot(flatVbeta, aes(x=subject, y=total, fill=TRBV)) +
#   geom_bar(stat="identity", position = "fill")+theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5, size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   scale_fill_manual(values = c("#FB9A99", "hotpink", palette, "royalblue", "navy", "#8A2BE2", "#9683EC")) 
# ggsave(plot=p_vbeta, filename = "plot/barplot_trbv_usage_by_patient_all_cells.pdf", device="pdf", width=8, height=3)
