## =============================================================================
## TCR-seq Technology Comparisons
## Consolidated from:
##   01_TIRTL_omniscope_depth_compare.R
##   02_plot_TIRTL_vs_omniscope_depth_comparison.R
##
## This file contains TWO INDEPENDENT analyses, kept in clearly separated
## sections (they compare different technologies on different metrics and
## do not share any data):
##
##   PART 1: CAIRR-seq vs TIRTL-seq -- total read/UMI depth, TRA vs TRB
##            (data source: "Omniscope" directory, which is the vendor/folder
##            name for the CAIRR-seq assay -- same technology as Part 2)
##            (data prep from script 01 + stats/plot from script 02)
##   PART 2: CAIRR-seq vs TIRTL-seq -- unique sequence counts, Alpha/Beta/Paired
##            (originally the second half of script 01)
##
## IMPORTANT: Omniscope and CAIRR-seq are the SAME technology, referred to by
## two different names in the original scripts (Part 1 called it "Omniscope",
## Part 2 called it "CAIRR-seq"). Labels and colors below have been harmonized
## to "CAIRR-seq" / "TIRTL-seq" throughout so both parts are visually and
## terminologically consistent. Part 1 and Part 2 remain two SEPARATE plots
## (one per metric: total reads vs. unique sequences) rather than being
## merged into a single combined analysis.
##
## Fix applied: script 01 read "all_alpha_data.csv" into BOTH
## all_samples_alpha and all_samples_beta (copy-paste bug), so every "TRB"
## row in the original read_depth_omni_vs_tirtl.csv was actually alpha-chain
## data mislabeled as beta. This version reads "all_beta_data.csv" for the
## beta variable instead.
## =============================================================================

## ---- 0. Setup ------------------------------------------------------------------
library(tidyverse)
library(scales)

## #############################################################################
## PART 1: CAIRR-seq vs TIRTL-seq -- total read depth comparison (TRA vs TRB)
## #############################################################################

## ---- 1.1 Load & aggregate Omniscope data ----------------------------------------
omni_datadir <- "/storage1/fs1/paleym/Active/SeqData/Omniscope/"
omni_contig_files <- list.files(omni_datadir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
omni_raw <- lapply(omni_contig_files, read.csv)

# total reads per sample, split by locus (alpha/beta)
omni_total_per_sample <- lapply(omni_raw, function(df) {
  df %>%
    group_by(locus) %>%
    summarise(totalCounts = sum(umi_count), .groups = "drop")
})
for (i in 1:length(omni_total_per_sample)) {
  omni_total_per_sample[[i]]$sample <- sapply(strsplit(omni_contig_files[i], "//|.csv"), "[", 2)
}
omni_bound <- bind_rows(omni_total_per_sample)
omni_bound$tech <- "CAIRR-seq"   # HARMONIZED: was "omniscope" -- same assay as Part 2's CAIRR-seq

## ---- 1.2 Load & aggregate TIRTL data ---------------------------------------------
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP028")

all_samples_alpha <- read.csv("data/outputs/all_alpha_data.csv")
all_samples_beta  <- read.csv("data/outputs/all_beta_data.csv")   

tirtl_total_alpha_per_sample <- all_samples_alpha %>%
  group_by(sample) %>%
  summarise(totalCounts = sum(readCount), .groups = "drop") %>%
  mutate(tech = "TIRTL-seq", locus = "TRA") %>%
  as.data.frame()

tirtl_total_beta_per_sample <- all_samples_beta %>%
  group_by(sample) %>%
  summarise(totalCounts = sum(readCount), .groups = "drop") %>%
  mutate(tech = "TIRTL-seq", locus = "TRB") %>%
  as.data.frame()

## ---- 1.3 Combine & export raw depth data -----------------------------------------
depth_bound <- bind_rows(
  omni_bound[, c("sample", "locus", "totalCounts", "tech")],
  tirtl_total_alpha_per_sample[, c("sample", "locus", "totalCounts", "tech")],
  tirtl_total_beta_per_sample[, c("sample", "locus", "totalCounts", "tech")]
)

depth_out_path <- "/storage1/fs1/paleym/Active/Isabel_Risch/MP036/data/umi_counts_omni_vs_tirtl.csv"
write.csv(depth_bound, file = depth_out_path, row.names = FALSE)

# deidentify the Omniscope data and then re-save the file for source data
keys <- read.csv("/storage1/fs1/paleym/Active/Isabel_Risch/dataset_uploading_for_publication/2026_risch_pathogenic_CD8T/deidentifying_patients_keys/keys.csv", sep = "\t")
keys <- keys[keys$Dataset=="Omniscope ",]
depth_bound_deid <- depth_bound; depth_bound_deid$sample <- toupper(depth_bound_deid$sample)
depth_bound_deid <- merge(keys, depth_bound_deid, by.x = "Original.Key", by.y = "sample", all.y = T)
depth_bound_deid$sample <- ifelse(is.na(depth_bound_deid$Publicly.Uploaded.Key), depth_bound_deid$Original.Key, depth_bound_deid$Publicly.Uploaded.Key)
depth_bound_deid <- depth_bound_deid[, c("sample", "locus", "totalCounts", "tech")]
depth_out_path <- "/storage1/fs1/paleym/Active/Isabel_Risch/MP036/data/umi_counts_omni_vs_tirtl_DEIDENTIFIED.csv"
write.csv(depth_bound_deid, file = depth_out_path, row.names = FALSE)

## ---- 1.4 Tidy for plotting --------------------------------------------------------
depth_d <- depth_bound %>%
  mutate(
    tech  = factor(tech, levels = c("CAIRR-seq", "TIRTL-seq")),
    locus = factor(locus, levels = c("TRA", "TRB"))
  )

## ---- 1.5 Summary statistics (log10 space) -----------------------------------------
## Geometric mean +/- SD, i.e. mean/SD of log10(totalCounts) then exponentiated
## for display -- appropriate since read counts span multiple orders of magnitude.
depth_summary_df <- depth_d %>%
  group_by(tech, locus) %>%
  summarise(
    n           = n(),
    mean        = mean(totalCounts),
    mean_log    = mean(log10(totalCounts)),
    sd_log      = sd(log10(totalCounts)),
    geo_mean    = 10^mean_log,
    geo_sd_low  = 10^(mean_log - sd_log),
    geo_sd_high = 10^(mean_log + sd_log),
    .groups = "drop"
  )

## ---- 1.6 Statistical testing (technology comparison, per locus) -------------------
## Wilcoxon rank-sum (Mann-Whitney): CAIRR-seq vs TIRTL-seq, within each locus
depth_stat_tests <- depth_d %>%
  group_by(locus) %>%
  summarise(p_value = wilcox.test(totalCounts ~ tech, exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(stars = ifelse(p_value < 0.05, "*", "ns"))

## ---- 1.7 Bracket geometry ---------------------------------------------------------
## Positions computed from the data itself (rather than the hard-coded y-values
## in the original script 02) so the brackets stay correctly placed regardless
## of the actual read-depth range.
depth_dodge_width <- 0.75
depth_bar_offset  <- 0.19   # approx. half-distance between the two dodged bars
depth_locus_levels <- levels(depth_d$locus)

depth_bracket_geom <- map_dfr(seq_along(depth_locus_levels), function(i) {
  lv <- depth_locus_levels[i]
  local_max <- max(
    depth_d$totalCounts[depth_d$locus == lv],
    depth_summary_df$geo_sd_high[depth_summary_df$locus == lv]
  )
  log_max   <- log10(local_max)
  bracket_y <- log_max + 0.12
  label_y   <- log_max + 0.22
  tibble(locus = lv, x = i - depth_bar_offset, xend = i + depth_bar_offset,
         y = bracket_y, yend = bracket_y, label_x = i, label_y = label_y)
}) %>%
  left_join(depth_stat_tests, by = "locus")

depth_exp_low    <- floor(log10(min(depth_d$totalCounts)))
depth_exp_high   <- ceiling(max(depth_bracket_geom$label_y))
depth_exp_range  <- depth_exp_low:depth_exp_high
depth_upper_limit <- 10^depth_exp_high

## ---- 1.8 Plot -----------------------------------------------------------------------
depth_p <- ggplot(depth_d, aes(x = locus, y = totalCounts, fill = tech)) +
  geom_col(data = depth_summary_df, aes(y = geo_mean), position = position_dodge(depth_dodge_width),
           width = 0.7, color = "grey20") +
  geom_errorbar(data = depth_summary_df,
                aes(y = geo_mean, ymin = geo_sd_low, ymax = geo_sd_high),
                position = position_dodge(depth_dodge_width), width = 0.15) +
  geom_point(aes(color = tech),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = depth_dodge_width),
             shape = 21, fill = "white", size = 2, show.legend = FALSE) +
  ggplot2::annotate("segment", x = depth_bracket_geom$x, xend = depth_bracket_geom$xend,
                    y = 10^depth_bracket_geom$y, yend = 10^depth_bracket_geom$yend) +
  ggplot2::annotate("text", x = depth_bracket_geom$label_x, y = 10^depth_bracket_geom$label_y,
                    label = depth_bracket_geom$stars, size = 4.5, fontface = "bold") +
  scale_y_log10(breaks = 10^depth_exp_range,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks = NULL,
                limits = c(NA, depth_upper_limit)) +
  scale_fill_manual(values = c("CAIRR-seq" = "#3B7DB8", "TIRTL-seq" = "#D9782D")) +
  scale_color_manual(values = c("CAIRR-seq" = "grey20", "TIRTL-seq" = "grey20")) +
  labs(title = "TCR Read Depth: TRA vs TRB by Technology",
       subtitle = "CAIRR-seq (Omniscope assay) vs TIRTL-seq",
       x = NULL, y = "Total TCR Read Counts", fill = "Technology",
       caption = "Bars: geometric mean \u00b1 SD (log scale). Brackets: Wilcoxon rank-sum test\n(* p<0.05, ns = not significant)"
  ) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.minor.ticks.length = unit(0.075, "cm"),
        axis.ticks = element_line(color = "grey30"))

print(depth_p)

## ---- 1.9 Save outputs ----------------------------------------------------------------
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP036")
dir.create("data/outputs/plots", showWarnings = FALSE, recursive = TRUE)

ggsave("data/outputs/plots/omni_vs_tirtl_read_depth_summary.png", depth_p, width = 5, height = 4, dpi = 600)
ggsave("data/outputs/plots/omni_vs_tirtl_read_depth_summary.pdf", depth_p, width = 5, height = 4)

depth_summary_out <- depth_summary_df %>%
  left_join(depth_stat_tests, by = "locus")
write.csv(depth_summary_out, file = "data/outputs/summary_omni_vs_tirtl.csv", row.names = FALSE)

cat("PART 1 done: saved omni_vs_tirtl_read_depth_summary.png/.pdf and summary_omni_vs_tirtl.csv\n")
cat("Wilcoxon rank-sum p-values (CAIRR-seq vs TIRTL-seq):\n")
print(depth_stat_tests)


## #############################################################################
## PART 2: CAIRR-seq vs TIRTL-seq -- unique sequence count comparison
## (independent of Part 1: different technologies, different metric --
##  unique clones rather than total reads)
## #############################################################################

uniq_out_dir <- "/storage1/fs1/paleym/Active/Isabel_Risch/MP036/plots"
dir.create(uniq_out_dir, showWarnings = FALSE)

## ---- 2.1 Load & tidy data ---------------------------------------------------------
uniq_cairr <- read.csv("/storage1/fs1/paleym/Active/Isabel_Risch/MP002/data/CAIRRseq_unique_sequence_counts_per_sample.csv") %>%
  mutate(technology = "CAIRR-seq")
uniq_tirtl <- read.csv("/storage1/fs1/paleym/Active/Isabel_Risch/MP028/data/TIRTLseq_unique_sequence_counts_per_sample.csv") %>%
  mutate(technology = "TIRTL-seq")

uniq_combined_wide <- bind_rows(uniq_tirtl, uniq_cairr) %>%
  mutate(tech = factor(technology, levels = c("CAIRR-seq", "TIRTL-seq")))

uniq_d <- uniq_combined_wide %>%
  pivot_longer(
    cols = c(numberUniqueAlphas, numberUniqueBetas, numberUniquePaired),
    names_to = "locus", values_to = "totalCounts"
  ) %>%
  mutate(
    locus = recode(locus,
                   numberUniqueAlphas = "Alpha",
                   numberUniqueBetas  = "Beta",
                   numberUniquePaired = "Paired"),
    locus = factor(locus, levels = c("Alpha", "Beta", "Paired"))
  )

# deidentify the Omniscope data and then re-save the file for source data
keys <- read.csv("/storage1/fs1/paleym/Active/Isabel_Risch/dataset_uploading_for_publication/2026_risch_pathogenic_CD8T/deidentifying_patients_keys/keys.csv", sep = "\t")
keys <- keys[keys$Dataset=="Omniscope ",]
uniq_d_deid <- uniq_d; uniq_d$sample <- toupper(uniq_d$sample)
uniq_d_deid <- merge(keys, uniq_d_deid, by.x = "Original.Key", by.y = "sample", all.y = T)
uniq_d_deid$sample <- ifelse(is.na(uniq_d_deid$Publicly.Uploaded.Key), uniq_d_deid$Original.Key, uniq_d_deid$Publicly.Uploaded.Key)
uniq_d_deid <- uniq_d_deid[, c("sample", "locus", "totalCounts", "tech")]

write_csv(uniq_d_deid, file.path(uniq_out_dir, "source_data_per_sample_long.csv"))
write_csv(uniq_combined_wide, file.path(uniq_out_dir, "source_data_per_sample_wide.csv"))

## ---- 2.2 Summary statistics (log10 space) -----------------------------------------
uniq_summary_df <- uniq_d %>%
  group_by(locus, tech) %>%
  summarise(
    n           = n(),
    mean        = mean(totalCounts),
    mean_log    = mean(log10(totalCounts)),
    sd_log      = sd(log10(totalCounts)),
    geo_mean    = 10^mean_log,
    geo_sd    = 10^sd_log,
    geo_sd_low  = 10^(mean_log - sd_log),
    geo_sd_high = 10^(mean_log + sd_log),
    .groups = "drop"
  )

write_csv(uniq_summary_df, file.path(uniq_out_dir, "alpha_beta_clonewise_summary_stats.csv"))

## ---- 2.3 Statistical testing (technology comparison, per locus) -------------------
uniq_stat_tests <- uniq_d %>%
  group_by(locus) %>%
  summarise(p_value = wilcox.test(totalCounts ~ tech)$p.value, .groups = "drop") %>%
  mutate(stars = ifelse(p_value < 0.05, "*", "ns"))

write_csv(uniq_stat_tests, file.path(uniq_out_dir, "source_data_wilcox.csv"))

print(uniq_summary_df)
print(uniq_stat_tests)

## ---- 2.4 Bracket geometry ------------------------------------------------------
uniq_dodge_width <- 0.75
uniq_bar_offset  <- 0.19
uniq_locus_levels <- levels(uniq_d$locus)

uniq_bracket_geom <- map_dfr(seq_along(uniq_locus_levels), function(i) {
  lv <- uniq_locus_levels[i]
  local_max <- max(
    uniq_d$totalCounts[uniq_d$locus == lv],
    uniq_summary_df$geo_sd_high[uniq_summary_df$locus == lv]
  )
  log_max   <- log10(local_max)
  bracket_y <- log_max + 0.12
  label_y   <- log_max + 0.22
  tibble(locus = lv, x = i - uniq_bar_offset, xend = i + uniq_bar_offset,
         y = bracket_y, yend = bracket_y, label_x = i, label_y = label_y)
}) %>%
  left_join(uniq_stat_tests, by = "locus")

uniq_exp_low    <- 0
uniq_exp_high   <- ceiling(max(uniq_bracket_geom$label_y))
uniq_exp_range  <- uniq_exp_low:uniq_exp_high
uniq_upper_limit <- 10^uniq_exp_high

## ---- 2.5 Plot --------------------------------------------------------------------
## NOTE: the significance-bracket annotate() calls were commented out in the
## original script (01_TIRTL_omniscope_depth_compare.R, lines 179-182), so
## this plot currently does NOT show brackets, matching that script's actual
## behavior. Uncomment the two ggplot2::annotate() calls below to re-enable them.
uniq_p <- ggplot(uniq_d, aes(x = locus, y = totalCounts, fill = tech)) +
  geom_col(data = uniq_summary_df, aes(y = geo_mean), position = position_dodge(uniq_dodge_width),
           width = 0.7, color = "grey20") +
  geom_errorbar(data = uniq_summary_df,
                aes(y = geo_mean, ymin = geo_sd_low, ymax = geo_sd_high),
                position = position_dodge(uniq_dodge_width), width = 0.15) +
  geom_point(aes(color = tech),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = uniq_dodge_width),
             shape = 21, fill = "white", size = 2, show.legend = FALSE) +
  # ggplot2::annotate("segment", x = uniq_bracket_geom$x, xend = uniq_bracket_geom$xend,
  #                   y = 10^uniq_bracket_geom$y, yend = 10^uniq_bracket_geom$yend) +
  # ggplot2::annotate("text", x = uniq_bracket_geom$label_x, y = 10^uniq_bracket_geom$label_y,
  #                   label = uniq_bracket_geom$stars, size = 4.5, fontface = "bold") +
  scale_y_log10(breaks = 10^uniq_exp_range,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks = NULL,
                limits = c(NA, uniq_upper_limit)) +
  scale_fill_manual(values = c("TIRTL-seq" = "#D9782D", "CAIRR-seq" = "#3B7DB8")) +
  scale_color_manual(values = c("TIRTL-seq" = "grey20", "CAIRR-seq" = "grey20")) +
  labs(title = "Unique TCR Sequence Counts: TIRTL-seq vs CAIRR-seq",
       x = NULL, y = "Unique Sequence Count", fill = "Technology",
       caption = "Bars: geometric mean \u00b1 SD (log scale). Brackets: Wilcoxon rank-sum test\n(* p<0.05, ns = not significant)"
  ) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.minor.ticks.length = unit(0.075, "cm"),
        axis.ticks = element_line(color = "grey30"))

print(uniq_p)

## ---- 2.6 Save plot -----------------------------------------------------------
ggsave(file.path(uniq_out_dir, "barplot_TIRTL_vs_CAIRR.pdf"), uniq_p, width = 6, height = 4)
ggsave(file.path(uniq_out_dir, "barplot_TIRTL_vs_CAIRR.png"), uniq_p, width = 6, height = 4, dpi = 600)

cat("PART 2 done: saved barplot_TIRTL_vs_CAIRR.pdf/.png and source data CSVs\n")

## =============================================================================
## Summary of outputs:
##   PART 1 (CAIRR-seq vs TIRTL-seq, total read depth, MP036/data/outputs/):
##     read_depth_omni_vs_tirtl.csv           - raw per-sample totals (data/, not outputs/)
##     plots/omni_vs_tirtl_read_depth_summary.png / .pdf
##     summary_omni_vs_tirtl.csv
##   PART 2 (CAIRR-seq vs TIRTL-seq, MP036/plots/):
##     source_data_per_sample_long.csv / _wide.csv
##     source_data_summary_stats.csv
##     source_data_wilcox.csv
##     barplot_TIRTL_vs_CAIRR.pdf / .png
## =============================================================================
