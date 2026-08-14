# Check how many reads are in the Omniscope and TIRTL data and compare sequencing depths

##################################################
############### TOTAL UMI COUNT ##################
##################################################

# ------- Omniscope --------- #

datadir <- "/storage1/fs1/paleym/Active/SeqData/Omniscope/"
contig.files <- list.files(datadir, paste0(".csv$", collapse = "|"), recursive = TRUE, full.names = TRUE)
omni <- lapply(contig.files, read.csv)

# get the number of reads total per sample (split by alpha/beta)

omni_total_per_sample <- list()
for (i in 1:length(omni)) {
  omni_total_per_sample[[i]] <- omni[[i]] %>% 
    group_by(locus) %>%
    summarise(totalCounts = sum(umi_count))
}

omni_bound <- bind_rows(omni_total_per_sample)
omni_bound$tech <- "omniscope"

# ------- TIRTL --------- #

setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP028")

all_samples_alpha <- read.csv("data/outputs/all_alpha_data.csv")
all_samples_beta <- read.csv("data/outputs/all_alpha_data.csv")

# get the number of reads total per sample (split by alpha/beta)
tirtl_total_alpha_per_sample <- all_samples_alpha %>% 
  group_by(sample) %>%
  summarise(totalCounts = sum(readCount))
tirtl_total_beta_per_sample <- all_samples_beta %>% 
  group_by(sample) %>%
  summarise(totalCounts = sum(readCount))

tirtl_total_alpha_per_sample$tech <- "TIRTL"
tirtl_total_alpha_per_sample$locus <- "TRA"
tirtl_total_beta_per_sample$tech <- "TIRTL"
tirtl_total_beta_per_sample$locus <- "TRB"

tirtl_total_alpha_per_sample <- as.data.frame(tirtl_total_alpha_per_sample)
tirtl_total_beta_per_sample <- as.data.frame(tirtl_total_beta_per_sample)
total_bound <- bind_rows(omni_bound, tirtl_total_alpha_per_sample[,c("locus", "totalCounts", "tech")]) %>% 
  bind_rows(tirtl_total_beta_per_sample[,c("locus", "totalCounts", "tech")])

write.csv(total_bound, file = "/storage1/fs1/paleym/Active/Isabel_Risch/MP036/data/read_depth_omni_vs_tirtl.csv", row.names = F)


## =============================================================================
## Compare TIRTL-seq vs CAIRR-seq unique sequence counts
## PRISM (GraphPad)-style barplots with individual points, mean +/- SEM,
## and Wilcoxon rank-sum significance testing
## =============================================================================

library(tidyverse)
library(ggpubr)

setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP036")

out_dir <- "plots"
dir.create(out_dir, showWarnings = FALSE)

## ---- 1. Load & tidy data -----------------------------------------------------
cairr <- read.csv("/storage1/fs1/paleym/Active/Isabel_Risch/MP002/data/CAIRRseq_unique_sequence_counts_per_sample.csv") %>%
  mutate(technology = "CAIRR-seq")
tirtl <- read.csv("/storage1/fs1/paleym/Active/Isabel_Risch/MP028/data/TIRTLseq_unique_sequence_counts_per_sample.csv") %>%
  mutate(technology = "TIRTL-seq")

combined_wide <- bind_rows(tirtl, cairr) %>%
  mutate(tech = factor(technology, levels = c("CAIRR-seq","TIRTL-seq")))

# Long format -- one row per sample x metric ("locus"-equivalent grouping)
d <- combined_wide %>%
  pivot_longer(
    cols = c(numberUniqueAlphas, numberUniqueBetas, numberUniquePaired),
    names_to = "locus",
    values_to = "totalCounts"
  ) %>%
  mutate(
    locus = recode(locus,
                   numberUniqueAlphas = "Alpha",
                   numberUniqueBetas  = "Beta",
                   numberUniquePaired = "Paired"),
    locus = factor(locus, levels = c("Alpha", "Beta", "Paired"))
  )

## Source data #1: raw per-sample values used in the plot
write_csv(d, file.path(out_dir, "source_data_per_sample_long.csv"))
write_csv(combined_wide, file.path(out_dir, "source_data_per_sample_wide.csv"))

## ---- 2. Summary statistics (log10 space) --------------------------------------
## Geometric mean +/- SD, i.e. mean/SD of log10(count) then exponentiated for
## display -- appropriate since these counts span multiple orders of magnitude.
summary_df <- d %>%
  group_by(locus, tech) %>%
  summarise(
    n           = n(),
    mean    = mean(totalCounts),
    # sd          = sd(totalCounts),
    mean_log    = mean(log10(totalCounts)),
    sd_log      = sd(log10(totalCounts)),
    geo_mean    = 10^mean_log,
    geo_sd_low  = 10^(mean_log - sd_log),
    geo_sd_high = 10^(mean_log + sd_log),
    .groups = "drop"
  )

## Source data #2: summary stats table (geometric mean, log-space mean/SD, n)
write_csv(summary_df, file.path(out_dir, "source_data_summary_stats.csv"))

## ---- 3. Statistical testing (technology comparison, per locus) ----------------
## Small, unequal, likely non-normal sample sizes -> Wilcoxon rank-sum (Mann-Whitney)
stat_tests <- d %>%
  group_by(locus) %>%
  summarise(
    p_value = wilcox.test(totalCounts ~ tech)$p.value,
    .groups = "drop"
  ) %>%
  mutate(stars = ifelse(p_value < 0.05, "*", "ns"))

write_csv(stat_tests, file.path(out_dir, "source_data_wilcox.csv"))

print(summary_df)
print(stat_tests)

## ---- 4. Bracket geometry ------------------------------------------------------
## Each locus gets its own bracket, positioned just above that locus's own
## local max (bar top or point), since Alpha/Beta and Paired counts differ by
## orders of magnitude and a single shared height would misplace most brackets.
dodge_width  <- 0.75
bar_offset   <- 0.19          # approx. half-distance between the two dodged bars

locus_levels <- levels(d$locus)
n_locus      <- length(locus_levels)

bracket_geom <- map_dfr(seq_along(locus_levels), function(i) {
  lv <- locus_levels[i]
  local_max <- max(
    d$totalCounts[d$locus == lv],
    summary_df$geo_sd_high[summary_df$locus == lv]
  )
  log_max   <- log10(local_max)
  bracket_y <- log_max + 0.12
  label_y   <- log_max + 0.22
  tibble(
    locus   = lv,
    x       = i - bar_offset,
    xend    = i + bar_offset,
    y       = bracket_y,
    yend    = bracket_y,
    label_x = i,
    label_y = label_y
  )
}) %>%
  left_join(stat_tests, by = "locus")

## Overall axis range: span from the smallest to the largest value across all
## loci/technologies, plus headroom for the highest bracket label
exp_low  <- 0
exp_high <- ceiling(max(bracket_geom$label_y)) 
exp_range <- exp_low:exp_high
upper_limit <- 10^(exp_high)

## ---- 5. Plot --------------------------------------------------------------------

p <- ggplot(d, aes(x = locus, y = totalCounts, fill = tech)) +
  geom_col(data = summary_df, aes(y = geo_mean), position = position_dodge(dodge_width),
           width = 0.7, color = "grey20") +
  geom_errorbar(data = summary_df,
                aes(y = geo_mean, ymin = geo_sd_low, ymax = geo_sd_high),
                position = position_dodge(dodge_width), width = 0.15) +
  geom_point(aes(color = tech),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width),
             shape = 21, fill = "white", size = 2, show.legend = FALSE) +
  # ggplot2::annotate("segment", x = bracket_geom$x, xend = bracket_geom$xend,
  #                   y = 10^bracket_geom$y, yend = 10^bracket_geom$yend) +
  # ggplot2::annotate("text", x = bracket_geom$label_x, y = 10^bracket_geom$label_y,
  #                   label = bracket_geom$stars, size = 4.5, fontface = "bold") +
  scale_y_log10(breaks = 10^exp_range,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks = NULL,
                limits = c(NA, upper_limit)) +
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

print(p)

## ---- 6. Save plot -----------------------------------------------------------
ggsave(file.path(out_dir, "barplot_TIRTL_vs_CAIRR.pdf"), p, width = 6, height = 4)
ggsave(file.path(out_dir, "barplot_TIRTL_vs_CAIRR.png"), p, width = 6, height = 4, dpi = 600)

## =============================================================================
## Outputs written to ./prism_output/:
##   source_data_per_sample_long.csv   - tidy per-sample values (as plotted)
##   source_data_per_sample_wide.csv   - original per-sample values, tech-labeled
##   source_data_summary_stats.csv     - geometric mean, log-space mean/SD, n
##   source_data_statistics.csv        - Wilcoxon p-values and significance labels
##   barplot_TIRTL_vs_CAIRR.pdf / .png - combined dodged barplot, all 3 loci
## =============================================================================




