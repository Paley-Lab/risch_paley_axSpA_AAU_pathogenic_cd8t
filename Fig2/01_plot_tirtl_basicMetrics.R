## =============================================================================
## TIRTL-seq TCR Metrics (TIRTL-seq only)
##
## PART 1: TIRTL-seq -- total read/UMI depth, TRA vs TRB
## PART 2: TIRTL-seq -- unique sequence counts, Alpha/Beta/Paired
## =============================================================================

## ---- 0. Setup ------------------------------------------------------------------
library(tidyverse)
library(scales)

TIRTL_COLOR <- "#D9782D"   # kept consistent with TIRTL-seq's color in the prior comparison plots

## #############################################################################
## PART 1: TIRTL-seq -- total read depth (TRA vs TRB)
## #############################################################################

## ---- 1.1 Load & aggregate TIRTL data ---------------------------------------------
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

## ---- 1.2 Combine & export raw depth data -----------------------------------------
depth_bound <- bind_rows(
  tirtl_total_alpha_per_sample[, c("sample", "locus", "totalCounts", "tech")],
  tirtl_total_beta_per_sample[, c("sample", "locus", "totalCounts", "tech")]
)

depth_out_path <- "/storage1/fs1/paleym/Active/Isabel_Risch/MP036/data/umi_counts_tirtl.csv"
write.csv(depth_bound, file = depth_out_path, row.names = FALSE)

## ---- 1.3 Tidy for plotting --------------------------------------------------------
depth_d <- depth_bound %>%
  mutate(locus = factor(locus, levels = c("TRA", "TRB")))

## ---- 1.4 Summary statistics (log10 space) -----------------------------------------
## Geometric mean +/- SD, i.e. mean/SD of log10(totalCounts) then exponentiated
## for display -- appropriate since read counts span multiple orders of magnitude.
depth_summary_df <- depth_d %>%
  group_by(locus) %>%
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

## ---- 1.5 Axis range (fixed 10^0 to 10^9, per request) -----------------------------
depth_exp_range <- 0:9
depth_lower_limit <- 10^0
depth_upper_limit <- 10^9

## ---- 1.6 Plot -----------------------------------------------------------------------
## Bars are drawn with geom_rect() and an explicit ymin = 10^0, NOT geom_col()
## (which has an implicit ymin = 0). Under scale_y_log10(), log10(0) = -Inf,
## which would fall outside an explicit lower "limits" bound and get the
## whole bar censored/dropped -- the same "bars disappear" issue from
## earlier. Giving the bar an explicit, non-zero floor avoids that entirely.
depth_d <- depth_d %>% mutate(x_num = as.numeric(locus))
depth_summary_df <- depth_summary_df %>% mutate(x_num = as.numeric(locus))

depth_p <- ggplot(depth_d, aes(x = x_num, y = totalCounts)) +
  geom_rect(data = depth_summary_df,
            aes(xmin = x_num - 0.3, xmax = x_num + 0.3, ymin = depth_lower_limit, ymax = geo_mean),
            inherit.aes = FALSE, fill = TIRTL_COLOR, colour = "grey20") +
  geom_errorbar(data = depth_summary_df,
                aes(x = x_num, y = geo_mean, ymin = geo_sd_low, ymax = geo_sd_high),
                inherit.aes = FALSE, width = 0.15) +
  geom_jitter(width = 0.1, shape = 21, fill = "white", colour = "grey20", size = 2) +
  scale_x_continuous(breaks = seq_along(levels(depth_d$locus)), labels = levels(depth_d$locus),
                     limits = c(0.5, length(levels(depth_d$locus)) + 0.5), expand = c(0, 0)) +
  scale_y_log10(breaks = 10^depth_exp_range,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks = NULL,
                limits = c(depth_lower_limit, depth_upper_limit)) +
  labs(title = "TIRTL-seq TCR Read Depth: TRA vs TRB",
       x = NULL, y = "Unique Reads Per Sample",
       caption = "Bars: geometric mean \u00b1 SD (log scale)."
  ) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.minor.ticks.length = unit(0.075, "cm"),
        axis.ticks = element_line(color = "grey30"))

print(depth_p)

## ---- 1.7 Save outputs ----------------------------------------------------------------
setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP036")
dir.create("data/outputs/plots", showWarnings = FALSE, recursive = TRUE)

ggsave("data/outputs/plots/tirtl_read_depth_summary.png", depth_p, width = 3, height = 4, dpi = 600)
ggsave("data/outputs/plots/tirtl_read_depth_summary.pdf", depth_p, width = 3, height = 4)

write.csv(depth_summary_df, file = "data/outputs/summary_tirtl_read_depth.csv", row.names = FALSE)

cat("PART 1 done: saved tirtl_read_depth_summary.png/.pdf and summary_tirtl_read_depth.csv\n")


## #############################################################################
## PART 2: TIRTL-seq -- unique sequence count comparison (Alpha/Beta/Paired)
## #############################################################################

uniq_out_dir <- "/storage1/fs1/paleym/Active/Isabel_Risch/MP036/plots"
dir.create(uniq_out_dir, showWarnings = FALSE)

## ---- 2.1 Load & tidy data ---------------------------------------------------------
uniq_tirtl <- read.csv("/storage1/fs1/paleym/Active/Isabel_Risch/MP028/data/TIRTLseq_unique_sequence_counts_per_sample.csv") %>%
  mutate(tech = "TIRTL-seq")

uniq_d <- uniq_tirtl %>%
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

write_csv(uniq_d,     file.path(uniq_out_dir, "source_data_per_sample_long_tirtl.csv"))
write_csv(uniq_tirtl, file.path(uniq_out_dir, "source_data_per_sample_wide_tirtl.csv"))

## ---- 2.2 Summary statistics (log10 space) -----------------------------------------
uniq_summary_df <- uniq_d %>%
  group_by(locus) %>%
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

write_csv(uniq_summary_df, file.path(uniq_out_dir, "unique_seq_summary_stats_tirtl.csv"))

print(uniq_summary_df)

## ---- 2.3 Axis range (computed from the data, for a log10 y-axis) -----------------
uniq_exp_low    <- 0
uniq_exp_high   <- ceiling(log10(max(uniq_d$totalCounts, uniq_summary_df$geo_sd_high)) + 0.2)
uniq_exp_range  <- uniq_exp_low:uniq_exp_high
uniq_upper_limit <- 10^uniq_exp_high

## ---- 2.4 Plot --------------------------------------------------------------------
uniq_p <- ggplot(uniq_d, aes(x = locus, y = totalCounts)) +
  geom_col(data = uniq_summary_df, aes(y = geo_mean), width = 0.6,
           fill = TIRTL_COLOR, colour = "grey20") +
  geom_errorbar(data = uniq_summary_df,
                aes(y = geo_mean, ymin = geo_sd_low, ymax = geo_sd_high),
                width = 0.15) +
  geom_jitter(width = 0.1, shape = 21, fill = "white", colour = "grey20", size = 2) +
  scale_y_log10(breaks = 10^uniq_exp_range,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks = NULL,
                limits = c(NA, uniq_upper_limit)) +
  labs(title = "TIRTL-seq Unique TCR Sequence Counts",
       x = NULL, y = "Number of Unique Clones",
       caption = "Bars: geometric mean \u00b1 SD (log scale)."
  ) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.minor.ticks.length = unit(0.075, "cm"),
        axis.ticks = element_line(color = "grey30"))

print(uniq_p)

## ---- 2.5 Save plot -----------------------------------------------------------
ggsave(file.path(uniq_out_dir, "barplot_TIRTL_only.pdf"), uniq_p, width = 3.5, height = 4)
ggsave(file.path(uniq_out_dir, "barplot_TIRTL_only.png"), uniq_p, width = 3.5, height = 4, dpi = 600)

cat("PART 2 done: saved barplot_TIRTL_only.pdf/.png and source data CSVs\n")

## =============================================================================
## Summary of outputs:
##   PART 1 (TIRTL-seq, total read depth, MP036/data/outputs/):
##     umi_counts_tirtl.csv                  - raw per-sample totals (data/, not outputs/)
##     plots/tirtl_read_depth_summary.png / .pdf
##     summary_tirtl_read_depth.csv
##   PART 2 (TIRTL-seq, unique sequence counts, MP036/plots/):
##     source_data_per_sample_long_tirtl.csv / _wide_tirtl.csv
##     unique_seq_summary_stats_tirtl.csv
##     barplot_TIRTL_only.pdf / .png
## =============================================================================