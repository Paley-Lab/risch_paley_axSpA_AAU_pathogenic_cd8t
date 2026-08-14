setwd("/Volumes/Active/Isabel_Risch/MP028/")

library(readxl)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)   # Prism-style beeswarm dot plots
library(pROC)
library(gridExtra)
library(openxlsx)
library(patchwork)

## ---------------------------------------------------------------------------
## 0. Setup
## ---------------------------------------------------------------------------

input_file  <- "data/outputs/betaOnly_table_TRBV_percentage_by_subject.xlsx"
out_dir     <- "data/outputs/plots"
dir.create(out_dir, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. Load & tidy the data
## ---------------------------------------------------------------------------

raw_list <- read_excel(input_file, sheet = 1)

dat <- bind_rows(raw_list) %>%
  mutate(
    diagnosis          = factor(Diagnosis, levels = c("HC", "AI")),
    pct_pathogenic     = `TRBV9/5-5/5-4 Proportion`   # already expressed as a percentage
  )

## ---------------------------------------------------------------------------
## 2. Prism-style dot plots
##    One column per test, dots colored by test, x = HC/AI, y = %
##    pathogenic cells. Mean +/- SD bars overlaid on jittered points, each
##    panel with its own independent y-axis (small values don't get
##    squashed by larger ones).
## ---------------------------------------------------------------------------

summary_stats <- dat %>%
  group_by(diagnosis) %>%
  summarise(
    mean = mean(pct_pathogenic),
    sd   = sd(pct_pathogenic),
    n    = n(),
    .groups = "drop"
  )

## Non-parametric test (Wilcoxon rank-sum / Mann-Whitney U) comparing
## HC vs AI within each test.
wilcox_results <- dat %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(pct_pathogenic ~ diagnosis)$p.value,
      error = function(e) NA_real_
    ),
    y_max = max(pct_pathogenic),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value)      ~ "",
      p_value < 0.001      ~ "***",
      p_value < 0.01       ~ "**",
      p_value < 0.05       ~ "*",
      TRUE                ~ "ns"
    )
  )

write.csv(wilcox_results %>% mutate(p_value = signif(p_value, 3)),
          file.path(out_dir, "wilcoxon_test_results_TRBV.csv"), row.names = FALSE)
print(wilcox_results)

# only draw brackets for statistically significant comparisons (p < 0.05)
# Built with plain geom_segment/geom_text (rather than ggsignif) for
# robustness across ggplot2/ggsignif versions.
sig_brackets <- wilcox_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  mutate(
    x          = 1,
    xend       = 2,
    y_position = y_max * 1.12,
    tick       = y_max * 0.05,
    annotation = signif_label
  )

prism_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.8),
    strip.background   = element_rect(fill = "grey90", color = "black"),
    strip.text         = element_text(face = "bold"),
    axis.text           = element_text(color = "black"),
    axis.ticks           = element_line(color = "black"),
    legend.position     = "none",
    plot.title          = element_text(face = "bold", hjust = 0.5)
  )

dot_plot <- ggplot(dat, aes(x = diagnosis, y = pct_pathogenic, color = diagnosis)) +
  geom_errorbar(
    data = summary_stats,
    aes(x = diagnosis, ymin = mean - sd, ymax = mean + sd, y = mean),
    inherit.aes = FALSE, width = 0.25, linewidth = 0.7, color = "black"
  ) +
  geom_errorbar(
    data = summary_stats,
    aes(x = diagnosis, ymin = mean, ymax = mean, y = mean),
    inherit.aes = FALSE, width = 0.4, linewidth = 0.9, color = "black"
  ) +
  geom_beeswarm(cex = 2.2, size = 2.4, alpha = 0.9) +
  { if (nrow(sig_brackets) > 0)
    list(
      geom_segment(
        data = sig_brackets,
        aes(x = x, xend = xend, y = y_position, yend = y_position),
        inherit.aes = FALSE, linewidth = 0.6, color = "black"
      ),
      geom_segment(
        data = sig_brackets,
        aes(x = x, xend = x, y = y_position - tick, yend = y_position),
        inherit.aes = FALSE, linewidth = 0.6, color = "black"
      ),
      geom_segment(
        data = sig_brackets,
        aes(x = xend, xend = xend, y = y_position - tick, yend = y_position),
        inherit.aes = FALSE, linewidth = 0.6, color = "black"
      ),
      geom_text(
        data = sig_brackets,
        aes(x = (x + xend) / 2, y = y_position, label = annotation),
        inherit.aes = FALSE, vjust = -0.3, size = 5
      )
    )
  } +
  scale_color_manual(values = c('AI' = "#EC3B93", 'HC' = "#499B78")) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = NULL,
    x = NULL,
    y = "UMIs with TRBV9/5-5/5-4 (%)"
  ) +
  prism_theme

ggsave(file.path(out_dir, "4_prism_dot_plots_TRBV.pdf"), dot_plot, width = 2.5, height = 2.5)
ggsave(file.path(out_dir, "4_prism_dot_plots_TRBV.png"), dot_plot, width = 2.5, height = 2.5, dpi = 600)
