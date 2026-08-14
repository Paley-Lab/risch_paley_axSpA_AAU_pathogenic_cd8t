setwd("/Volumes/Active/Isabel_Risch/MP002/")
library(readxl)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)   # Prism-style beeswarm dot plots
library(ggsignif)     # significance brackets
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

df <- read_excel(input_file, sheet = 1)
df$Sample <- toupper(df$Sample)

# add in some metadata stuff
tmp <- read.csv("data/outputs/plots/tidy_combined_data.csv")

df <- merge(df, unique(tmp[,c('sample', "YeiHStatus")]), by=1)


write.csv(dat, file.path(out_dir, "summary_data_trbv_genes.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 2. Prism-style dot plots
##    Rows = Test, Columns = Prior probability, x = Healthy/Disease,
##    y = % pathogenic cells. Mean +/- SD bars overlaid on jittered points,
##    styled to resemble a GraphPad Prism scatter/dot plot.
## ---------------------------------------------------------------------------

# summary stats (mean +/- SD) per Test x YeiHStatus x group, used to draw bars
summary_stats <- df %>%
  group_by(YeiHStatus, Diagnosis) %>%
  summarise(
    mean = mean(`TRBV9/5-5/5-4 Percentage`),
    sd   = sd(`TRBV9/5-5/5-4 Percentage`),
    n    = n(),
    .groups = "drop"
  )

## Non-parametric test (Wilcoxon rank-sum / Mann-Whitney U) comparing
## Healthy vs Disease within each Test x YeiHStatus combination.
wilcox_results <- df %>%
  group_by(YeiHStatus) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(`TRBV9/5-5/5-4 Percentage` ~ Diagnosis)$p.value,
      error = function(e) NA_real_
    ),
    y_max = max(`TRBV9/5-5/5-4 Percentage`),
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
          file.path(out_dir, "wilcoxon_test_results_trbv.csv"), row.names = FALSE)
print(wilcox_results)

# only draw brackets for statistically significant comparisons (p < 0.05)
# Built with plain geom_segment/geom_text (rather than ggsignif) so bracket
# placement inside faceted, patchwork-combined plots is robust across
# ggplot2/ggsignif versions.
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


## Build one row of panels (facet_wrap, scales = "free_y") per prior-probability
## level so every single panel gets its own, fully independent y-axis - unlike
## facet_grid, where "free_y" only frees the scale per row of the grid and
## panels within a row still share one axis (which squashes small values).
make_dot_row <- function(prior_level, show_x_labels) {
  
  sub_dat      <- df           %>% filter(YeiHStatus == prior_level)
  sub_summary  <- summary_stats %>% filter(YeiHStatus == prior_level)
  sub_sig      <- sig_brackets  %>% filter(YeiHStatus == prior_level)
  
  p <- ggplot(sub_dat, aes(x = Diagnosis, y = `TRBV9/5-5/5-4 Percentage`,
                           color = Diagnosis)) +
    geom_errorbar(
      data = sub_summary,
      aes(x = Diagnosis, ymin = mean - sd, ymax = mean + sd, y = mean),
      inherit.aes = FALSE, width = 0.25, linewidth = 0.7, color = "black"
    ) +
    geom_errorbar(
      data = sub_summary,
      aes(x = Diagnosis, ymin = mean, ymax = mean, y = mean),
      inherit.aes = FALSE, width = 0.4, linewidth = 0.9, color = "black"
    ) +
    geom_beeswarm(cex = 2.2, size = 2.4, alpha = 0.9) +
    { if (nrow(sub_sig) > 0)
      list(
        # horizontal bar of the bracket
        geom_segment(
          data = sub_sig,
          aes(x = x, xend = xend, y = y_position, yend = y_position),
          inherit.aes = FALSE, linewidth = 0.6, color = "black"
        ),
        # left tick
        geom_segment(
          data = sub_sig,
          aes(x = x, xend = x, y = y_position - tick, yend = y_position),
          inherit.aes = FALSE, linewidth = 0.6, color = "black"
        ),
        # right tick
        geom_segment(
          data = sub_sig,
          aes(x = xend, xend = xend, y = y_position - tick, yend = y_position),
          inherit.aes = FALSE, linewidth = 0.6, color = "black"
        ),
        # asterisk label
        geom_text(
          data = sub_sig,
          aes(x = (x + xend) / 2, y = y_position, label = annotation),
          inherit.aes = FALSE, vjust = -0.3, size = 5
        )
      )
    } +
    scale_color_manual(values = c('AI' = "#EC3B93", 'HC' = "#499B78")) +
    expand_limits(y = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = paste0(prior_level, "\nUMIs with TRBV9/5-5/5-4 (%)")) +
    prism_theme +
    theme(
      strip.text = if (show_x_labels) element_text(face = "bold") else element_blank(),
      strip.background = if (show_x_labels) element_rect(fill = "grey90", color = "black") else element_blank(),
      axis.text.x = if (show_x_labels) element_text() else element_blank(),
      axis.ticks.x = if (show_x_labels) element_line() else element_blank()
    )
  
  p
}

row_low  <- make_dot_row("YeiH Low",  show_x_labels = TRUE)
row_high <- make_dot_row("YeiH High", show_x_labels = TRUE)

dot_plot <- row_low / row_high +
  plot_annotation(
    title = NULL,
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
  )

ggsave(file.path(out_dir, "4_prism_dot_plots_TRBV.pdf"), dot_plot, width = 2.5, height = 6)
ggsave(file.path(out_dir, "4_prism_dot_plots_TRBV.png"), dot_plot, width = 2.5, height = 6, dpi = 300)
