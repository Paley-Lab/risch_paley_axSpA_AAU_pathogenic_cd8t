## ============================================================================
## Comparing three tests (CDR3b, CDR3b + TRBV, CDR3b + TRBV + TRAV) for their
## ability to distinguish AI (disease) vs HC (healthy control) samples.
##
## This is the same analysis as test_comparison_analysis.R, adapted for a
## dataset with no prior-probability stratification (single cohort).
##
## Outputs:
##   1. Prism-style dot plots (HC vs AI), one panel per test, each with its
##      own independent y-axis; significant Wilcoxon comparisons get a
##      bracket
##   2. One ROC plot showing all three tests
##   3. A sensitivity / specificity table for each test (Youden-optimal
##      threshold)
##
## Input:  all_pathogenic_percentage_by_subject.xlsx (3 sheets, one per
##         test: "CDR3b", "CDR3b + TRBV", "CDR3b + TRBV + TRAV"; columns:
##         sample, overallPathProportion, diagnosis [HC / AI])
## ============================================================================
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

input_file  <- "data/outputs/all_pathogenic_percentage_by_subject.xlsx"
out_dir     <- "data/outputs/plots"
dir.create(out_dir, showWarnings = FALSE)

test_sheets <- c("CDR3b", "CDR3b + TRBV", "CDR3b + TRBV + TRAV")

## Three-color palette for the three tests: distinct against a white
## background and clearly distinguishable from one another (based on the
## Okabe-Ito colorblind-safe palette; blue darkened slightly for extra
## contrast against the bluish green).
test_palette <- c(
  "CDR3b"                = "#009E73",  # bluish green
  "CDR3b + TRBV"          = "#004C73",  # blue (darkened from Okabe-Ito #0072B2)
  "CDR3b + TRBV + TRAV"   = "#CC79A7"   # Okabe-Ito pink (reddish purple)
)

## ---------------------------------------------------------------------------
## 1. Load & tidy the data
## ---------------------------------------------------------------------------

raw_list <- lapply(test_sheets, function(s) {
  df <- read_excel(input_file, sheet = s)
  df <- df[, 1:3]
  names(df) <- c("sample", "diagnosis","overallPathProportion")
  df$Test <- s
  df
})

dat <- bind_rows(raw_list) %>%
  mutate(
    Test               = factor(Test, levels = test_sheets),
    diagnosis          = factor(diagnosis, levels = c("HC", "AI")),
    pct_pathogenic     = overallPathProportion   # already expressed as a percentage
  )

write.csv(dat, file.path(out_dir, "tidy_combined_data.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 2. Prism-style dot plots
##    One column per test, dots colored by test, x = HC/AI, y = %
##    pathogenic cells. Mean +/- SD bars overlaid on jittered points, each
##    panel with its own independent y-axis (small values don't get
##    squashed by larger ones).
## ---------------------------------------------------------------------------

summary_stats <- dat %>%
  group_by(Test, diagnosis) %>%
  summarise(
    mean = mean(pct_pathogenic),
    sd   = sd(pct_pathogenic),
    n    = n(),
    .groups = "drop"
  )

## Non-parametric test (Wilcoxon rank-sum / Mann-Whitney U) comparing
## HC vs AI within each test.
wilcox_results <- dat %>%
  group_by(Test) %>%
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
      p_value <= 0.05       ~ "*",
      p_value <= 0.01       ~ "**",
      TRUE                ~ "ns"
    )
  )

write.csv(wilcox_results %>% mutate(p_value = signif(p_value, 3)),
          file.path(out_dir, "wilcoxon_test_results.csv"), row.names = FALSE)
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

dot_plot <- ggplot(dat, aes(x = diagnosis, y = pct_pathogenic, color = Test)) +
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
  facet_wrap(~ Test, nrow = 1, scales = "free_y") +
  scale_color_manual(values = test_palette) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Pathogenic cell frequency by test (HC vs AI)",
    x = NULL,
    y = "Pathogenic cells (%)"
  ) +
  prism_theme

ggsave(file.path(out_dir, "1_prism_dot_plots.pdf"), dot_plot, width = 6.5, height = 3)
ggsave(file.path(out_dir, "1_prism_dot_plots.png"), dot_plot, width = 6.5, height = 3, dpi = 600)

## ---------------------------------------------------------------------------
## 3. ROC curve - all three tests overlaid. Positive class = "AI".
## ---------------------------------------------------------------------------

roc_list <- lapply(levels(dat$Test), function(t) {
  d <- dat %>% filter(Test == t)
  roc(d$diagnosis, d$overallPathProportion,
      levels = c("HC", "AI"), direction = "<")
})
names(roc_list) <- levels(dat$Test)

roc_df <- bind_rows(lapply(names(roc_list), function(t) {
  r <- roc_list[[t]]
  data.frame(
    Test        = t,
    Sensitivity = rev(r$sensitivities),
    Specificity = rev(r$specificities)
  )
}))
roc_df$Test <- factor(roc_df$Test, levels = levels(dat$Test))

aucs <- sapply(roc_list, function(r) as.numeric(auc(r)))
legend_labels <- paste0(names(aucs), " (AUC = ", sprintf("%.2f", aucs), ")")
names(legend_labels) <- names(aucs)

roc_plot <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Test)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = test_palette, labels = legend_labels) +
  coord_equal() +
  labs(
    title = "ROC: AI vs HC",
    x = "1 - Specificity (False positive rate)",
    y = "Sensitivity (True positive rate)",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.background = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(out_dir, "2_ROC.pdf"), roc_plot, width = 4, height = 4)
ggsave(file.path(out_dir, "2_ROC.png"), roc_plot, width = 4, height = 4, dpi = 600)

## ---------------------------------------------------------------------------
## 4. Sensitivity / specificity table (Youden-optimal threshold), by test
## ---------------------------------------------------------------------------

sens_spec_rows <- lapply(levels(dat$Test), function(t) {
  d <- dat %>% filter(Test == t)
  r <- roc_list[[t]]
  best <- coords(r, x = "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity",
                         "ppv", "npv", "accuracy"),
                 transpose = FALSE)
  zerothresh <- coords(r, x = 1e-100, 
                       ret = c("sensitivity", "specificity", "accuracy"),
                       transpose = FALSE)
  # coords() can return >1 row if there are ties; keep the first
  best <- best[1, , drop = FALSE]
  zerothresh <- zerothresh[1, , drop = FALSE]
  
  data.frame(
    Test                   = t,
    N_HC                   = sum(d$diagnosis == "HC"),
    N_AI                   = sum(d$diagnosis == "AI"),
    AUC                    = as.numeric(auc(r)),
    Optimal_threshold_pct  = best$threshold,
    Sensitivity            = best$sensitivity,
    Specificity            = best$specificity,
    PPV                    = best$ppv,
    NPV                    = best$npv,
    Accuracy               = best$accuracy,
    Sensitivity_ThreshOfZero         = zerothresh$sensitivity,
    Specificity_ThreshOfZero         = zerothresh$specificity,
    Accuracy_ThreshOfZero         = zerothresh$accuracy
  )
})

sens_spec_all <- bind_rows(sens_spec_rows) %>%
  mutate(across(c(AUC, Sensitivity, Specificity, PPV, NPV, Accuracy), ~ signif(.x, 3)),
         Optimal_threshold_pct = signif(Optimal_threshold_pct, 3))

print(sens_spec_all)

write.csv(sens_spec_all, file.path(out_dir, "3_sensitivity_specificity_table.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Sensitivity_specificity")
writeData(wb, "Sensitivity_specificity", sens_spec_all)
saveWorkbook(wb, file.path(out_dir, "3_sensitivity_specificity_table.xlsx"), overwrite = TRUE)

## ---------------------------------------------------------------------------
## 5. Session info (for reproducibility)
## ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\nAll outputs written to:", normalizePath(out_dir), "\n")
