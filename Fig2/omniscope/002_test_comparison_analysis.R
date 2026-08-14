## ============================================================================
## Comparing three tests (cdr3b, cdr3b_trbv, cdr3b_trbv_trav) for their
## ability to distinguish Disease vs Healthy samples, stratified by the
## prior probability (YeiHStatus: low / high) of the samples tested.
##
## Outputs:
##   1. Prism-style dot plots (Healthy vs Disease), faceted by test and by
##      prior probability
##   2. Two ROC plots (one for "low" prior-probability samples, one for
##      "high" prior-probability samples), each showing all three tests
##   3. Sensitivity / specificity tables for each test, split by prior
##      probability (evaluated at the Youden-optimal threshold)
##
## Input:  claude_input.xlsx  (3 sheets, one per test: cdr3b, cdr3b_trbv,
##         cdr3b_trbv_trav; columns: sample, proportion, Patient.description,
##         YeiHStatus)
## ============================================================================

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

input_file  <- "data/outputs/omniscope_pathogenic_proportions_allAnalyses.xlsx"
out_dir     <- "data/outputs/plots"
dir.create(out_dir, showWarnings = FALSE)

test_sheets <- c("cdr3b", "cdr3b_trbv", "cdr3b_trbv_trav")
test_labels <- c(
  cdr3b            = "CDR3b",
  cdr3b_trbv       = "CDR3b + TRBV",
  cdr3b_trbv_trav  = "CDR3b + TRBV + TRAV"
)

## ---------------------------------------------------------------------------
## 1. Load & tidy the data
## ---------------------------------------------------------------------------

raw_list <- lapply(test_sheets, function(s) {
  df <- read_excel(input_file, sheet = s)
  df$Test <- s
  df
})

dat <- bind_rows(raw_list) %>%
  mutate(
    Test               = factor(Test, levels = test_sheets, labels = test_labels[test_sheets]),
    Patient.description = factor(Patient.description, levels = c("Healthy", "AAU", "axSpa, AAU", "axSpA, AAU", "axSpa"),
                                 labels = c("Healthy", "Disease", "Disease", "Disease", "Disease")),
    YeiHStatus          = factor(Sample.condition, levels = c("YEIH low control", "YEIH high control",
                                                              "YEIH high", "YEIH low"),
                                labels = c("YeiH Low", "YeiH High", "YeiH High", "YeiH Low")),
    pct_pathogenic     = proportion   # already expressed as a percentage - no conversion needed
  )

write.csv(dat, file.path(out_dir, "tidy_combined_data.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 2. Prism-style dot plots
##    Rows = Test, Columns = Prior probability, x = Healthy/Disease,
##    y = % pathogenic cells. Mean +/- SD bars overlaid on jittered points,
##    styled to resemble a GraphPad Prism scatter/dot plot.
## ---------------------------------------------------------------------------

# summary stats (mean +/- SD) per Test x YeiHStatus x group, used to draw bars
summary_stats <- dat %>%
  group_by(Test, YeiHStatus, Patient.description) %>%
  summarise(
    mean = mean(pct_pathogenic),
    sd   = sd(pct_pathogenic),
    n    = n(),
    .groups = "drop"
  )

## Non-parametric test (Wilcoxon rank-sum / Mann-Whitney U) comparing
## Healthy vs Disease within each Test x YeiHStatus combination.
wilcox_results <- dat %>%
  group_by(Test, YeiHStatus) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(pct_pathogenic ~ Patient.description)$p.value,
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
          file.path(out_dir, "wilcoxon_test_results.csv"), row.names = FALSE)
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

## Three-color palette for the three tests: distinct against a white
## background and clearly distinguishable from one another (based on the
## Okabe-Ito colorblind-safe palette; blue darkened slightly for extra
## contrast against the bluish green).
test_palette <- c(
  "CDR3b"                = "#009E73",  # bluish green
  "CDR3b + TRBV"          = "#004C73",  # blue (darkened from Okabe-Ito #0072B2)
  "CDR3b + TRBV + TRAV"   = "#CC79A7"   # Okabe-Ito pink (reddish purple)
)

## Build one row of panels (facet_wrap, scales = "free_y") per prior-probability
## level so every single panel gets its own, fully independent y-axis - unlike
## facet_grid, where "free_y" only frees the scale per row of the grid and
## panels within a row still share one axis (which squashes small values).
make_dot_row <- function(prior_level, show_x_labels) {
  
  sub_dat      <- dat           %>% filter(YeiHStatus == prior_level)
  sub_summary  <- summary_stats %>% filter(YeiHStatus == prior_level)
  sub_sig      <- sig_brackets  %>% filter(YeiHStatus == prior_level)
  
  p <- ggplot(sub_dat, aes(x = Patient.description, y = pct_pathogenic,
                           color = Test)) +
    geom_errorbar(
      data = sub_summary,
      aes(x = Patient.description, ymin = mean - sd, ymax = mean + sd, y = mean),
      inherit.aes = FALSE, width = 0.25, linewidth = 0.7, color = "black"
    ) +
    geom_errorbar(
      data = sub_summary,
      aes(x = Patient.description, ymin = mean, ymax = mean, y = mean),
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
    facet_wrap(~ Test, nrow = 1, scales = "free_y") +
    scale_color_manual(values = test_palette) +
    expand_limits(y = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = paste0(prior_level, "\nPathogenic cells (%)")) +
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
    title = "Pathogenic cell frequency by test, group, and prior probability",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
  )

ggsave(file.path(out_dir, "1_prism_dot_plots.pdf"), dot_plot, width = 7, height = 6)
ggsave(file.path(out_dir, "1_prism_dot_plots.png"), dot_plot, width = 7, height = 6, dpi = 300)

## ---------------------------------------------------------------------------
## 3. ROC curves
##    One plot per prior-probability stratum, each containing all 3 tests.
##    Positive class = "Disease".
## ---------------------------------------------------------------------------

make_roc_plot <- function(prior_level, title_txt) {
  
  sub <- dat %>% filter(YeiHStatus == prior_level)
  
  roc_list <- lapply(levels(sub$Test), function(t) {
    d <- sub %>% filter(Test == t)
    roc(d$Patient.description, d$proportion,
        levels = c("Healthy", "Disease"), direction = "<")
  })
  names(roc_list) <- levels(sub$Test)
  
  # build a tidy data frame of the ROC coordinates for ggplot
  roc_df <- bind_rows(lapply(names(roc_list), function(t) {
    r <- roc_list[[t]]
    data.frame(
      Test        = t,
      Sensitivity = rev(r$sensitivities),
      Specificity = rev(r$specificities)
    )
  }))
  roc_df$Test <- factor(roc_df$Test, levels = levels(sub$Test))
  
  aucs <- sapply(roc_list, function(r) as.numeric(auc(r)))
  legend_labels <- paste0(names(aucs), " (AUC = ", sprintf("%.2f", aucs), ")")
  names(legend_labels) <- names(aucs)
  
  p <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Test)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = test_palette, labels = legend_labels) +
    coord_equal() +
    labs(
      title = title_txt,
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
  
  list(plot = p, roc_list = roc_list, aucs = aucs)
}

roc_low  <- make_roc_plot("YeiH Low",  "ROC: YeiH Low samples")
roc_high <- make_roc_plot("YeiH High", "ROC: YeiH High samples")

ggsave(file.path(out_dir, "2_ROC_YeiH_Low.pdf"),  roc_low$plot,  width = 6, height = 7)
ggsave(file.path(out_dir, "2_ROC_YeiH_Low.png"),  roc_low$plot,  width = 6, height = 7, dpi = 300)
ggsave(file.path(out_dir, "2_ROC_YeiH_High.pdf"), roc_high$plot, width = 6, height = 7)
ggsave(file.path(out_dir, "2_ROC_YeiH_High.png"), roc_high$plot, width = 6, height = 7, dpi = 300)

# combined side-by-side figure for convenience
combined_roc <- gridExtra::arrangeGrob(roc_low$plot, roc_high$plot, ncol = 2)
ggsave(file.path(out_dir, "2_ROC_combined.pdf"), combined_roc, width = 6, height = 4)
ggsave(file.path(out_dir, "2_ROC_combined.png"), combined_roc, width = 6, height = 4, dpi = 300)

## ---------------------------------------------------------------------------
## 4. Sensitivity / specificity tables (Youden-optimal threshold), by test
##    and prior-probability stratum
## ---------------------------------------------------------------------------

build_sens_spec_table <- function(prior_level) {
  sub <- dat %>% filter(YeiHStatus == prior_level)
  
  rows <- lapply(levels(sub$Test), function(t) {
    d <- sub %>% filter(Test == t)
    r <- roc(d$Patient.description, d$proportion,
             levels = c("Healthy", "Disease"), direction = "<")
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
      YeiHStatus          = prior_level,
      Test                = t,
      N_Healthy           = sum(d$Patient.description == "Healthy"),
      N_Disease           = sum(d$Patient.description == "Disease"),
      AUC                 = as.numeric(auc(r)),
      Optimal_threshold_pct = best$threshold,
      Sensitivity         = best$sensitivity,
      Specificity         = best$specificity,
      PPV                 = best$ppv,
      NPV                 = best$npv,
      Accuracy            = best$accuracy,
      Sensitivity_ThreshOfZero         = zerothresh$sensitivity,
      Specificity_ThreshOfZero         = zerothresh$specificity,
      Accuracy_ThreshOfZero         = zerothresh$accuracy
    )
  })
  
  bind_rows(rows)
}

sens_spec_low  <- build_sens_spec_table("YeiH Low")
sens_spec_high <- build_sens_spec_table("YeiH High")
sens_spec_all  <- bind_rows(sens_spec_low, sens_spec_high) %>%
  mutate(across(c(AUC, Sensitivity, Specificity, PPV, NPV, Accuracy), ~ signif(.x, 3)),
         Optimal_threshold_pct = signif(Optimal_threshold_pct, 3))

print(sens_spec_all)

write.csv(sens_spec_all, file.path(out_dir, "3_sensitivity_specificity_table.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "YeiH_Low")
writeData(wb, "YeiH_Low",
          sens_spec_all %>% filter(YeiHStatus == "YeiH Low"))
addWorksheet(wb, "YeiH_High")
writeData(wb, "YeiH_High",
          sens_spec_all %>% filter(YeiHStatus == "YeiH High"))
saveWorkbook(wb, file.path(out_dir, "3_sensitivity_specificity_table.xlsx"), overwrite = TRUE)

## ---------------------------------------------------------------------------
## 5. Session info (for reproducibility)
## ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\nAll outputs written to:", normalizePath(out_dir), "\n")
