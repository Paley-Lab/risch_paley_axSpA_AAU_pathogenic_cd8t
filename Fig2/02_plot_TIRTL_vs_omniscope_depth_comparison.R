#!/usr/bin/env Rscript
# Summary barplot of TCR-seq read depth (TRA vs TRB) for Omniscope vs TIRTL

d <- read.csv("data/outputs/read_depth_omni_vs_tirtl.csv",
              stringsAsFactors = FALSE)

d$tech  <- factor(d$tech,  levels = c("omniscope", "TIRTL"),
                  labels = c("Omniscope", "TIRTL"))
d$locus <- factor(d$locus, levels = c("TRA", "TRB"))

# ---- summary stats for bars/error bars (computed on log10 scale) ----
summary_df <- aggregate(log10(totalCounts) ~ tech + locus, data = d,
                        FUN = function(x) c(mean = mean(x), sd = sd(x)))
summary_df <- do.call(data.frame, summary_df)
names(summary_df)[3:4] <- c("mean_log", "sd_log")

# ---- Wilcoxon rank-sum test: Omniscope vs TIRTL, within each locus ----
p_values <- sapply(levels(d$locus), function(loc) {
  wilcox.test(totalCounts ~ tech, data = d[d$locus == loc, ], exact = FALSE)$p.value
})
stars <- symnum(p_values, corr = FALSE, na = FALSE,
                cutpoints = c(0, 0.05, 1),
                symbols   = c("*", "ns"))

# ---- plot ----
exp_range <- 0:9  # adjust to match your data's actual range of powers of 10

p <- ggplot(d, aes(x = locus, y = totalCounts, fill = tech)) +
  geom_col(data = summary_df, aes(y = 10^mean_log), position = position_dodge(0.75),
           width = 0.7, color = "grey20") +
  geom_errorbar(data = summary_df,
                aes(y = 10^mean_log, ymin = 10^(mean_log - sd_log), ymax = 10^(mean_log + sd_log)),
                position = position_dodge(0.75), width = 0.15) +
  geom_point(aes(color = tech),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
             shape = 21, fill = "white", size = 2, show.legend = FALSE) +
  ggplot2::annotate("segment", x = c(0.8, 1.8), xend = c(1.2, 2.2),
                    y = c(10^8.35, 10^8.55), yend = c(10^8.35, 10^8.55)) +
  ggplot2::annotate("text", x = c(1, 2), y = c(10^8.42, 10^8.62),
                    label = as.character(stars), size = 4.5, fontface = "bold") +
  scale_y_log10(breaks = 10^exp_range,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks = NULL,
                limits = c(NA, 10^8.75)) +
  scale_fill_manual(values = c("Omniscope" = "#3B7DB8", "TIRTL" = "#D9782D")) +
  scale_color_manual(values = c("Omniscope" = "grey20", "TIRTL" = "grey20")) +
  labs(title = "TCR Read Depth: TRA vs TRB by Technology",
       x = NULL, y = "Total TCR Read Counts", fill = "Technology",
       caption = "Bars: mean \u00b1 SD (log scale). Brackets: Wilcoxon rank-sum test\n(* p<0.05, ns = not significant)"
  ) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.minor.ticks.length = unit(0.075, "cm"),
        axis.ticks = element_line(color = "grey30"))


ggsave("data/outputs/plots/omni_vs_tirtl_read_depth_summary.png", p, width = 5, height = 4, dpi = 600)
ggsave("data/outputs/plots/omni_vs_tirtl_read_depth_summary.pdf", p, width = 5, height = 4)


means <- d %>% group_by(locus, tech) %>% summarise("mean" = mean(totalCounts))

summary_df <- merge(summary_df, means)

write.csv(summary_df, file = "data/outputs/summary_omni_vs_tirtl.csv")

cat("Saved PNG and PDF.\n")
cat("\nWilcoxon rank-sum p-values (Omniscope vs TIRTL):\n")
print(data.frame(locus = levels(d$locus), p_value = p_values, sig = as.character(stars)))