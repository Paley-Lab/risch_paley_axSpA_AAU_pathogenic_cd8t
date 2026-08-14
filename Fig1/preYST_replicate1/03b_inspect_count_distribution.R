library(ggplot2)
library(tidyr)
library(dplyr)
library(Biostrings)


# ── 1. Read count matrix ──────────────────────────────────────────────────
count_matrix <- read.csv("outputs/collapsed_count_matrix.csv", row.names = 1)

# ── 2. Define sample groups ───────────────────────────────────────────────────
# Adjust these vectors to match your actual column names
control_cols     <- c('Pre_YST_library_rep1_preSort', 'Pre_YST_library_rep1_no_peptide')
experiment_cols  <- c("Pre_YST_library_rep1_GPER1", 
                       "Pre_YST_library_rep1_PRPF3", 
                      "Pre_YST_library_rep1_PSG5", "Pre_YST_library_rep1_RNASEH2b", 
                      "Pre_YST_library_rep1_yeih")

# ── 3. Reshape to long format ─────────────────────────────────────────────────
long_counts <- count_matrix |>
  pivot_longer(everything(), names_to = "sample", values_to = "count") |>
  mutate(
    group = case_when(
      sample %in% control_cols    ~ "Control",
      sample %in% experiment_cols ~ "Experiment",
      TRUE ~ "Unknown"
    ),
    log2_count = log2(count + 1)   # +1 pseudo-count to handle zeros
  )

# ── 4. Histogram, faceted by group ─────────────────────────────────────────────
ggplot(long_counts, aes(x = log2_count, fill = sample)) +
  geom_histogram(
    aes(y = after_stat(density)),   # rescale counts to density so samples are comparable
    position = "identity",
    alpha = 0.4,
    bins  = 50,
    colour = "grey30",
    linewidth = 0.1
  ) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    name   = expression(log[2](count + 1)),
    limits = c(0, NA)
  ) +
  scale_y_continuous(name = "Density") +
  labs(
    title = "Read Count Histogram",
    fill  = "Sample"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# ── 5. Density plot, faceted by group ──────────────────────────────────────────
ggplot(long_counts, aes(x = log2_count, colour = sample, fill = sample)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    name   = expression(log[2](count + 1)),
    limits = c(0, NA)
  ) +
  scale_y_continuous(name = "Density") +
  labs(
    title  = "Read Count Density",
    colour = "Sample",
    fill   = "Sample"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# ── 6. Combined: histogram + density overlay, faceted by group ─────────────────
ggplot(long_counts, aes(x = log2_count)) +
  geom_histogram(
    aes(y = after_stat(density), fill = sample),
    position  = "identity",
    alpha     = 0.35,
    bins      = 50,
    colour    = "grey30",
    linewidth = 0.1
  ) +
  geom_density(aes(colour = sample), linewidth = 0.9, fill = NA) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    name   = expression(log[2](count + 1)),
    limits = c(0, NA)
  ) +
  scale_y_continuous(name = "Density") +
  labs(
    title  = "Read Count Distribution (Histogram + Density)",
    fill   = "Sample",
    colour = "Sample"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# ── 7. Side-by-side panel: histogram and density next to each other ────────────
library(patchwork)  # install.packages("patchwork") if needed

p_hist <- ggplot(long_counts, aes(x = log2_count, fill = sample)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity",
                 alpha = 0.4, bins = 50, colour = "grey30", linewidth = 0.1) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(title = "Histogram", x = expression(log[2](count + 1)), y = "Density") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"), legend.position = "none")

p_dens <- ggplot(long_counts, aes(x = log2_count, colour = sample, fill = sample)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(title = "Density", x = expression(log[2](count + 1)), y = "Density") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"))

p_hist + p_dens

# ── 8. Save all plots explicitly ────────────────────────────────────────────
p_histogram <- ggplot(long_counts, aes(x = log2_count, fill = sample)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity",
                 alpha = 0.4, bins = 50, colour = "grey30", linewidth = 0.1) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(title = "Read Count Histogram", x = expression(log[2](count + 1)), y = "Density", fill = "Sample") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"))

p_density <- ggplot(long_counts, aes(x = log2_count, colour = sample, fill = sample)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(title = "Read Count Density", x = expression(log[2](count + 1)), y = "Density", colour = "Sample", fill = "Sample") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"))

p_combined <- ggplot(long_counts, aes(x = log2_count)) +
  geom_histogram(aes(y = after_stat(density), fill = sample), position = "identity",
                 alpha = 0.35, bins = 50, colour = "grey30", linewidth = 0.1) +
  geom_density(aes(colour = sample), linewidth = 0.9, fill = NA) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(title = "Read Count Distribution (Histogram + Density)",
       x = expression(log[2](count + 1)), y = "Density", fill = "Sample", colour = "Sample") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"))

ggsave("outputs/histogram_by_group.png", plot = p_histogram, width = 8, height = 6, dpi = 300)
ggsave("outputs/density_by_group.png",   plot = p_density,   width = 8, height = 6, dpi = 300)
ggsave("outputs/combined_by_group.png",  plot = p_combined,  width = 8, height = 6, dpi = 300)

