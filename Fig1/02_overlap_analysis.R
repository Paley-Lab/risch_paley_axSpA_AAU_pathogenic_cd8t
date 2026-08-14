## =============================================================================
## Differential Motif Selection: Volcano Plots + YeiH/GPER1/PRPF3/RNASEH2B Venn
## =============================================================================
## Input:  <Condition>_enrichRatio_Results.csv with columns:
##         triplet, p.val, avg_log2FC, n_obs, padj
## Output: plots/volcano_plots_all_conditions.pdf, plots/overlap_venn.pdf,
##         plots/overlap_upset.pdf, plots/yeiH_union_motifs_noPSG5.csv
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------

# install.packages(c("tidyverse", "ggrepel", "UpSetR", "VennDiagram"))  # run once if needed
library(tidyverse)
library(ggrepel)
library(UpSetR)
library(VennDiagram)

setwd("/Volumes/Active/Isabel_Risch/MP033/outputs/differential_selection_results/")
dir.create("plots", showWarnings = FALSE)

PADJ_CUT  <- 0.05
LFC_CUT   <- 1

# Every condition file to read in. "AllAntigens" is kept separate from the
# conditions used in the YeiH Venn diagram. PSG5 excluded (per current design)
# so that exactly 4 sets remain for a classic 4-set Venn: YeiH, GPER1, PRPF3,
# RNASEH2B.
files <- c(
  GPER1       = "GPER1_enrichRatio_Results.csv",
  PRPF3       = "PRPF3_enrichRatio_Results.csv",
  # PSG5      = "PSG5_enrichRatio_Results.csv",
  RNASEH2B    = "RNASEH2B_enrichRatio_Results.csv",
  YeiH        = "YeiH_enrichRatio_Results.csv",
  AllAntigens = "allAntigens_enrichRatio_Results.csv"
)

## ---- 1. Load data -----------------------------------------------------------

res_list <- imap(files, ~ read_csv(.x, show_col_types = FALSE) %>%
                   mutate(condition = .y))

res_all <- bind_rows(res_list)

## ---- 2. Classify significance & volcano plots --------------------------------

res_all <- res_all %>%
  mutate(sig = case_when(
    padj < PADJ_CUT & avg_log2FC >  LFC_CUT ~ "Enriched",
    padj < PADJ_CUT & avg_log2FC < -LFC_CUT ~ "Depleted",
    TRUE                                    ~ "NS"
  ))

sig_colors <- c(Enriched = "firebrick3", Depleted = "steelblue3", NS = "grey75")

make_volcano <- function(df, cond_name) {
  # Label the top few hits in each direction for readability
  top_labels <- df %>%
    filter(sig != "NS") %>%
    arrange(padj) %>%
    group_by(sig) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  ggplot(df, aes(avg_log2FC, -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_labels, aes(label = triplet),
                    size = 2.8, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = sig_colors) +
    labs(title = paste("Volcano plot:", cond_name),
         x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"adj. p-value"),
         color = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

conditions <- names(files)

pdf("plots/volcano_plots_all_conditions.pdf", width = 4, height = 4)
for (cond in conditions) {
  p <- make_volcano(filter(res_all, condition == cond), cond)
  print(p)
}
dev.off()

## ---- 3. Positive-significant motif sets per condition -------------------------

pos_sig_sets <- res_all %>%
  filter(sig == "Enriched") %>%
  group_by(condition) %>%
  summarise(motifs = list(triplet), .groups = "drop") %>%
  deframe()

# quick summary of how many enriched motifs each condition has
enrich_counts <- map_int(pos_sig_sets, length)
print(enrich_counts)

## ---- 4. YeiH + other antigen sets (no pairwise intersection step) -------------

# Drop AllAntigens so exactly 4 sets remain: GPER1, PRPF3, RNASEH2B, YeiH
venn_sets <- pos_sig_sets[names(pos_sig_sets) != "AllAntigens"]

# Union of enriched motifs across those 4 conditions
enriched_union <- reduce(venn_sets, union)

# Count how many of the 4 sets each motif appears in
occurrence_df <- tibble(triplet = enriched_union) %>%
  rowwise() %>%
  mutate(n_sets_present = sum(map_lgl(venn_sets, ~ triplet %in% .x)),
         found_in = paste(names(venn_sets)[map_lgl(venn_sets, ~ triplet %in% .x)],
                          collapse = ";")) %>%
  ungroup() %>%
  arrange(desc(n_sets_present), triplet)

write_csv(occurrence_df, "plots/yeiH_union_motifs_noPSG5.csv")
cat("Union of enriched motifs across YeiH/GPER1/PRPF3/RNASEH2B:", nrow(occurrence_df), "\n")

## ---- 5. Overlap visualizations -------------------------------------------------

## 5a. Classic 4-set Venn diagram (equal-sized ellipses, standard layout).
##     VennDiagram::venn.diagram computes all 15 intersection regions
##     automatically from the 4 input sets.
venn_colors <- c(GPER1 = "#F8D686", PRPF3 = "#74FA79",
                 RNASEH2B = "#8DD4FB", YeiH = "#FF2F93")[names(venn_sets)]

venn_plot <- venn.diagram(
  x = venn_sets,
  category.names = names(venn_sets),
  filename = NULL,
  fill = venn_colors,
  alpha = 0.5,
  col = "black",
  cex = 3,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.cex = 2,
  cat.fontface = "bold",
  main = "Positive-significant motif overlap: YeiH, GPER1, PRPF3, RNASEH2B"
)

pdf("plots/overlap_venn.pdf", width = 7, height = 7)
grid.newpage()
grid.draw(venn_plot)
dev.off()

# VennDiagram litters the working directory with a .log file - clean it up
invisible(file.remove(list.files(pattern = "^VennDiagram.*\\.log$")))

## 5b. UpSet plot of the same 4 sets (complements the Venn - easier to read
##     exact intersection sizes off the bars)
pdf("plots/overlap_upset.pdf", width = 7, height = 5)
print(upset(fromList(venn_sets), order.by = "freq",
            mainbar.y.label = "Motif Intersection Size", sets.x.label = "Enriched Motifs"))
dev.off()

# ---- What % of DAMs are negatively enriched? ------------------------
res <- res_all[-which(res_all$condition %in% "AllAntigens"),]
res_summary <- res %>% group_by(condition, sig) %>%
  summarise(count=n())
sum(res_summary[res_summary$sig=="Depleted", "count"])
# 1749
sum(res_summary[res_summary$sig=="Enriched", "count"])
# 393
sum(res_summary[res_summary$sig=="NS", "count"])
# 34898

sum(res_summary[res_summary$sig=="Depleted", "count"]) / sum(res_summary[res_summary$sig!="NS", "count"])

## ---- 6. Summary --------------------------------------------------------------

cat("\n===== SUMMARY =====\n")
cat("Enriched motif counts per condition:\n"); print(enrich_counts)
cat("\nSize of union across YeiH/GPER1/PRPF3/RNASEH2B:", length(enriched_union), "\n")
cat("\nFiles written: plots/volcano_plots_all_conditions.pdf, plots/overlap_venn.pdf,\n",
    "plots/overlap_upset.pdf, plots/yeiH_union_motifs_noPSG5.csv\n")
