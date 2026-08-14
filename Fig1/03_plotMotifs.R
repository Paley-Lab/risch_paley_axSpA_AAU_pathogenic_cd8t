setwd("/Volumes/Active/Isabel_Risch/MP033/")
library(ggseqlogo)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggalluvial)
library(dplyr)


# ---- read in the outputs from the previous analysis ------------------------------------------
# allAntigensResults <- read.csv("outputs/differential_selection_results/allAntigens_enrichRatio_Results.csv")
logEnrichScores <- read.csv("outputs/log2_enrichment_ratios.csv", row.names = 1)
allResults <- readRDS("outputs/differential_selection_results/allSignificanceResults.Rds")
allResults_combined <- allResults[!(names(allResults) %in% "allAntigens")]
for (i in 1:length(allResults_combined)) {
  allResults_combined[[i]]$Antigen <- names(allResults_combined)[i]
}
allResults_combined <- dplyr::bind_rows(allResults_combined)
write.csv(allResults_combined, file = "outputs/all_enrichment_results.csv", row.names = F)

# ---- Output significant results for a SeqLogo plot...
enriched_up <- read.csv("outputs/differential_selection_results/plots/yeiH_union_motifs_noPSG5.csv")
enriched_up <- enriched_up[grep("YeiH", enriched_up$found_in),]
enriched_up <- enriched_up[enriched_up$found_in != "YeiH",]

png(filename = "outputs/differential_selection_results/plots/ggseqlogo_enriched_overlaps.png", 
    width = 2, height = 2.5,units = "in", res = 300)
ggseqlogo(enriched_up$triplet, method="bits")
dev.off()

new_preYST_motif <- c(paste(unique(stringr::str_sub(enriched_up$triplet, 1, 1)), collapse = "|"), 
               paste(unique(stringr::str_sub(enriched_up$triplet, 2, 2)), collapse = "|"), 
               paste(unique(stringr::str_sub(enriched_up$triplet, 3, 3)), collapse = "|"))

new_motif <- paste0("^.{4,5}[", 
                    new_preYST_motif[2],"][", 
                    new_preYST_motif[3],
                    "][Y|F]ST.{5}")
print(new_motif)


# enriched_up_allAnt <- allAntigensResults[(allAntigensResults$padj < padj_thre) & 
#                                     (allAntigensResults$avg_log2FC > logFC_thre),]
# ggseqlogo(enriched_up_allAnt$triplet, method="bits") 
# 
# 
# enriched_down <- allAntigensResults[(allAntigensResults$padj < padj_thre) & 
#                                       (allAntigensResults$avg_log2FC < -logFC_thre),]
# ggseqlogo(enriched_down$triplet, method="bits") + scale_y_continuous(limits = c(0,2))



# ---- ...and a heatmap of amino acid proportion per position

plotAAproportions <- function(heatmap) {
  AA_ORDER <- c('*', 'G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y',
                'N', 'Q', 'D', 'E', 'K', 'R', 'H')
  heatmap <- heatmap[,-ncol(heatmap)] 
  heatmap <- heatmap[,-grep("eptid", colnames(heatmap))] # get rid of the no peptide control column
  heatmap <- t(heatmap) %>% melt()
  
  heatmap$pos1 <- stringr::str_sub(heatmap[[2]], start = 1, end = 1)
  heatmap$pos2 <- stringr::str_sub(heatmap[[2]], start = 2, end = 2)
  heatmap$pos3 <- stringr::str_sub(heatmap[[2]], start = 3, end = 3)
  
  pos1 <- heatmap %>% group_by(pos1) %>% summarise(#`score` = mean(value), 
    `proportion` = n())
  pos1$proportion <- pos1$proportion/sum(pos1$proportion)
  
  pos2 <- heatmap %>% group_by(pos2) %>% summarise(#`score` = sum(value),
    `proportion` = n())
  pos2$proportion <- pos2$proportion/sum(pos2$proportion)
  
  pos3 <- heatmap %>% group_by(pos3) %>% summarise(#`score` = sum(value),
    `proportion` = n())
  pos3$proportion <- pos3$proportion/sum(pos3$proportion)
  
  
  totalHeatmap <- merge(as.data.frame(AA_ORDER), pos1, all = T, by = 1)
  colnames(totalHeatmap)[2] <- "pos1"
  totalHeatmap <- merge(totalHeatmap, pos2, by=1, all = T)
  colnames(totalHeatmap)[3] <- "pos2"
  totalHeatmap <- merge(totalHeatmap, pos3, by=1, all = T)
  colnames(totalHeatmap)[4] <- "pos3"
  
  totalHeatmap[is.na(totalHeatmap)] <- 0
  
  # merge() always re-sorts by the join column - restore the custom AA_ORDER
  totalHeatmap <- totalHeatmap[match(AA_ORDER, totalHeatmap$AA_ORDER), ]
  
  return(totalHeatmap)
}

plotAlluvial <- function(alluvial) {
  alluvial$pos1 <- stringr::str_sub(alluvial[["triplet"]], start = 1, end = 1) %>% unlist()
  alluvial$pos2 <- stringr::str_sub(alluvial[["triplet"]], start = 2, end = 2) %>% unlist()
  alluvial$pos3 <- stringr::str_sub(alluvial[["triplet"]], start = 3, end = 3) %>% unlist()
  
  # Count combinations of pos1/pos2/pos3 (needed for alluvial "freq")
  df_counts <- alluvial %>% group_by(pos1, pos2, pos3) %>%
    summarise("Freq" = n())
  
  # Convert to long ("lodes") format for ggalluvial
  df_long <- ggalluvial::to_lodes_form(df_counts,
                           axes = c("pos1", "pos2", "pos3"),
                           id = "Cohort")
  
  # Build the alluvial plot
  p <- ggplot(df_long, aes(x = x, stratum = stratum, alluvium = Cohort,
                           y = Freq, fill = stratum, label = stratum)) +
    geom_alluvium(alpha = 0.7) +
    geom_stratum() +
    geom_text(stat = "stratum", size = 3) +
    scale_x_discrete(limits = c("pos1", "pos2", "pos3"), expand = c(.05, .05)) +
    labs(x = "Position", y = "Count", fill = "Amino Acid") +
    theme_minimal()
  
  # Save as PDF
  ggsave("outputs/alluvial_plot.pdf", plot = p, width = 8, height = 6)
  
  return(p)
  
}

plotAlluvial(enriched_up)

heatmap_up <- plotAAproportions(logEnrichScores[logEnrichScores$aa %in% enriched_up$triplet,])
# heatmap_down <- plotAAproportions(logEnrichScores[logEnrichScores$aa %in% enriched_down$triplet,])

write.csv(heatmap_up, file = "outputs/differential_selection_results/plots/heatmap_enriched_up.csv", row.names = F)
# write.csv(heatmap_down, file = "outputs/heatmap_enriched_down.csv", row.names = F)


plot_heatmap <- function(df, plot_title, colorscale="darkred", limits=c(0, 0.5)) {
  mat <- as.matrix(df)
  
  eps <- 1e-6
  lo <- limits[1]
  hi <- limits[2]
  
  ramp_colors <- colorRampPalette(c("white", colorscale))(99)
  colors <- c("grey80", ramp_colors)                    # bin 1 = exact 0 -> light grey
  breaks <- c(-eps, seq(max(eps, lo), hi, length.out = 100))  # bin 1 spans (-eps, eps], catches 0
  
  print(pheatmap(
    mat,
    color = colors,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = as.character(plot_title),
    border_color = "grey20"
  ))
}


df_up <- heatmap_up
rownames(df_up) <- df_up[[1]]
df_up <- df_up[,-1]

# df_dn <- heatmap_down
# rownames(df_dn) <- df_dn[[1]]
# df_dn <- df_dn[,-1]

pdf(file = "outputs/differential_selection_results/plots/heatmap_enriched_overlaps.pdf", 
    width = 2, height = 3)
plot_heatmap(df_up,"Enriched UP", limits = c(0,0.3))
dev.off()

# plot_heatmap(df_dn, "Enriched DOWN", colorscale = "navy")
