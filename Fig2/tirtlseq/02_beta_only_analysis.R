setwd("/Volumes/Active/Isabel_Risch/MP028")
library(openxlsx)
library(dplyr)

pathogenicCDR3bSequence <- "^.{4,5}[A|G|D|S|E|Q|R][F|R|L|Y|H|V|I|T|S|M][Y|F]ST.{5}" # "^.{6,7}[Y|F]ST.{5}"


# read in all the samples
all_samples_beta <- list()
input_dir <- "data/inputs/SRR011138/beta_chain_only_data/"
sample_list <- list.files(input_dir)
sample_list <- sapply(strsplit(sample_list, "_"), "[", 1) %>% unique()
for(i in sample_list) {
  all_samples_beta[[i]] <- read.delim(gzfile(paste0(input_dir, i, "_pseudobulk_TRB.tsv.gz")))
  all_samples_beta[[i]]$sample <- i
}
input_dir <- "data/inputs/SRR009265/samples/"
sample_list <- list.files(input_dir)
for(i in sample_list) {
  all_samples_beta[[paste0(i, "_2")]] <- read.delim(gzfile(paste0(input_dir, i,"/", i, "_pseudobulk_TRB.tsv.gz")))
  all_samples_beta[[paste0(i, "_2")]]$sample <- paste0(i, "_2")
}

## collapse everything into one dataframe
all_samples_beta <- dplyr::bind_rows(all_samples_beta)

## Annotate likely pathogenic cells
all_samples_beta$pathogenic1 <- ifelse(grepl(pathogenicCDR3bSequence, all_samples_beta$aaSeqCDR3), TRUE, FALSE)
all_samples_beta$pathogenic2 <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", all_samples_beta$v), all_samples_beta$pathogenic1, FALSE)

## turns out ID01 and ID08_2 are the same sample, merge those, just rename them ID01
all_samples_beta$sample <- ifelse(all_samples_beta$sample %in% "IDO8_2",
                                  "IDO1", all_samples_beta$sample)

## manually input the healthy vs. disease patients
diagnosis <- t(data.frame(`IDO1` = "AI", # this is a duplicate of ID08_2
                          `IDO2` = "HC",
                          `IDO3` = "HC",
                          `IDO4` = "HC",
                          `IDO5` = "AI",
                          `IDO6` = "HC",
                          `IDO7` = "HC",
                          `IDO8` = "HC",
                          `IDO9` = "HC",
                          `IDO10` = "AI",
                          `IDO11` = "AI",
                          `IDO12` = "AI",
                          `IDO6_2` = "AI",
                          `IDO7_2` = "HC"))
colnames(diagnosis) <- "diagnosis"

all_samples_beta <- merge(all_samples_beta, diagnosis, by.x="sample", by.y="row.names")
all_samples_beta <- all_samples_beta[!is.na(all_samples_beta$sample),]

## Re-derive the read fraction for each sample (because we have merged two samples together)
all_samples_beta_rederived <- all_samples_beta %>%
  group_by(sample) %>%
  mutate(readFraction_derived = readCount / sum(readCount)) %>%
  ungroup()

all_samples_beta_rederived$readFraction <- all_samples_beta_rederived$readFraction_derived
all_samples_beta_rederived$readFraction_derived <- NULL

## save point
write.csv(all_samples_beta_rederived, file = "data/inputs/all_samples_BETA_ONLY.csv", row.names = F)


# initialize excel notebook to save results
wb <- createWorkbook()

## get a plot showing percentage of pathogenic cells per patient per location
for (i in 1:2) {
  all_samples_beta_rederived$pathogenic <- all_samples_beta_rederived[[paste0("pathogenic",i)]]
  plot <- all_samples_beta_rederived %>% 
    group_by(sample, pathogenic, diagnosis) %>% 
    summarise(overallPathProportion = sum(readFraction))
  plot <- plot[!plot$pathogenic,]
  plot$overallPathProportion <- signif(plot$overallPathProportion, 6)
  plot$overallPathProportion <- 1 - plot$overallPathProportion
  plot$overallPathProportion <- 100 * plot$overallPathProportion
  plot$pathogenic <- NULL
  
  sheetName <- ifelse(i == 1, "CDR3b", "CDR3b + TRBV")
  addWorksheet(wb, sheetName)
  writeDataTable(wb, sheetName, x = plot, rowNames = F)

}

saveWorkbook(wb, "data/outputs/betaOnly_table_pathogenic_percentage_by_subject.xlsx", overwrite = TRUE)



## UPDATE Jul 22 2026: plot % of sample with the selected TRBV genes only (no CDR3b)
all_samples_beta$selectedTrbv <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", all_samples_beta$v), TRUE, FALSE)

# initialize excel notebook to save TRBV results
wb <- createWorkbook()

## get a plot showing percentage of pathogenic cells per patient per location
plot <- all_samples_beta %>% 
  group_by(sample, selectedTrbv) %>% 
  summarise(overallPathProportion = sum(readFraction))
plot <- plot[plot$selectedTrbv,]
plot$overallPathProportion <- signif(plot$overallPathProportion, 6)
plot$overallPathProportion <- 100 * plot$overallPathProportion
plot$selectedTrbv <- NULL

plot <- merge(plot, diagnosis, by.x="sample", by.y="row.names")
colnames(plot) <- c("Sample", "TRBV9/5-5/5-4 Proportion", "Diagnosis")
sheetName <- "pathogenic TRBV proportion"
addWorksheet(wb, sheetName)
writeDataTable(wb, sheetName, x = plot, rowNames = F)

saveWorkbook(wb, "data/outputs/betaOnly_table_TRBV_percentage_by_subject.xlsx", overwrite = TRUE)

