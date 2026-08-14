library(openxlsx)

setwd("/Volumes/Active/Isabel_Risch/MP028")
paired_1 <- read.delim("data/inputs/SRR011138/all_samples_combined.tsv")
paired_2 <- read.delim("data/inputs/SRR009265/combined/all_samples_combined.tsv")
paired_2$sample <- paste0(paired_2$sample, "_2")
paired <- rbind(paired_1, paired_2)

## Deduplicate the matrix so the two pairing algorithms don't produce redundant clones-- give preference to T-SHELL
##
## !!! REVIEW (NB 2026-07-22 — CORRECTNESS BUG, CORRECTED BELOW !!!
## Original line deduplicated on `alpha_beta` ALONE, across the ENTIRE combined table (all patients),
## and it runs BEFORE the per-sample grouping. `duplicated()` keeps only the first occurrence of each
## clone value globally, so any alpha/beta pairing shared by more than one patient survives in exactly
## ONE sample and is silently deleted from every other sample.
##
## Original line, preserved for provenance:
##   paired <- paired[!duplicated(paired$alpha_beta),]
paired <- paired[order(paired$method, decreasing = T),]
paired <- paired[!duplicated(paired[, c("sample", "alpha_beta")]),]   # CORRECTED: dedup within sample

## turns out ID01 and ID08_2 are the same sample, merge those, just rename them ID01
paired$sample <- ifelse(paired$sample %in% "IDO8_2",
                        "IDO1", paired$sample)

# data frame with the disease/healthy information, taken from email communications
diagnosis <- t(data.frame(`IDO1` = "AI", # this sample is a duplicate of ID08_2
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

## Annotate likely pathogenic cells in a successive manner, starting with just CDR3b, then adding a TRBV criterion, then TRAV21
pathogenicCDR3bSequence <- "^.{4,5}[A|G|D|S|E|Q|R][F|R|L|Y|H|V|I|T|S|M][Y|F]ST.{5}" # "^.{6,7}[Y|F]ST.{5}"
paired$pathogenic1 <- ifelse(grepl(pathogenicCDR3bSequence, paired$cdr3b), TRUE, FALSE)
paired$pathogenic2 <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", paired$vb), paired$pathogenic1, FALSE)
paired$pathogenic3 <- ifelse(grepl("TRAV21", paired$va), paired$pathogenic2, FALSE)

paired <- merge(paired, diagnosis, by.x="sample", by.y="row.names")
# remove the sample that has no label (ID01), that's a duplicate of the original (ID08_2)
paired <- paired[!is.na(paired$sample),]

# save the entire paired object so we can come back to it later
write.csv(paired, file = "data/inputs/all_samples_paired.csv", row.names = F)

# initialize excel notebook to save results
wb <- createWorkbook()

## calculate percentage of pathogenic cells per patient using our various pathogenic criteria 
for (i in 1:3) {
  paired$pathogenic <- paired[[paste0("pathogenic",i)]]
  plot <- paired %>% 
    group_by(sample, pathogenic, diagnosis) %>% 
    summarise(totalNumber = sum(wij)) %>% 
    group_by(sample) %>%
    mutate(totalDenominator = sum(totalNumber)) %>% 
    mutate(proportion = totalNumber / totalDenominator)
  plot <- plot[!plot$pathogenic,]
  plot$proportion <- signif(plot$proportion, 6)
  plot$proportion <- 1 - plot$proportion
  plot$proportion <- 100 * plot$proportion
  plot$pathogenic <- NULL
  
  sheetName <- ifelse(i == 1, "CDR3b", ifelse(i == 2, "CDR3b + TRBV",  "CDR3b + TRBV + TRAV"))
  addWorksheet(wb, sheetName)
  writeDataTable(wb, sheetName, x = plot, rowNames = F)
  
}

saveWorkbook(wb, "data/outputs/table_pathogenic_percentage_by_subject_PAIRED.xlsx", overwrite = TRUE)




