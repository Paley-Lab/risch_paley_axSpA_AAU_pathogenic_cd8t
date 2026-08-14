.libPaths(c("/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/R_libraries/4.4.0", 
            "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

setwd("/storage1/fs1/paleym/Active/Isabel_Risch/MP002/")

# library(Seurat)
library(openxlsx)

datadir <- "/storage1/fs1/paleym/Active/SeqData/Omniscope/"
# datadir <- 'nick_analysis/2026_nick_outputs_from_GoogleDrive/SequencingRuns/'

pathogenicCDR3bSequence <- "^.{3,4}[A|G|D|S|E|Q|R][F|R|L|Y|H|V|I|T|S|M][Y|F]ST.{4}" # "^.{5,6}[Y|F]ST.{4}"


####### ------- Make our beta-only and paired dataframes ------- #######
contig.files <- list.files(datadir, paste0(".csv$", collapse = "|"), recursive = TRUE, full.names = TRUE)
df <- lapply(contig.files, read.csv) 

# define a filtering function for the paired data. It only allows each cell to 
# have one alpha and one beta, just for simplicity.
filter_paired_tcr_base <- function(df) {
  counts <- table(df$cell_id[df$locus %in% c("TRA", "TRB")],
                  df$locus[df$locus %in% c("TRA", "TRB")])
  
  valid_cells <- rownames(counts)[counts[, "TRA"] == 1 & counts[, "TRB"] == 1]
  
  df[df$cell_id %in% valid_cells, ]
}

beta <- list()
alpha <- list()
paired <- list()
for (i in seq_along(df)) {
  # only consider the beta sequences, filter the alphas
  betaRows <- grep("B", df[[i]]$locus)
  beta[[i]] <- df[[i]][betaRows,]
  # add a column that corresponds to the sample
  beta[[i]]$sample <- sapply(strsplit(contig.files, "//|\\."), "[", 2)[i]
  
  # only consider the beta sequences, filter the alphas
  alpha[[i]] <- df[[i]][-betaRows,]
  # add a column that corresponds to the sample
  alpha[[i]]$sample <- sapply(strsplit(contig.files, "//|\\."), "[", 2)[i]
  
  # for the paired data, only take stuff with ONE alpha and ONE beta
  paired[[i]] <- filter_paired_tcr_base(df[[i]])
  paired[[i]]$sample <- sapply(strsplit(contig.files, "//|\\."), "[", 2)[i]
}

# merge this into one table
beta <- dplyr::bind_rows(beta)
alpha <- dplyr::bind_rows(alpha)
paired <- dplyr::bind_rows(paired)

# save those objects, it took a while
saveRDS(beta, file = "data/beta_only_chains.Rds")
saveRDS(alpha, file = "data/alpha_only_chains.Rds")


# Make some edits to the paired table to get alpha and beta chains on same row
tra <- paired[grep("A", paired$locus),]
trb <- paired[grep("B", paired$locus),]
colnames(tra) <- paste0(colnames(tra), "_alpha")
colnames(trb) <- paste0(colnames(trb), "_beta")
paired_reshaped <- merge(tra, trb, by=c(1,67), all=T)
colnames(paired_reshaped)[1] <- "cell_id"
colnames(paired_reshaped)[2] <- "sample"
paired_reshaped$mean_umi_count <- (paired_reshaped$umi_count_alpha + paired_reshaped$umi_count_beta)/2

saveRDS(paired_reshaped, file = "data/paired_object.Rds")

####### ------- Annotate pathogenic cells ------- #######

beta <- readRDS("data/beta_only_chains.Rds")
paired_reshaped <- readRDS("data/paired_object.Rds")


# CDR3b-only requirement (beta only)
beta$pathogenic_cdr3b <- ifelse(grepl(pathogenicCDR3bSequence, beta$cdr3_aa), TRUE, FALSE)

# CDR3b plus TRBV requirement (beta only object)
beta$pathogenic_cdr3b_trbv <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", beta$v_call), beta$pathogenic_cdr3b, FALSE)

# CDR3b plus TRBV plus TRAV requirement on paired object
paired_reshaped$pathogenic_cdr3b <- grepl(pathogenicCDR3bSequence, paired_reshaped$cdr3_aa_beta)
paired_reshaped$pathogenic_cdr3b_trbv <- grepl("TRBV9|TRBV5-5|TRBV5-4", paired_reshaped$v_call_beta) & paired_reshaped$pathogenic_cdr3b
paired_reshaped$pathogenic_cdr3b_trbv_trav21 <- grepl("TRAV21", paired_reshaped$v_call_alpha) & paired_reshaped$pathogenic_cdr3b_trbv


# define a really hard-coded function to get the pathogenic %age for both of our pathogenic columns
summarizePathogenicInfo <- function(dataframe=NULL, pathogenicColumn=NULL, mode="beta") {
  total <- dataframe
  total$pathogenic <- total[[pathogenicColumn]]
  # get the percentage of reads that are pathogenic per sample
  if (mode=="beta") {
    pathogenic_pct <- total %>% 
      group_by(sample, pathogenic) %>% 
      summarise("count" = sum(umi_count))
  } else if (mode=="paired") {
    pathogenic_pct <- total %>% 
      group_by(sample, pathogenic) %>% 
      summarise("count" = sum(mean_umi_count))
  }
  
  for (i in unique(pathogenic_pct$sample)) {
    pathogenic_true_row <- pathogenic_pct[(pathogenic_pct$sample==i) & (pathogenic_pct$pathogenic),]
    if(nrow(pathogenic_true_row) < 1) {
      pathogenic_pct <- rbind(pathogenic_pct, data.frame("sample" = i, "pathogenic" = TRUE, "count"= 0))
    }
  }  
  
  pathogenic_pct <- pathogenic_pct %>%
    group_by(sample, pathogenic) %>% 
    summarise(totalNumber = sum(count)) %>% 
    mutate(totalDenominator = sum(totalNumber)) %>% 
    mutate(proportion = totalNumber / totalDenominator)
  
  ## get a plot showing percentage of pathogenic cells per patient per location
  plot <- pathogenic_pct[pathogenic_pct$pathogenic,]
  plot$proportion <- signif(plot$proportion, 6)
  plot$proportion <- 100 * plot$proportion
  plot$pathogenic <- NULL
  
  # add in the metadata about YeiH binding 
  meta <- read.delim("nick_analysis/omniscope_metadata.tsv")
  
  plot$sample <- toupper(plot$sample) #make sample names uppercase to match meta
  merged <- merge(plot[c(1,4)], meta[,c(1,14,16)], by=1, all = T)
  
  return(merged)
}

# run this for both definitions of pathogenicity
plot_cdr3b <- summarizePathogenicInfo(beta, "pathogenic_cdr3b", mode="beta")
plot_cdr3b_trbv <- summarizePathogenicInfo(beta, "pathogenic_cdr3b_trbv", mode="beta")
plot_cdr3b_trbv_trav <- summarizePathogenicInfo(paired_reshaped, "pathogenic_cdr3b_trbv_trav21", mode="paired")

# initialize excel notebook to save results
wb <- createWorkbook()
addWorksheet(wb, "cdr3b")
writeDataTable(wb, "cdr3b", x = plot_cdr3b, rowNames = F)
addWorksheet(wb, "cdr3b_trbv")
writeDataTable(wb, "cdr3b_trbv", x = plot_cdr3b_trbv, rowNames = F)
addWorksheet(wb, "cdr3b_trbv_trav")
writeDataTable(wb, "cdr3b_trbv_trav", x = plot_cdr3b_trbv_trav, rowNames = F)
saveWorkbook(wb, "data/outputs/omniscope_pathogenic_proportions_allAnalyses.xlsx", overwrite = TRUE)



## ---------- UPDATE Jul 22 2026: plot % of sample with the selected TRBV genes only (no CDR3b) ----------- ##
beta <- readRDS("data/beta_only_chains.Rds")

beta$selectedTrbv <- ifelse(grepl("TRBV9|TRBV5-5|TRBV5-4", beta$v_call), TRUE, FALSE)

# initialize excel notebook to save TRBV results
wb <- createWorkbook()

## get a plot showing percentage of pathogenic cells per patient per location
plot <- beta %>% 
  group_by(sample, selectedTrbv) %>% 
  summarise("count" = sum(umi_count)) %>%
  group_by(sample) %>% 
  mutate(totalDenominator = sum(count)) %>% 
  mutate(proportion = count / totalDenominator)

plot <- plot[plot$selectedTrbv,]
plot$proportion <- signif(plot$proportion, 6)
plot$proportion <- 100 * plot$proportion
plot$selectedTrbv <- NULL

plot$Diagnosis <- ifelse(grepl("hc", plot$sample), "HC", "AI")
plot <- plot[,c('sample', 'proportion', 'Diagnosis')]
colnames(plot) <- c("Sample", "TRBV9/5-5/5-4 Percentage", "Diagnosis")
sheetName <- "pathogenic TRBV proportion"
addWorksheet(wb, sheetName)
writeDataTable(wb, sheetName, x = plot, rowNames = F)

saveWorkbook(wb, "data/outputs/betaOnly_table_TRBV_percentage_by_subject.xlsx", overwrite = TRUE)


