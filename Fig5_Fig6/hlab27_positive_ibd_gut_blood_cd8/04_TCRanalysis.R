
###########################################################################
########## Getting a dataframe of TCR sequences to check ###########
###########################################################################

# only look at the cells in the KLRB1+ cluster
tmp <- subset(d, idents = "KLRB1_CD8T")

# pull the metadata
meta <- tmp@meta.data

# get the mean pathogenic score for each clone
means <- meta %>%
  group_by(subject, tissue_merged) %>%
  dplyr::summarize(Mean_pathogenic_score = 
                     mean(Core_JCI_MP156_c4_vs_c3_noTCR_20g, na.rm=TRUE))
