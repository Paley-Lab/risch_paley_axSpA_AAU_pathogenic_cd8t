setwd("/Volumes/Active/Isabel_Risch/MP003")

library(Seurat)
library(scRepertoire)
library(dplyr)
library(turboGliph)

d <- readRDS("data/robjects/9_TCRgenesRemoved_sexRegressed.Rds")
combined.TCR <- readRDS("data/robjects/combined_tcrs.Rds")

###########################################
#### GLIPH v2 analysis of entire dataset #### Mar 10, 2026
###########################################

dir <- paste0("GLIPH2/2026/gliph_entire_dataset_by_condition/")
if(!dir.exists(dir)) dir.create(dir) 

# To use turbo_gliph, gliph2 and gliph_combined, CDR3 beta sequences are required. 
# These must be specified either by a character vector or in a dataframe in a column named CDR3b. 
# Additional information is not required for clustering, but is recommended for automatic cluster scoring. 
# This includes the following information:
# column TRBV : V gene information of the particular sequence
# column patient : index number or similar of the patient from which the clone was isolated
# column counts : number of clones in the sample
# column HLA : HLA alleles of the particular patient, separated by commas

colnames(all) <- c("CDR3b","TRBV","TRBJ","patient","counts")
all$HLA <- "HLA-B27"
all$disease_status <- ifelse(grepl("isease", all$patient), "Disease", "HC")

res_gliph2 <- turboGliph::gliph2(cdr3_sequences = all, n_cores = 1, cluster_min_size = 3, result_folder = paste0(dir, "/turboGliph"))
saveRDS(res_gliph2, file = paste0(dir, "/turboGliph/gliph2_res.Rds"))

res_gliph2 <- readRDS(paste0(dir, "/turboGliph/gliph2_res.Rds"))
# res_gliph_combined_2 <- turboGliph::gliph_combined(cdr3_sequences = all, n_cores = 1, cluster_min_size = 3, 
#                                                    result_folder = paste0(dir, "/turboGliph"), local_method = 'rrs')

p <- plot_network(clustering_output = res_gliph2, n_cores = 1, 
                  color_info = "disease_status", 
                  color_palette = function(n){return(sample(c( "#189D75", "#FF69B4"), n))})
# plot_network(clustering_output = res_gliph_combined_2, n_cores = 1, color_info = "disease_status")

p[["x"]][["edges"]]$color <- ifelse(p[["x"]][["edges"]]$color %in% "orange", "lightgrey", "skyblue")
p[["x"]][["edges"]]$length <- 100
p[["x"]][["nodes"]]$size <- 25
p

# The TurboGliph plotting function is a little stochastic, so save this particular
# plot for reproducibility
saveRDS(p, file = "data/turboGliphPlot.Rds")

p <- readRDS(file = "data/turboGliphPlot.Rds")
cluster_info <- res_gliph2[["cluster_properties"]]
cluster_info <- cluster_info[order(cluster_info$total.score, decreasing = F),]
cluster_info$cluster_ID <- rownames(cluster_info)
write.csv(cluster_info, row.names = F, file = paste0(dir, "/turboGliph/cluster_summary_info.csv"))

# The plotting function sucks. Export the raw information to be plotted in Cytoscape
write.csv(p[["x"]][["nodes"]], row.names = F, file = paste0(dir, "/turboGliph/network_node_info.csv"))
write.csv(p[["x"]][["edges"]], row.names = F, file = paste0(dir, "/turboGliph/network_edge_info.csv"))

