
setwd("/Volumes/Active/Isabel_Risch/MP021")

library(Seurat)
library(VISION)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)

dds <- readRDS("/Volumes/Active/Isabel_Risch/MP021/data/de/DESeq2/pathogeniccluster_vs_MAIT_subject_tissue/final_dds.Rds")
res <- results(dds, contrast = c("cluster", "pathogenic_CD8T", "MAIT"))
de <- as.data.frame(res)


#########################################
########## Fisher's Exact Test ##########
#########################################

# read in the genesets
gs2gene <- list("H" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/h.all.v2023.2.Hs.symbols.gmt"),
                "KEGG" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"),
                "Reactome" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/c2.cp.reactome.v2023.2.Hs.symbols.gmt"),
                "CP" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/c2.cp.v2023.2.Hs.symbols.gmt"),
                "GOBP" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/c5.go.bp.v2023.2.Hs.symbols.gmt"),
                "GOCC" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/c5.go.cc.v2023.2.Hs.symbols.gmt"),
                "GOMF" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/c5.go.mf.v2023.2.Hs.symbols.gmt"),
                "GO" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/c5.go.v2023.2.Hs.symbols.gmt"),
                "C8" = qusage::read.gmt("~/Documents/resources/MSigDB/genesets/c8.all.v2023.2.Hs.symbols.gmt")
)

de_up <- de[de$log2FoldChange > 1.5 & de$padj < 0.05,]
de_up <- rownames(de_up)
de_dn <- de[de$log2FoldChange < -1.5 & de$padj < 0.05,]
de_dn <- rownames(de_dn)

# Use the Jiyang Yu Lab's package NetBID2 
res_up <- NetBID2::funcEnrich.Fisher(input_list=de_up, bg_list=rownames(de), Pv_thre=0.5,Pv_adj = 'BH', gs2gene = gs2gene, 
                                     use_gs = "H", min_gs_size = 10, max_gs_size = 200)
res_dn <- NetBID2::funcEnrich.Fisher(input_list=de_dn, bg_list=rownames(de), Pv_thre=0.5,Pv_adj = 'BH', gs2gene = gs2gene, 
                                     use_gs = "H", min_gs_size = 10, max_gs_size = 200)
NetBID2::draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=10,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,Pv_thre = 0.5,
                                 # pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
                                 cluster_gs=TRUE,cluster_gene = TRUE,h=0.95, Pv_col = "Adj_P")


res_up <- NetBID2::funcEnrich.Fisher(input_list=de_up, bg_list=rownames(de), Pv_thre=0.5,Pv_adj = 'BH', gs2gene = gs2gene, 
                                     use_gs = "GO", min_gs_size = 10, max_gs_size = 200)
res_dn <- NetBID2::funcEnrich.Fisher(input_list=de_dn, bg_list=rownames(de), Pv_thre=0.5,Pv_adj = 'BH', gs2gene = gs2gene, 
                                     use_gs = "GO", min_gs_size = 10, max_gs_size = 200)
NetBID2::draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=10,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,Pv_thre = 0.5,
                                 # pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
                                 cluster_gs=TRUE,cluster_gene = TRUE,h=0.95, Pv_col = "Adj_P")


res_up <- NetBID2::funcEnrich.Fisher(input_list=de_up, bg_list=rownames(de), Pv_thre=0.5,Pv_adj = 'BH', gs2gene = gs2gene, 
                                     use_gs = "KEGG", min_gs_size = 10, max_gs_size = 200)
res_dn <- NetBID2::funcEnrich.Fisher(input_list=de_dn, bg_list=rownames(de), Pv_thre=0.5,Pv_adj = 'BH', gs2gene = gs2gene, 
                                     use_gs = "KEGG", min_gs_size = 10, max_gs_size = 200)
NetBID2::draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=10,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,Pv_thre = 0.5,
                                 # pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
                                 cluster_gs=TRUE,cluster_gene = TRUE,h=0.95, Pv_col = "Adj_P")

res_up <- NetBID2::funcEnrich.Fisher(input_list=de_up, bg_list=rownames(de), Pv_thre=0.05,Pv_adj = 'BH', gs2gene = gs2gene, 
                                     use_gs = "C8", min_gs_size = 10, max_gs_size = 200)
res_dn <- NetBID2::funcEnrich.Fisher(input_list=de_dn, bg_list=rownames(de), Pv_thre=0.05,Pv_adj = 'BH', gs2gene = gs2gene, 
                                     use_gs = "C8", min_gs_size = 10, max_gs_size = 200)
NetBID2::draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=10,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,Pv_thre = 0.5,
                                 # pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
                                 cluster_gs=TRUE,cluster_gene = TRUE,h=0.95, Pv_col = "Adj_P")
