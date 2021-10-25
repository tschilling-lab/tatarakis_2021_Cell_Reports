library(SingleCellExperiment)
library(cellassign)
library(Seurat)
library(scran)
library(loomR)
setwd("/home/cang/Dropbox/Projects_UCI/ZebraFishNeuralCrest/dev/publication")

# CELLASSIGN FOR WT
zebraNC <- readRDS(file = './data/WT/Arch1_new_aggregate_filtered_noNT.RDS')
zebraNC.sce <- as.SingleCellExperiment(zebraNC)
zebraNC.sce <- computeSumFactors(zebraNC.sce)
s <- sizeFactors(zebraNC.sce)
earlyNC.genes <- as.character( read.csv('./data/gene_lists/early_NC_genes.csv')[,2] )
pigment.genes <- as.character( read.csv('./data/gene_lists/pigment_genes.csv')[,2] )
skeletal.genes <- as.character( read.csv('./data/gene_lists/skeletal_genes.csv')[,2] )
marker_gene_list <- list(
  earlyNC = earlyNC.genes,
  neural_glial = c("neurog1","ascl1a","ascl1b","eya1","six1b","phox2bb","phox2a","gap43","pou3f1","sox2","egr2b","mbpb","s100b","foxd3","sox10","wnt11r","her15.2","her4.1","her2","her4.2"),
  pigment = pigment.genes,
  skeletal = skeletal.genes
)
marker_gene_mat <- marker_list_to_mat(marker_gene_list)
fit <- cellassign(exprs_obj = zebraNC.sce[rownames(marker_gene_mat),], 
                  marker_gene_info = marker_gene_mat, 
                  s = s, 
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)
pheatmap::pheatmap(cellprobs(fit))
write.csv(cellprobs(fit), file="./data/cellassign/cellprobs.csv")
write.csv(celltypes(fit), file="./data/cellassign/celltypes.csv")