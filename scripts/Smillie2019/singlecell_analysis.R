# Script to process single-cell data from Smillie, et al. (2019) Cell

# Load libraries

library(Seurat)
library(Matrix)
library(singlecellmethods)
library(dplyr)
library(harmony)
library(uwot)

# Load data

data <- readMM("/path/gene_sorted-Imm.matrix.mtx")
genes <- read.table("/path/Imm.genes.tsv")
row.names(data) <- genes$V1

cells <- read.table("/path/Imm.barcodes2.tsv")
colnames(data) <- cells$V1

meta_data <- read.table("/path/all.meta2.txt", sep = "\t", header = T)
meta_data <- meta_data[-1,]
row.names(meta_data) <- meta_data$NAME
meta_data$Subject <- factor(meta_data$Subject, levels = levels(meta_data$Subject)[-1])

# Select T cells

meta_data$Tcell = meta_data$Cluster %in% c("CD4+ Activated Fos-hi", "CD4+ Activated Fos-lo", "CD4+ Memory", "CD4+ PD1+", "CD8+ IELs", 
                                             "CD8+ IL17+", "CD8+ LP", "Cycling T", "Tregs")
meta_data <- meta_data[meta_data$Tcell,]
data <- data[,row.names(meta_data)]

# QC cells (> 300 genes, < 20% MT UMIs)

meta_data$percent_mito <- colSums(data[grepl("^MT-", row.names(data)),])/colSums(data)
meta_data$nGene <- as.numeric(as.character(meta_data$nGene))
meta_data <- meta_data[meta_data$nGene > 300 & meta_data$percent_mito < .2,]
data <- data[,row.names(meta_data)]

# Normalize and scale

data <- as(data, "dgCMatrix")
var_genes <- vargenes_vst(data, meta_data$Subject, 200)
exprs_norm <- normalizeData(data, method = "log")
cc_genes <- c(Seurat::cc.genes$s.genes, Seurat::cc.genes$g2m.genes)
exprs_cosine <- exprs_norm[var_genes[!var_genes %in% cc_genes],] %>% scaleData %>% cosine_normalize(2)

# PCA

pca_res_filter <- irlba::prcomp_irlba(t(exprs_cosine), 20)

# Batch correction with Harmony

harmony_res <- HarmonyMatrix(pca_res_filter$x, meta_data, c("Subject"), theta = c(1),
                               plot_convergence = TRUE, max.iter.harmony = 10, epsilon.cluster = -Inf, epsilon.harmony = -Inf,
                               max.iter.cluster = 10, do_pca = F, verbose = T, return_object = FALSE)

# UMAP for visualization

umap_res_harmony <- umap(harmony_res, n_neighbors = 30L, metric = "euclidean", min_dist = .1)
