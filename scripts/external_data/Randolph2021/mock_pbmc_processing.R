# Scripts to process mock-treated PBMC single-cell data from Randolph, et al. (2021) Science

# Load libraries

library(Seurat)
library(Matrix)
library(singlecellmethods)
library(dplyr)
library(harmony)
library(uwot)

# Load data
ge_counts <- readRDS('/path/randolph_raw.rds') # 19248 236993 (genes x cells)
meta.data <- readRDS('/path/randolph_meta.rds') #  236993     17

# Remove IAV genes from Gene count matrix
iav_genes <- c('PB2','PB1','PA','HA','NP','NA','M2','M1','NEP','NS1')
ge_counts <- ge_counts[!rownames(ge_counts) %in% iav_genes,]

# Select MOCK
meta.data <- meta.data[meta.data$SOC_infection_status=='NI',] # 124976 cells 
ge_counts <- ge_counts[,colnames(ge_counts) %in% rownames(meta.data)]

# Select cells with > 500 genes
meta.data <- meta.data[meta.data$nFeature_RNA>500,] # 108,209    18
ge_counts <- ge_counts[ ,colnames(ge_counts) %in% rownames(meta.data)]

# Normalize and scale variable genes
mrnaNorm <- singlecellmethods::normalizeData(ge_counts, method = "log")
var_genes <- singlecellmethods::vargenes_vst(object = ge_counts, topn = 3000)
mrnaCosScaled <- mrnaNorm[rownames(mrnaNorm) %in% var_genes, ] %>% scaleData() %>% cosine_normalize(2)

# PCA
pca_mrnaScaled <- irlba::prcomp_irlba(t(mrnaCosScaled ), 20)

# Harmony batch correction
harmony <- HarmonyMatrix(pca_mrnaScaled$x, meta.data, c("SOC_indiv_ID","batchID"), theta = c(1,1), lambda = c(1,1),
                               plot_convergence = TRUE, nclust = 100, max.iter.harmony = 20,
                               max.iter.cluster = 20, do_pca = F, verbose = T)

# UMAP
umap_mrnaScaled <- uwot::umap(harmony, n_neighbors = 30, metric = "euclidean", min_dist = .1)
