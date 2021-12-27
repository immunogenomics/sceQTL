# Scripts to conduct trajectory analysis on monocytes from Randolph, et al. (2021) Science

# Load libraries

library(Seurat)
library(Matrix)
library(singlecellmethods)
library(dplyr)
library(harmony)
library(uwot)
library(DDRTree)
library(princurve)

# Load post-QC data

exprs_raw <- readRDS("/path/randolph_raw.rds")
exprs_raw <- exprs_raw[!row.names(exprs_raw) %in% c("PB2", "PB1", "PA", "HA", "NP", "NA", "M2", "M1", "NEP", "NS1"),]
meta.data <- readRDS("/path/randolph_meta.rds")

# Normalize and scale

exprs_raw <- exprs_raw[,row.names(meta.data)][,grepl("monocyte", meta.data$celltype)]
mono <- CreateSeuratObject(counts = exprs_raw)
mono <- PercentageFeatureSet(mono, pattern = "^MT-", col.name = "percent.mt")
mono <- SCTransform(mono, vars.to.regress = "percent.mt", verbose = FALSE)

# PCA

pca_res_filter <- irlba::prcomp_irlba(t(mono@assays$SCT@scale.data), 20)

# Batch correction with Harmony

harmony_res <- HarmonyMatrix(pca_res_filter$x, meta.data[grepl("monocyte", meta.data$celltype),], c("SOC_indiv_ID", "batchID", "SOC_infection_status"), theta = c(1,1,1),
                               plot_convergence = TRUE, max.iter.harmony = 10, epsilon.cluster = -Inf, epsilon.harmony = -Inf,
                               max.iter.cluster = 10, do_pca = F, verbose = T, return_object = FALSE)

# UMAP for visualization

umap_res_harmony <- umap(harmony_res[,1:10], n_neighbors = 30L, metric = "euclidean", min_dist = .1)

# Plot cells by treatment condition

fig.size(4,5)
ggplot(as.data.frame(umap_res_harmony), aes(V1, V2, color = meta.data[grepl("monocyte", meta.data$celltype), "SOC_infection_status"])) + 
    geom_point(shape = ".") + theme_classic() + theme(legend.title = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 12)) +
    xlab("UMAP1") + ylab("UMAP2")

# Trajectory analysis with DDRTree

ncells <- nrow(harmony_res)
ncenter <- round(2 * 100 * log(ncells) / (log(ncells) + log(100)))

ddr_args <- c(list(
        X = t(harmony_res), 
        dimensions = 2, ## LOW DIMENSIONALITY
        ncenter = ncenter, ## number of nodes allowed in the regularization graph
        param.gamma = 10, ## param.gamma regularization parameter for k-means 
        maxIter = 20,
        tol = 1e-3,
        sigma = 0.0001,
        verbose = FALSE))

ddrtree_res <<- do.call(DDRTree, ddr_args)
pseud_integr <- t(ddrtree_res$Z)
pseud_integr <- as.data.frame(pseud_integr)

pc_res <- principal_curve(t(ddrtree_res$Z))
pseud_integr$pseudo <- pc_res$lambda

# Plot trajectory across cells

fig.size(4,5)
ggplot(data = as.data.frame(umap_res_harmony), aes(x = V1, y = V2, color = pseud_integr$pseudo)) + 
    geom_point(shape = ".") + 
    theme_classic() +
    theme(legend.text=element_text(size=12)) + scale_color_viridis(option = "plasma") +
    theme(legend.title=element_blank(), axis.text = element_blank()) +
    xlab("UMAP1") + ylab("UMAP2")

# Plot deciles by treatment condition composition
data.frame(V3 = meta.data.mono$SOC_infection_status) %>% mutate(bin = cut(pseud_integr$pseudo, breaks = quantile(pseud_integr$pseudo, probs = seq(0, 1, .1), labels = 1:10, include.lowest = TRUE))) %>%                                 
#    group_by(bin, V3) %>% summarise(count = n()) %>% filter(!is.na(bin)) %>%
    filter(!is.na(bin)) %>% ggplot(aes(bin, fill = V3)) + geom_bar(stat = "count", position = "fill") + theme_classic() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
