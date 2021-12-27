# Script to cluster eGenes and test gene set enrichments

# Load libraries
library(dplyr)
library(tidyr)
library(Seurat)
library(msigdbr)
library(RColorBrewer)
library(pheatmap)
library(parallel)
library(patchwork)

# Load eQTL output (data frame from concatenating output of single-cell eQTL models; 195,330 rows = 6,511 eGenes x 30 fixed effects estimated for each)
sceqtl_results <- readRDS("/path/sceqtl_results.rds")

# Extract genotype betas, and interaction betas and p values for 7 eQTL interaction terms per eGene

g_mat <- cvall %>% filter(term == "G") %>% select(gene, Estimate)
row.names(g_mat) <- g_mat$gene

beta_mat <- sceqtl_results %>% filter(grepl("G:CV", term)) %>% select(gene, term, Estimate) %>% spread(term, Estimate)
row.names(beta_mat) <- beta_mat$gene
beta_mat <- beta_mat[,-1]

p_mat <- cvall %>% filter(grepl("G:CV", term)) %>% select(gene, term, pval) %>% spread(term, pval)
row.names(p_mat) <- p_mat$gene
p_mat <- p_mat[,-1]

# Identify significant interactions

sig_mat <- beta_mat
sig_mat[p_mat >= .05/7/nrow(p_mat)] = 0

# Adjust beta signs

beta_mat <- beta_mat*sign(g_mat[row.names(beta_mat),"Estimate"])

# Louvain clustering of re-scaled betas

snn_ref <- BuildSNNSeurat(apply(beta_mat[rowSums(sig_mat != 0) > 1,], 1, function(x){
    lower = ifelse(any(x < 0), min(x), 0)
    upper = ifelse(any(x > 0), max(x), 0)
    extreme = ifelse(abs(lower) > upper, abs(lower), upper)

    return(x/extreme)
}) %>% t, nn.eps = 0)
resolution_list <- c(2,3,4,5)
ids_ref <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
    Seurat:::RunModularityClustering(SNN = snn_ref, modularity = 1,
        resolution = res_use, algorithm = 1, n.start = 20,
        n.iter = 20, random.seed = 100, print.output = FALSE,
        temp.file.location = NULL, edge.file.name = NULL)
}, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list))))
ids_ref <- data.frame(ids_ref)
rm(snn_ref)
gc()

# Plot heatmap of mean eGene interaction effects in clusters

pheatmap(sapply(0:7, function(x){colMeans((apply(beta_mat[rowSums(sin_mat != 0) > 1,], 1, function(x){
    lower = ifelse(any(x < 0), min(x), 0)
    upper = ifelse(any(x > 0), max(x), 0)
    extreme = ifelse(abs(lower) > upper, abs(lower), upper)

    return(x/extreme)
}) %>% t)[ids_ref$V2 == x,])}),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(200),
         breaks = seq(-1, 1, .01),
         cluster_cols = F, cluster_rows = F, scale = "none", filename = "plot.pdf", height = 6, width = 12)

# Plot eGene clusters

fig.size(6,12)
Reduce(`+`, lapply(0:7, function(x){
    ggplot(apply(beta_mat[rowSums(sig_mat != 0) > 1,], 1, function(x){
    lower = ifelse(any(x < 0), min(x), 0)
    upper = ifelse(any(x > 0), max(x), 0)
    extreme = ifelse(abs(lower) > upper, abs(lower), upper)

    return(x/extreme)
}) %>% t %>% data.frame %>% filter(ids_ref$V2 == x) %>% rownames_to_column("gene") %>% gather(CV, direction, -gene), aes(x = CV, y = direction, group = gene)) + geom_line(color = "lightgrey") + geom_line(data = data.frame(V1 = paste0("GxCV", 1:7), V2 = colMeans((apply(beta_mat[rowSums(sig_mat != 0) > 1,], 1, function(x){
    lower = ifelse(any(x < 0), min(x), 0)
    upper = ifelse(any(x > 0), max(x), 0)
    extreme = ifelse(abs(lower) > upper, abs(lower), upper)

    return(x/extreme)
}) %>% t)[ids_ref$V2 == x,])), aes(V1, V2, group = NA), color = "red") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1
))
})) + plot_layout(ncol = 4)

# Test enrichment of GO term gene sets from MSigDB

gene_sets = msigdbr(species = "Homo sapiens", category = "C5")
msigdbr_list = split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
enrich_p = lapply(msigdbr_list, function(gene_set){
	x = length(intersect(gene_set, row.names(beta_mat[rowSums(sig_mat != 0) > 1,])[ids_ref$V2 == 0]))
	y = length(intersect(gene_set, row.names(beta_mat[rowSums(sig_mat != 0) > 1,])[ids_ref$V2 != 0]))

	return(fisher.test(matrix(c(x, y, sum(ids_ref$V2 == 0)-x, sum(ids_ref$V2 != 0)-y), nrow = 2))$p.value)
})
enrich_OR = lapply(msigdbr_list, function(gene_set){
        x = length(intersect(gene_set, row.names(beta_mat[rowSums(sig_mat != 0) > 1,])[ids_ref$V2 == 0]))
        y = length(intersect(gene_set, row.names(beta_mat[rowSums(sig_mat != 0) > 1,])[ids_ref$V2 != 0]))

        return(fisher.test(matrix(c(x, y, sum(ids_ref$V2 == 0)-x, sum(ids_ref$V2 != 0)-y), nrow = 2))$estimate)
})
