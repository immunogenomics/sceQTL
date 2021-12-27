# Script to calculate single-cell betas from single-cell eQTL model interaction effects and CV scores

# Load libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Load eQTL results, CV scores, and UMAP coordinates

sceqtl_results <- readRDS("/path/sceqtl_results.rds")
cv_scores <- readRDS("/path/cca_donor_batch.rds")[,1:7]
umap_res <- readRDS("/path/umap_donor_batch.rds")

# Calculate single-cell betas

betas <- sceqtl_results %>% filter(grepl("^G", term) & gene == "GZMA") %>% select(Estimate) %>% unlist
sc_betas <- t(data.frame(X1 = rowSums(sweep(cbind(1, cv_scores[,1:7]), MARGIN=2, betas, `*`))))

# Plot eQTL betas per cell

plot_eqtl <- function (ab, umap, exprs, pct = 0.95, geno)
{
    max.cutoff = quantile(exprs[ab, ], pct)
    min.cutoff = quantile(exprs[ab, ], 1 - pct)
    tmp <- sapply(X = exprs[ab, ], FUN = function(x) {
        return(ifelse(test = x > max.cutoff, yes = max.cutoff,
            no = x))
    })
    tmp <- sapply(X = tmp, FUN = function(x) {
        return(ifelse(test = x < min.cutoff, yes = min.cutoff,
            no = x))
    })
    umap_res_plot <- cbind(umap, tmp)
    return(ggplot(data = as.data.frame(umap_res_plot)[sample(nrow(umap_res_plot)),
        ], aes(x = V1, y = V2)) + geom_point_rast(mapping = aes(color = tmp),
        shape = ".") + scale_color_distiller(palette = "RdYlBu",
        limits = c(geno - max(abs(geno - tmp)), geno + max(abs(geno -
            tmp)))) + theme_classic() + theme(axis.text = element_blank(),
        axis.title = element_blank()))
}

plot_eqtl("X1", umap_res, sc_betas, 1, betas[1])
