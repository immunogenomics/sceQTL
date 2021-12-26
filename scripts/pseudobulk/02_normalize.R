# Script to normalize pseudobulk expression

# Load libraries
library(ggplot2)
library(Seurat)
library(data.table)
library(ggrepel)
library(peer)
library(edgeR)

# Load donor-level covariates (from dbGaP: phs002467)
pheno <- read.delim("/path/donor_covariates.txt", sep = " ")
pheno <- pheno[!duplicated(pheno$donor),]
row.names(pheno) <- pheno$donor

# Load genotype PCs
pcs <- read.delim("/path/geno_PCs.txt", sep = " ", header = F)
row.names(pcs) <- pcs$V2
pcs <- pcs[,3:ncol(pcs)]

# PEER normalization: regress out latent variables (Stegle 2010, PLoS Comp Bio)

# Based on GTEx 2017 Nature:
K = 45

# Replace with reading in your pseudobulk expression matrix
sum_exp <- fread("/path/pseudobulk_exprs_mat.txt.gz")
sum_exp <- as.data.frame(sum_exp)
row.names(sum_exp) <- sum_exp[,1]
sum_exp <- sum_exp[,-1]

# Remove genes that with non-zero expression in <= half of the samples
sum_exp <- sum_exp[rowSums(sum_exp > 0) > .5*ncol(sum_exp),]

# Log2 CPM normalization
norm_exp <- log2(cpm(sum_exp)+1)

# Inverse normal transformation
rn<-apply(norm_exp,1,function(x){
        qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) )
})
rn<-t(rn)

# Select known covariates to regress out (age, sex, and 5 genotype PCs)
covs <- cbind(pheno[colnames(norm_exp),"age"],
                  pheno[colnames(norm_exp),"Sex"] == "F",
               pcs[colnames(norm_exp),1:5])
colnames(covs) <- c("age", "female", "PC1", "PC2", "PC3", "PC4", "PC5")

model = PEER()
PEER_setPhenoMean(model, as.matrix(t(rn)))
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model,K)
PEER_getNk(model)
PEER_setCovariates(model, as.matrix(covs))
PEER_setNmax_iterations(model,10000)
PEER_update(model)

# Save residuals
residuals = PEER_getResiduals(model)
colnames(residuals)  <- row.names(norm_exp)
row.names(residuals) <- colnames(norm_exp)
dump <- data.frame( ids = row.names(norm_exp),t(residuals))
gz1 <- gzfile("/path/peer_residuals.txt.gz","w")
write.table(dump, gz1, sep = "\t", quote = F, row.names = FALSE)
close(gz1)
