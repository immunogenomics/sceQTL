# Script to run single-cell Poisson mixed effects eQTL model with one cell state interaction (e.g., CV, PC, CD4+, trajectory, disease status) and simulated null (based on Buzkova, et al. 2011)

# Load libraries
library(lme4)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)

gene <- args[1]
snp <- args[2]

# Load raw UMI counts matrix and metadata (from GEO: GSE158769) and gene expression PCs
exprs_raw <- fread("/path/GSE158769_exprs_raw.tsv.gz")
meta <- fread("/path/GSE158769_meta_data.txt.gz")
pca_res <- readRDS("/path/pca_irlba.rds")
pca_res <- pca_res$x

# Load genotype dosages (from dbGaP: phs002025) and PCs
geno <- fread("/path/geno_imputed_dosage_maf05_hg38_QC.vcf.gz")
pcs <- read.delim("/path/geno_PCs.txt", sep = " ", header = F)
row.names(pcs) <- pcs$V2
pcs <- pcs[,3:ncol(pcs)]

# Make data frame of variables for model
E <- as.numeric(exprs_raw[gene,]) 
G <- t(subset(geno, ID == snp)[,as.character(meta$donor)])[,1]
IND <- factor(meta$donor)
B <- meta$batch
AGE <- scale(meta$tbru_age)
SEX <- meta$Sex
nUMI <- scale(log(meta$nUMI))
MT <- meta$percent_mito
PC <- pcs[as.character(meta$donor),1:5]
expPC <- pca_res[,1:5]

# Load interaction covariate (e.g., batch-corrected CV from CCA)
cov <- readRDS("/path/cca_donor_batch.rds")[,1]

data <- data.frame(E,G,IND,B,AGE,SEX,nUMI,MT,
                   PC1 = PC[,1],PC2 = PC[,2],PC3 = PC[,3],PC4 = PC[,4],PC5 = PC[,5],
                   expPC1 = expPC[,1], expPC2 = expPC[,2], expPC3 = expPC[,3], expPC4 = expPC[,4], expPC5 = expPC[,5], cov)
data$G <- as.numeric(as.character(data$G))

null_model <- lme4::glmer(formula = E ~ G + (1 | B) + (1 | IND) + AGE + SEX + nUMI + MT + PC1 + PC2 + PC3 + PC4 + PC5 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + cov, 
                          family = "poisson", nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap"))
null_pred_y <- predict(null_model, data,type="response")

set.seed(1010)
null_betas <- sapply(1:1000, function(x) {
        y_vec <- sapply(1:500089, function(x){rpois(1, null_pred_y[x])})
        data$y <- y_vec

        full_model <- lme4::glmer(formula = y ~ G + (1 | IND) + (1 | B) + AGE + SEX + nUMI + MT + PC1 + PC2 + PC3 + PC4 + PC5 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + cov + G*cov, 
                                  family = "poisson", nAGQ=0, data= data, control = glmerControl(optimizer = "nloptwrap"))

        return(fixef(full_model)[18])
})
    
full_model <- lme4::glmer(formula = E ~ G + (1 | B) + (1 | IND) + AGE + SEX + nUMI + MT + PC1 + PC2 + PC3 + PC4 + PC5 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + cov + G*cov, 
                          family = "poisson", nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap"))

sum(abs(fixef(full_model)[18]) > abs(null_betas))
