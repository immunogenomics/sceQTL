# Script to process Blueprint data similarly to pseudobulk memory T cell data by PEER normalizing and regressing out covariates

library(data.table)
library(readr)
library(stringr)
library(peer)

K = 30

### phenotype data downloaded from here: http://dcc.blueprint-epigenome.eu/#/datasets/EGAD00001002671
pheno <- read.delim("/path/EGAD00001002663.pheno.txt", sep = "\t")
pheno <- pheno[c("DONOR_ID","DONOR_AGE","DONOR_SEX")]
pheno <- pheno[!duplicated(pheno$DONOR_ID),]
row.names(pheno) <- pheno$DONOR_ID
pheno <- pheno[c("DONOR_AGE","DONOR_SEX")]
pheno <- as.data.frame(pheno)

### Expression data (batch corrected and normalized)
norm_exp <- fread(paste0("/path/tcel_gene_nor_combat_20151109.txt.gz"))
norm_exp <- as.data.frame(norm_exp)
rownames(norm_exp) <- str_split_fixed(norm_exp$ens.id, "[.]", 2)[,1]
norm_exp <- norm_exp[,-1]
dim(norm_exp)

### select the genes that are also present in single-cell data
tbru_exp <- read_tsv("/path/tbru_IDs.txt.gz")
tbru_exp <- as.data.frame(tbru_exp)
norm_exp <- subset(norm_exp, rownames(norm_exp) %in% tbru_exp$ids)

### inverse normal transformation
rn<-apply(norm_exp,1,function(x){
    		qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) )
	})

### Genotyping PCs
pcs <- read.delim("/path/EGAZ00001235598_blueprint06092016_qc_pruned_noMHC.pc.eigenvec", sep = " ", header = F)
row.names(pcs) <- pcs$V2
pcs <- pcs[,3:ncol(pcs)]
colnames(pcs) <- c(sprintf("gPC%d", seq(1,20)))
pcs <- pcs[c("gPC1", "gPC2", "gPC3")]

### age group and sex as covariates
covs <- merge(pheno, pcs, by='row.names')
rownames(covs) <- covs[,1]
covs <- covs[,-1]
covs <- covs[rownames(covs) %in% rownames(rn),]
covs <- covs[match(rownames(rn), rownames(covs)),]
levels(covs$DONOR_AGE)  <- 1:10
covs$DONOR_SEX <- covs$DONOR_SEX == "Female"
covs$DONOR_AGE <- as.numeric(covs$DONOR_AGE)
head(covs)

### correct for covariates and 30 peer factors
model = PEER()
PEER_setPhenoMean(model, as.matrix(rn))
dim(PEER_getPhenoMean(model))
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model,K)
PEER_getNk(model)
PEER_setCovariates(model, as.matrix(covs))
PEER_setNmax_iterations(model,10000)
    
#perform the inference
PEER_update(model)

### Output the peer factors    
factors = PEER_getX(model)
dump <- data.frame(id=colnames(norm_exp),factors)
gz1 <- gzfile(paste0("/path/Blueprint_INT_NoePCs_3gPCs_age_sex_peer_factors_K30.txt.gz"),"w")
write.table(dump, gz1, sep = "\t", quote = F, row.names = FALSE)
close(gz1)

### output the residulas
residuals = PEER_getResiduals(model)
colnames(residuals)  <- row.names(norm_exp)
row.names(residuals) <- colnames(norm_exp)
dump <- data.frame( ids = row.names(norm_exp),t(residuals))
gz1 <- gzfile(paste0("/path/Blueprint_INT_NoePCs_3gPCs_age_sex_peer_factors_K30_residuals.txt.gz"),"w")
write.table(dump, gz1, sep = "\t", quote = F, row.names = FALSE)
close(gz1)
