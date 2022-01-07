# Script to calculate enrichment of eQTL signal in ATAC-seq peaks from sorted T cell types in Calderon, et al. (2019)Nat Genetics

# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(edgeR)

args = commandArgs(trailingOnly=TRUE)
celltype = args[1]
n = args[2]

# Load table of eGenes and lead SNPs

genes <- read.table("/path/pseudobulk_sig_eqtls.txt", stringsAsFactors = F)

# Load keys to convert gene names to Ensembl IDs (from Cell Ranger refdata-cellranger-GRCh38-3.0.0) and GRCh38 variants to hg19 (based on liftOver)
ensg_ids <- fread("/path/features.tsv.gz")
conversion_key <- fread("/path/convert_grch38_hg19.txt.gz")

# Define open regions based on ATAC
atac <- fread("/data/srlab2/anathan/external_data/Calderon2019/GSE118189_ATAC_counts.txt.gz")
atac <- data.frame(atac, check.names = F)
atac <- atac[,grepl(paste(c("V1",unique(substring(colnames(atac), 6, nchar(colnames(atac))-2))[6:19]), collapse = "|"), colnames(atac))]

atac[,-1] <- cpm(atac[,-1])
atac_bed <- data.frame(V1 = atac$V1) %>% separate(V1, into = c("V1", "V2", "V3"), sep = "_")
atac_bed$V2 <- as.numeric(atac_bed$V2)
atac_bed$V3 <- as.numeric(atac_bed$V3)
atac_bed$V4 = rowMeans(atac[,grepl(paste0("-", celltype, "-U"), colnames(atac))])
atac_bed$size = as.numeric(atac_bed$V3) - as.numeric(atac_bed$V2)

genes <- genes %>% left_join(ensg_ids, by = c("V3"="V2")) %>% separate(V1, into = c("chr", "pos", "ref", "alt"), sep = "_")

results <- sapply(1:nrow(genes), function(x) {
	if(file.exists(paste0("/path/", genes$chr[x], "/", genes$chr[x], "_", genes$V1.y[x], "_SingleCausal.caviar_post"))) {
	cav <- read.table(paste0("/path/", genes$chr[x], "/", genes$chr[x], "_", genes$V2.y[x], "_SingleCausal.caviar_post"), header = T)
	if(file.exists(paste0("/path/", genes$chr[x], "/", genes$chr[x], "_", genes$V2.y[x], "_ConditionSecondSNP_SingleCausal.caviar_post"))) {
        	cav <- read.table(paste0("/path/", genes$chr[x], "/", genes$chr[x], "_", genes$V2.y[x], "_ConditionSecondSNP_SingleCausal.caviar_post"), header = T)}
	if(max(cav$Prob_in_pCausalSet) >=.5) {
		cav <- cav %>% left_join(conversion_key, by = c("SNP_ID"="V1")) %>% separate(V2, into = c("chr", "pos", "ref", "alt"), sep = "_")
		cav$overlap = Infi
		tmp <- atac_bed %>% filter(V1 == genes$chr[y])
                cav$overlap = sapply(as.numeric(cav$pos), function(x){sum(tmp$V2 <= x & tmp$V3 >= x & tmp$V4 > n)})
                e_overlap = sum(cav$overlap)
                e_pip = sum(cav$Prob_in_pCausalSet)
                e_num = sum(cav$Prob_in_pCausalSet*cav$overlap)
                return(c(x, e_overlap, e_pip, e_num, nrow(cav), sum(cumsum(sort(cav$Prob_in_pCausalSet, decreasing = T)) <= .95)))}
        else{
                return(c(NA,NA,NA,NA,NA,NA))
        }}

        else {
                return(c(NA,NA,NA,NA,NA,NA))	
}})

# Load eQTL output (data frame from concatenating output of single-cell eQTL models; 195,330 rows = 6,511 eGenes x 30 fixed effects estimated for each)
sceqtl_results <- readRDS("/path/sceqtl_results.rds")

genes_sig <- sceqtl_results %>% filter(term == "G") %>% mutate(qval = qvalue::qvalue(lrt_pval)$qvalue) %>% filter(qval < .05) %>% select(gene) %>% unlist

results <- results %>% filter(!is.na(X1))

# Find intersection with eGenes significant in Peruvian/Blueprint joint meta-analysis
genes_bp <- read.table("/path/joint_bp_sig_eqtls.txt")
genes_int <- genes[genes$V3 %in% genes_bp$V3,]

# Calculate mean enrichment stats
results %>% mutate(es = X4/X5/(X2/X5*X3/X5)) %>% filter(!is.na(es)) %>% mutate(gene = genes$V3[X1]) %>% filter(gene %in% genes_int$V3) %>% select(es) %>% unlist %>% mean
results %>% mutate(es = X4/X5/(X2/X5*X3/X5)) %>% filter(!is.na(es)) %>% mutate(gene = genes$V3[X1]) %>% filter(gene %in% genes_int$V3) %>% filter(gene %in% genes_sig) %>% select(es) %>% unlist %>% mean
results %>% mutate(es = X4/X5/(X2/X5*X3/X5)) %>% filter(!is.na(es)) %>% mutate(gene = genes$V3[X1]) %>% filter(gene %in% genes_int$V3) %>% filter(!gene %in% genes_sig) %>% select(es) %>% unlist %>% mean
