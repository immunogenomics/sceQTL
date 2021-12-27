# Script to calculate enrichment of eQTL signal in transcription start sites

# Load libraries
library(data.table)
library(dplyr)
library(tidyr)

# Load table of eGenes and lead SNPs

genes <- read.table("/path/pseudobulk_sig_eqtls.txt", stringsAsFactors = F)

# Load keys to convert gene names to Ensembl IDs (from Cell Ranger refdata-cellranger-GRCh38-3.0.0) and GRCh38 variants to hg19 (based on liftOver)
ensg_ids <- fread("/path/features.tsv.gz")
conversion_key <- fread("/path/convert_grch38_hg19.txt.gz")

# Load +/- 2kb windows around all 6,511 eGenes
my_tss <- read.table("/path/tss_track.txt", header = F)

genes <- genes %>% left_join(ensg_ids, by = c("V3"="V2")) %>% separate(V1, into = c("chr", "pos", "ref", "alt"), sep = "_")

results <- sapply(1:nrow(genes), function(x) {
	if(file.exists(paste0("/path/", genes$chr[x], "/", genes$chr[x], "_", genes$V1.y[x], "_SingleCausal.caviar_post"))) {
	cav <- read.table(paste0("/path/", genes$chr[x], "/", genes$chr[x], "_", genes$V2.y[x], "_SingleCausal.caviar_post"), header = T)
	if(file.exists(paste0("/path/", genes$chr[x], "/", genes$chr[x], "_", genes$V2.y[x], "_ConditionSecondSNP_SingleCausal.caviar_post"))) {
        	cav <- read.table(paste0("/path/", genes$chr[x], "/", genes$chr[x], "_", genes$V2.y[x], "_ConditionSecondSNP_SingleCausal.caviar_post"), header = T)}
	if(max(cav$Prob_in_pCausalSet) >=.5) {
		cav <- cav %>% left_join(conversion_key, by = c("SNP_ID"="V1")) %>% separate(V2, into = c("chr", "pos", "ref", "alt"), sep = "_")
		cav$overlap = Inf
		tmp <- my_tss %>% filter(V1 == genes$chr[x])
		cav$overlap = sapply(as.numeric(cav$pos), function(y){sum(tmp$V2 <= y & tmp$V3 >= y)})
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
