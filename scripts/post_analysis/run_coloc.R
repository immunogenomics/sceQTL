# Script to run coloc on GWAS summary statistics and pseudobulk eQTL results

# Load libraries

library(data.table)
library(coloc)
library(dplyr)
library(tidyr)

# Load table of eGenes and lead SNPs

genes <- read.table("/path/pseudobulk_sig_eqtls.txt", stringsAsFactors = F)
genes <- genes %>% mutate(snp = V1) %>% separate(V1, into = c("chr", "pos", "ref", "alt"), sep = "_") %>% mutate(chr = gsub("chr", "", chr))
genes$nsnps = NA
genes$PP.H0.abf = NA
genes$PP.H1.abf = NA
genes$PP.H2.abf = NA
genes$PP.H3.abf = NA
genes$PP.H4.abf = NA

# Load GWAS summary statistics

gwas_sumstat <- fread("/path/24076602-GCST005531-EFO_0003885.h.tsv.gz")
gwas_sumstat <- gwas_sumstat %>% filter(!is.na(hm_variant_id)) %>% filter(!duplicated(hm_variant_id)) %>% filter(!is.na(hm_beta))

geno <- fread("/path/geno_imputed_dosage_maf05_hg38_QC.vcf.gz")
meta <- fread("/path/GSE158769_meta_data.txt.gz")

gwas_dataset = list(type = "cc",
                      snp = paste0("chr", gwas_sumstat$hm_variant_id),
                      beta = gwas_sumstat$hm_beta,
                      varbeta = gwas_sumstat$standard_error^2,
                      position = gwas_sumstat$hm_pos)
check_dataset(gwas_dataset)

for(x in 1:22){
	eqtl <- fread(paste0("/path/fastqtl_nominal/chr", x, "/chr", x, "_output.nominal.tbru.txt.gz"))

	resid <- data.frame(fread(paste0("/path/chr", x, "_input.bed.gz")))

	for(i in which(genes$chr == x)){
		G <- t(subset(geno, ID == sub(":", "_", genes$snp[i]))[,unique(as.character(meta$donor))])[,1]
		select_eqtl <- eqtl %>% filter(V1 == genes$V2[i]) %>% mutate(z = qnorm(V4/2, lower.tail = F)*sign(V5), se = V5/z)

		eQTL_dataset = list(beta = select_eqtl$V5,
                	            varbeta = select_eqtl$se^2,
                        	    sdY = sd(unlist(resid[resid$ID == genes$V2[i],-c(1:4)])),
                          	  type = "quant",
         	                   snp = select_eqtl$V2)
		tryCatch({
			my.res <- coloc.abf(dataset1=eQTL_dataset,
                    			    dataset2=gwas_dataset)

			genes[i,8:13] = my.res$summary
		}, error=function(cond){return(NA)})
	}
}
