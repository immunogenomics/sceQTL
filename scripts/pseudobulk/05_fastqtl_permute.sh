#!/bin/sh

# Script to run FastQTL (Ongen 2015, Bioinformatics)

# VCF: post QC and imputation genotype dosage VCF (genotype data from dbGaP: phs002025)
# BED: output from previous script
# Region: select chromosome
# out: output file

chr=$1

mkdir -p /path/fastqtl_permute/chr$chr

/path/fastQTL \
   --vcf  /path/geno_imputed_dosage_maf05_hg38_QC.vcf.gz \
   --bed  /path/chr${chr}_input.bed.gz   \
   --region  chr$chr \
   --permute 1000 \
   --out  /path/fastqtl_nominal/chr$chr/chr${chr}_output.permute.tbru.txt

gzip -f /path/fastqtl_nominal/chr$chr/chr${chr}_output.permute.tbru.txt
