#!/bin/bash

# Script to merge genotypes from pseudobulk memory T cells and Blueprint

module load vcftools/0.1.15
module load tabix/0.2.6
module load bcftools/1.10.2

# Format pseudobulk genotype VCF
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' /path/geno_imputed_dosage_maf05_hg38_${chr}.vcf.gz | bgzip -c > /path/geno_imputed_dosage_maf05_hg38_${chr}_QC.vcf.gz
tabix -f -p /path/geno_imputed_dosage_maf05_hg38_${chr}_QC.vcf.gz

# Format Blueprint genotype VCF
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' /path/EGAZ00001235598_blueprint06092016_qc_hg38_PICARD_chr${chr}.vcf.gz | bgzip -c > /path/EGAZ00001235598_blueprint06092016_qc_hg38_PICARD_chr${chr}.QC.vcf.gz
tabix -f -p vcf BLUEPRINT/PlinkFiles/EGAZ00001235598_blueprint06092016_qc_hg38_PICARD_chr${chr}.QC.vcf.gz

# Merge VCFs
/path/vcftools/0.1.15/bin/vcf-merge /path/EGAZ00001235598_blueprint06092016_qc_hg38_PICARD_chr${chr}.QC.vcf.gz /path/geno_imputed_dosage_maf05_hg38_${chr}_QC.vcf.gz | bgzip -c > /path/merged_chr${chr}.vcf.gz
tabix -f -p vcf /path/merged_chr${chr}.vcf.gz

# Get pseudobulk variants
zcat /path/geno_imputed_dosage_maf05_hg38_${chr}_QC.vcf.gz | grep -v "#" | cut -f3,3 | sort | uniq -c | awk '{if (\$1 == 1) print \$2}'  > /path/pseudobulk_chr${chr}.VarIDs.txt

# Blueprint variants
zcat /path/EGAZ00001235598_blueprint06092016_qc_hg38_PICARD_chr${chr}.QC.vcf.gz | grep -v "#" | cut -f3,3 > /path/Blueprint_chr${chr}.VarIDs.txt

# Get shared variant IDs
cat /path/pseudobulk_chr${chr}.VarIDs.txt /path/Blueprint_chr${chr}.VarIDs.txt | sort | uniq -c | awk '{if (\$1 == 2) print \$2}' > /merged_chr${chr}.VarIDs.txt

# Merged VCF of shared variants
/path/vcftools/0.1.15/bin/vcftools --gzvcf /path/merged_chr${chr}.vcf.gz --snps /path/merged_chr${chr}.VarIDs.txt --recode --out /path/merged_chr${chr}.QC --keep /path/merged_QC.vcftools.keep
bgzip -f  /path/merged_chr${chr}.QC.vcf
tabix -f -p vcf  /path/merged_chr${chr}.QC.vcf.gz
