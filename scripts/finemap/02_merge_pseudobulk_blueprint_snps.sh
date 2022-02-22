#!/bin/bash

module load vcftools/0.1.15
module load tabix/0.2.6
module load bcftools/1.10.2

zcat /path/geno_imputed_dosage_maf05_hg38_${chr}_QC.vcf.gz | grep -v "#" | cut -f3,3 | sort | uniq -c | awk '{if (\$1 == 1) print \$2}'  > /path/merged_chr${chr}.IDs.txt
