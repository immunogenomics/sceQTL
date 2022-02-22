# Calculate in-sample LD between variants

chr=$1
gene=$2


/path/plinkv1.90b6.22 --vcf merged_${chr}.QC.vcf.gz \
--extract merged_${chr}_${gene}.TMP.snpIDs \
--keep-allele-order \
--r square \
--out merged_${chr}_${gene}.TMP

/path/plinkv1.90b6.22 --vcf merged_${chr}.QC.vcf.gz \
--extract merged_${chr}_${gene}.TMP.snpIDs \
--keep-allele-order \
--freq \
--out merged_${chr}_${gene}.TMP

sed -i 's/\t/ /g'  merged_${chr}_${gene}.TMP.ld
sed -i 's/nan/0/g'  merged_${chr}_${gene}.TMP.ld
