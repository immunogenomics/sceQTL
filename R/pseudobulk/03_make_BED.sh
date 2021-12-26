#!/bin/bash

# Script to make BED file from PEER residuals

# BED header
echo -e "#Chr\tStart\tEnd\tID" > tmp
zcat /path/peer_residuals.txt.gz | head -n1 | cut -f2- | sed 's/\./-/g' | paste tmp - > header

zcat /path/peer_residuals.txt.gz | sed -e "1d" |
sort -k 1b,1 > tmp01

# Needs to be by chromosome for FastQTL
for chr in $(seq 1 22);do
      echo chr$chr ...

      # Get gene TSS using GTFs from Cell Ranger (refdata-cellranger-GRCh38-3.0.0)
      sed -e "1d" /path/chr${chr}_tss_gene_basic.bed |
        awk -v chr=$chr '{if($1==chr){print $2,chr,$4 -1,$4}}' |
      awk '{ $1 = substr($1, 1, 15)} 1' |
      sort -k 1b,1 > tmp02

      outfile="/path/chr${chr}_input.txt"
      # Link ENSG to gene name using features.tsv.gz file from cellranger-3.1.0, GRCh38
      zcat /path/features.tsv.gz |
      grep ^ENSG |
      sort -k 1b,1 |
      join - tmp02 |
      awk '{print $2,$5,$6,$7,$2}' |
      sort -k 1b,1 |
      uniq -f 4 |
      join  - tmp01 |
      cut -d " " -f2- |
      sort -k 2,2n |
      perl -pe "s/ /\t/g" |
      awk '{print $0}' |
      cat header - > "$outfile"

      bgzip -f "$outfile"

      tabix -f -p bed "$outfile.gz"
done
