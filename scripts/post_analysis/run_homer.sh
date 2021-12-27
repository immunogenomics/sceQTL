#!/bin/bash

# Script to run HOMER on promoters of eGenes interacting with given CV
x=$1

perl /path/homer/bin/findMotifs.pl /path/homer_out/CV${x}_genes.txt human /path/homer_out/CV${x}/ -start -2000 -end 2000 -len 8,10 -p 4

# Script to run HOMER on BED file of +/- 20bp windows around eQTL variants interacting with given CV
x=$1

perl /path/homer/bin/findMotifsGenome.pl /path/homer_out/CV${x}_snps.bed hg19 /path/homer_out/CV${x}_snps_given/ -size given
