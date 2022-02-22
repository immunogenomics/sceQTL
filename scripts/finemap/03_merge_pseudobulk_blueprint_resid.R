library(data.table)
library(stringr)
library(readr)

## get GTF files of TSS for each gene in Cell Ranger reference (refdata-cellranger-GRCh38-3.0.0)
files <- Sys.glob(file.path("/path/chr*_tss_gene_basic.bed"))

tss <- rbindlist(lapply(files, function(filename) {
  read.table(filename, header=T)
}), fill=TRUE)

tss$CHR <- str_remove(tss$CHR, "chr")
tss$ENS <- str_split_fixed(tss$ID, "[.]", 2)[,1]
tss$start <- tss$POS-1

tss <- tss[,c(1,6,4,5)]
colnames(tss) <- c("CHR", "Start", "End", "ID")

## combine residuals from the two datasets -> 13978 genes
tbru_exp <- read_tsv("/path/peer_residuals.txt.gz")
tbru_exp <- as.data.frame(tbru_exp)
rownames(tbru_exp) <- tbru_exp$ids
tbru_exp <- tbru_exp[,-1]
blue_exp <- read_tsv("/path/Blueprint_INT_NoePCs_3gPCs_age_sex_peer_factors_K30_residuals.txt.gz")
blue_exp <- as.data.frame(blue_exp)
rownames(blue_exp) <- blue_exp$ids
blue_exp <- blue_exp[,-1]

merged_exp <- merge(tbru_exp, blue_exp, by='row.names')
merged_exp <- as.data.frame(merged_exp)
dim(merged_exp)
merged_exp[1:4,1:5]

### Merge with TSS data from Cell Ranger GTF -> 13963 genes
data <- merge(tss, merged_exp, by.x="ID", by.y="Row.names")
data <- data[with(data, order(as.numeric(CHR), as.numeric(End))),]
data <- as.data.frame(data)
data <- data[,c(2,3,4,1,5:ncol(data))]

### fix names to match vcf
colnames(data) <- str_replace_all(colnames(data), "[.]", "-")

### create a cohort covariate
#fam <- read.table("/path/SampleIDs.txt", header=F)
#fam <- fam %>% dplyr::mutate(cov=ifelse(grepl('-',V1), 1, 0))

## order expression like fam
samp_names <- c(colnames(data[,1:4]), as.character(fam$V1))
data <- data[,match(samp_names, colnames(data))]
data[,1] <- paste("chr", data[,1], sep="")
colnames(data)[1] <- "#CHR"

reg_exp <- sapply(1:nrow(data), function(x){resid(lm(unlist(data[x,-c(1:4)])~ fam$cov))})
res <- data.frame(data[,1:4],t(reg_exp))
colnames(res)[1] <- "#CHR"
colnames(res) <- str_replace_all(colnames(res), "[.]", "-")
write.table(res,"/path/merged.bed", sep = "\t", quote = F, row.names = FALSE)

##bgzip and index
