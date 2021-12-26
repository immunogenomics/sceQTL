# Script to make pseudobulk expression profiles

# Load libraries
library(data.table)
library(presto)
library(DESeq2)
library(edgeR)
library(dplyr)

# Load metadata (from GEO: GSE158769)
meta_data <- fread("/path/GSE158769_meta_data.txt.gz")
meta_data <- meta_data %>% mutate(id = paste(donor, batch, sep = "_"))

# Load RAW UMI COUNTS matrix (from GEO)
exprs_raw <- fread("/path/GSE158769_exprs_raw.tsv.gz")

# Remove replicates
reps <- meta_data %>% filter(!duplicated(id)) %>% arrange(batch) %>% filter(duplicated(donor)) %>% select(id) %>% unlist
meta_data <- meta_data %>% filter(!id %in% reps)
exprs_raw <- exprs_raw[,meta_data$cell_id]

# Make pseudobulk profiles; specify variable to make pseudobulk samples for
all_collapse <- collapse_counts(exprs_raw, meta_data, c("donor"))
colnames(all_collapse$counts_mat) <- all_collapse$meta_data$donor

# Save genes (rows) x samples (columns) matrix

out <- data.frame(id=row.names(all_collapse$counts_mat), all_collapse$counts_mat)

gz1 <- gzfile("/path/pseudobulk_exprs_mat.txt.gz","w")
write.table(out, gz1, sep = "\t", quote = F, row.names = FALSE)
close(gz1)
