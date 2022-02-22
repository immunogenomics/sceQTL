# Script to prepare files for fine-mapping from FastQTL nominal output

gene = "CTLA4"

# Read fastqtl nominal output
data <- read_delim(path1, delim=" ", col_names=F)
colnames(data) <- c("EnsemblID", "var", "dist", "P", "beta")
data <- data %>% filter(EnsemblID == gene)

# Read minor allele frequencies per variant
maf <- read.table(path2, header=T)

data <- merge(data, maf, by.x="var", by.y="SNP")

data$Z <- qnorm(data$P/2, lower.tail = FALSE)
data$Z <- ifelse(data$beta < 0, -1*data$Z, data$Z)
data$SE <- data$beta/data$Z

data <- data[,c("var", "Z")]
data <- unique(data)

write.table(data, paste(path3,"_", gene, ".TMP.z", sep=""), quote=F, col.names=F, row.names=F)
