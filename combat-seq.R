
# Install package Combat-seq from source
# Download the tar from sva Bioconductor. https://www.bioconductor.org/packages/release/bioc/html/sva.html
install.packages("C:/Users/Roger/Downloads/sva_3.56.0.tar.gz", repos = NULL, type = "source")


library(sva)

# Load count matrix
counts = read.csv("C:/Users/Roger/Downloads/kallisto.100/kallisto.100/DESeq2_TPM_values.tsv", sep="\t", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))  # remove X from columns

# Remove genes that make all the analysis fail
genes_to_remove <- c("MTND2P30", "RPL23AP95", "USP9YP10")
counts <- counts[!rownames(counts) %in% genes_to_remove, ]



# Load metadata, and add the vector from batch correction
metadata = read.csv("C:/Users/Roger/Downloads/metadata.tsv", sep="\t")
metadata <- subset(metadata, sampleID %in% colnames(counts))   # filter metadata

# re-order counts
counts <- counts[, metadata$sampleID]

# We re-assure that the names match in order
all(colnames(counts) == metadata$sampleID)

# We define our batch vector (experiment)
batch <- metadata$experiment

# We define our biological variables to not remove effect. 
metadata$treatment <- as.factor(metadata$treatment)
metadata$time <- as.factor(metadata$time)

cov1 <- metadata$time
cov2 <- metadata$treatment
covar_mat <- cbind(cov1, cov2)

# Adjust values 
adjusted_counts <- ComBat_seq(as.matrix(counts), batch=batch, covar_mod = covar_mat)


# Convert into .tsv
df <- as.data.frame(adjusted_counts)


write.table(df,
            file = "C:/Users/Roger/Downloads/adjusted_counts.tsv",
            quote = FALSE,
            sep = "\t",
            col.names = NA,
            row.names = TRUE)




