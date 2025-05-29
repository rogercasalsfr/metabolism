## Counts to TPM
remotes::install_github("IOBR/IOBR")

library(IOBR)


counts <- read.table("C:/Users/Roger/Downloads/adjusted_counts.tsv",
                      header = TRUE,
                      sep = "\t",
                      row.names = 1,
                      check.names = FALSE)

# Pass count to tpm
data_tpm <- count2tpm(countMat = counts, source = "local", idType = "Symbol")


write.csv(data_tpm, "C:/Users/Roger/Downloads/TPM_results.csv", row.names = TRUE)

write.table(data_tpm,
            file = "C:/Users/Roger/Downloads/TPM_results.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)


# Check median and mean, because we want robust models. 
# We just subset time = 4.

metadata_4 = subset(metadata, time %in% 'four')
data_tpm_final = data_tpm[colnames(data_tpm) %in% metadata_4$sampleID]


# 1. Ensure sample IDs in metadata match the column names in data_tpm_final
metadata_4$sampleID <- as.character(metadata_4$sampleID)
colnames(data_tpm_final) <- as.character(colnames(data_tpm_final))

# 2. Subset metadata to only samples present in the expression matrix
metadata_filtered <- metadata_4[metadata_4$sampleID %in% colnames(data_tpm_final), ]

# 3. Reorder columns of data_tpm_final to match the order in metadata
data_tpm_filtered <- data_tpm_final[, metadata_filtered$sampleID]

# 4. Split sample names by treatment
grouped_samples <- split(metadata_filtered$sampleID, metadata_filtered$treatment)

# 5. For each group, compute the row-wise median
median_by_treatment <- sapply(grouped_samples, function(samples) {
  apply(data_tpm_filtered[, samples, drop = FALSE], 1, median, na.rm = TRUE)
})


# 5. Compute row-wise mean per treatment
mean_by_treatment <- sapply(grouped_samples, function(samples) {
  apply(data_tpm_filtered[, samples, drop = FALSE], 1, mean, na.rm = TRUE)
})

# 6. Convert to data frames and set row names
median_df <- as.data.frame(median_by_treatment)
mean_df <- as.data.frame(mean_by_treatment)
rownames(median_df) <- rownames(data_tpm_final)
rownames(mean_df) <- rownames(data_tpm_final)

# 7. Optional: view results
head(median_df)
head(mean_df)



write.table(median_df,
            file = "C:/Users/Roger/Downloads/TPM_results_mean.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)


boxplot(log2(median_df+1))
