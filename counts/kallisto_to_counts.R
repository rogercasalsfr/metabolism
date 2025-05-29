library(DESeq2)

library(biomaRt)
library(tximport)



#
# 1. user-defined variables
#
setwd("C:/Users/Roger/Downloads/")
kallisto_dir = "C:/Users/Roger/Downloads/kallisto.100/kallisto.100"
results_dir = kallisto_dir
# options = ' -i {} -o kallisto_output_a -t {} -b 100 --rf-stranded --verbose '.format(transcriptome_index, number_threads)

#
# 1. generate gene to transcript mapping
#
df = read.csv("C:/Users/Roger/Downloads/homo_sapiens/homo_sapiens/transcripts_to_genes.txt", sep='\t')
t2g = df[, c(1, 3)]
dim(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('072', dirnames)]  # Keep only those with "072" in the name
paths = file.path(dirnames, 'abundance.h5')

labels = basename(dirnames)  # Extract "072_159", etc.

metadata = data.frame(labels)
metadata$path = paths
View(metadata)

#
# 3. read files
#
library(rhdf5)
# Remove transcript version
t2g[,1] <- sub("\\..*", "", t2g[,1])

txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

#
# 4. find abundance


tpm = round(txi$counts)
colnames(tpm) = metadata$labels
dim(tpm)
View(tpm)

#
# 5. store
#
store = paste(results_dir, '/DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)
