log <- file(snakemake@log[[1]], open="wt")
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(snakemake@input[["samples"]], header=TRUE, row.names="rep_id", check.names=FALSE)
coldata$group <- as.factor(coldata$group)

dds <- DESeqDataSetFromMatrix(countData=round(cts),
                              colData=coldata,
                              design=~ group)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])
