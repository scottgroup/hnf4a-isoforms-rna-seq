log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# Loading count matrix from file, converting to integer
counts <- as.matrix(
    read.table(
        snakemake@input[["counts"]], header=TRUE,
        row.names="gene", check.names=FALSE
    )
)
mode(counts) <- "integer"
counts[,'A3N1':=NULL]

# Loading samples information from file.
# IMPORTANT -> Replicate A3N1 was removed due to being an outlier
#             (low HNF4a expression)
samples <- read.table(
    snakemake@input[["samples"]], header=TRUE,
    row.names="sample", check.names=FALSE
)

# Calculating DESeq2
dds <- DESeqDataSetFromMatrix(
    countData=counts,
    colData=samples,
    design= ~condition
)
dds <- DESeq(dds)

# Compiling results for each condition
for (i in 1:12) {
    exp <- sprintf("A%s", i)
    res <- results(dds, contrast=c("condition","A0",exp))

    # Writing results to file
    fname <- paste(
        snakemake@output[["results"]],
        paste(exp, "csv", sep='.'), sep='/'
    )

    write.csv(
        as.data.frame(res),
        file=fname,
        quote=FALSE
    )


}
