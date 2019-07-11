log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

#Loading count matrix from file, converting to integer
counts <- read.table(
	snakemake@input[["counts"]], header=TRUE,
	row.names="gene", check.names=FALSE
)
counts$A3N1 <- NULL

counts <- as.matrix(counts)
mode(counts) <- "integer"


# Loading samples information from file.
# IMPORTANT -> Replicate A3N1 was removed due to being an outlier
#             (low HNF4a expression)
all_samples <- read.table(
    snakemake@input[["samples"]], header=TRUE,
    row.names="sample", check.names=FALSE
)


dir.create(snakemake@output[["results"]], showWarnings=FALSE)
# Looping through the samples
for (i in 1:12) {
	exp <- sprintf("A%s", i)
	
	# Slicing 
	samples <- subset(all_samples, condition=='A0' | condition==exp)
	count <- counts[,c(row.names(samples))]


	# Calculating DESeq2
	dds <- DESeqDataSetFromMatrix(
    		countData=count,
    		colData=samples,
    		design= ~condition
	)


	dds$condition <- relevel(dds$condition, ref='A0')
	dds <- DESeq(dds)
	results = results(dds, coef=exp)

    	# Writing results to file
	fname <- paste(
		snakemake@output[["results"]],
		paste(exp, "csv", sep='.'), sep='/'
	)



	write.csv(
	        as.data.frame(results),
	        file=fname,
	        quote=FALSE
	)

}
