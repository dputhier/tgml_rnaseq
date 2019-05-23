#!/usr/bin/env Rscript 

##------------------------------------------------------------------------------
## D puthier
## Created: Tue Jul  8 15:09:18 CEST 2014
## Last modification: Wed Aug 13 17:52:46 CEST 2014
##------------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

library("getopt")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("DESeq2")
library("edgeR")

## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
    'help',        'h', 0, "logical",	"Help about the program",
    'input_file',  'i', 1, "character",	"REQUIRED: tabulated flate file). E.g output from featureCounts.",
    'outdir',      'o', 1, "character", 	"Output directory. Default to current working directory.",
    'skip',			's', 1, "integer", 	"The number of line to skip in the matrix files.",
    'class1',			'c', 1, "character",	"The sample to select for class 1.",
	'class2',			'd', 1, "character",	"The sample to select for class 2."), byrow=TRUE, ncol=5);

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$class1) | is.null(opt$class2)) {
        cat("Perform differential expression call with DESeq2\n")
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## -----------------------------------------------------------------------------
## Processing data
## -----------------------------------------------------------------------------

# output
if(is.null(opt$outdir))
	opt$outdir <- getwd()

dir.create(opt$outdir, showWarnings=FALSE)

## phenotype
class1 <- unlist(strsplit(opt$class1, ","))
class2 <- unlist(strsplit(opt$class2, ","))
pheno <- c(rep("class1", length(class1)),
			rep("class2", length(class2)))

## Gete selected samples

m  <- read.table(opt$input_file, sep='\t', 
		head=TRUE, row=1, quote='', 
		comment.char="", check.names=FALSE)

if(!all(c(class1, class2) %in% colnames(m))){
	cat("ERROR: Check column names. Unknow sample selected.\n")
	q(status=1)
}

m <- m[, c(class1, class2)]


cat("Calling differentially expressed genes (DESeq2)\n")
des <- DESeqDataSetFromMatrix(m, as.data.frame(pheno), design=formula(~pheno))
dds <- DESeq(des)
res <- results(dds, cooksCutoff=FALSE)

## -----------------------------------------------------------------------------
## Diagnostic plots
## -----------------------------------------------------------------------------

png(file.path(opt$outdir,"DESeq2_diagnostic_disp.png"), 
		width = 6, height = 6, units = 'in', res = 300)

## Dispersion plot

plotDispEsts(dds)
dev.off()

## MA plot

cat("Producing MA plot.\n")
png(file.path(opt$outdir,"DESeq2_diagnostic_MA.png"),
		width = 6, height = 6, units = 'in', res = 300)

plotMA(	res, 
		main="MA plot of Two conditions",
		ylim=c(-2,2)
	)
dev.off()




## -----------------------------------------------------------------------------
## Output
## -----------------------------------------------------------------------------

print(colnames(m))
cat("Writing tables (raw counts).\n")
write.table(m, 
		file.path(opt$outdir,
		"DESeq2_raw_count_table.txt"), 
		sep='\t', quote=F, 
		col.names=NA) 

cat("Writing tables (normalized counts).\n")
norm.counts <- counts(dds, normalized=TRUE)
colnames(norm.counts) <- colnames(m)
write.table(norm.counts, 
		file.path(opt$outdir,
				"DESeq2_norm_count_table_lin.txt"), 
		sep='\t', quote=F, 
		col.names=NA) 

cat("Writing tables (log2-transformed normalized counts + 1).\n")
log2.counts <- log2(counts(dds, normalized=TRUE) + 1)
colnames(log2.counts) <- colnames(m)
write.table(log2.counts, 
		file.path(opt$outdir,
				"DESeq2_norm_count_table_log2.txt"), 
		sep='\t', quote=F, 
		col.names=NA) 

cat("Writing tables (rlog transformation). Can be used for clustering of samples\n.")
rld <- assay(rlog(dds))
colnames(rld) <- colnames(m)
write.table(rld, 
		file.path(opt$outdir,
				"DESeq2_norm_count_table_rld.txt"), 
		sep='\t', quote=F, 
		col.names=NA) 

cat("Writing tables (vsd transformation). Can be used for clustering of samples\n.")
vsd <- assay(varianceStabilizingTransformation(dds))
colnames(vsd) <- colnames(m)
write.table(vsd, 
		file.path(opt$outdir,
				"DESeq2_norm_count_table_vsd.txt"), 
		sep='\t', quote=F, 
		col.names=NA) 

resOrdered <- res[order(res$padj),]
#save.image("image.Rdata")
out <- as.data.frame(resOrdered)
tmp    <- sapply(c("BH","BY","holm","bonferroni","fdr"), 
		function(meth) p.adjust(out$pvalue, meth))
out <- data.frame(out,tmp)
colnames(out) <- paste("DESeq2_",colnames(out), sep="")

write.table(data.frame(out, log2.counts[rownames(out),], check.names = FALSE), 
			file.path(opt$outdir,
			"DESeq2_pval_and_norm_count_log2.txt"), sep="\t", col.names=NA, quote=F)


## -----------------------------------------------------------------------------
## EdgeR
## -----------------------------------------------------------------------------

cat("-- EdgeR analysis")
y <- DGEList(counts=round(m,0),group=pheno)
y <- calcNormFactors(y)
design <- model.matrix(~pheno)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)

# OUtput
qlf.sort <- qlf$table[rownames(out),]
tmp    <- sapply(c("BH","BY","holm","bonferroni","fdr"), 
		function(meth) p.adjust(qlf.sort$PValue, meth))
qlf.sort <- data.frame(qlf.sort, tmp)
colnames(qlf.sort) <- paste("EdgeR_",colnames(qlf.sort), sep="")


out <- data.frame(out, qlf.sort, log2.counts[rownames(out),], check.names = FALSE)

write.table(out, 
		file.path(opt$outdir,
				"DESeq2_EdgeR_pval_and_norm_count_log2.txt"), sep="\t", col.names=NA, quote=F)

