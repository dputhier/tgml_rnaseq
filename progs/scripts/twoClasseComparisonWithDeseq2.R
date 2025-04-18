#!/usr/bin/Rscript 

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
    'help',           	'h', 0, "logical",	"Help about the program",
    'input_file',     	'i', 1, "character",	"REQUIRED: tabulated flate file). E.g output from featureCounts.",
    'pheno',		'p', 1, "character", 	"REQUIRED: comma separated list of phenotypes.",
    'outdir',           'o', 1, "character", 	"Output directory. Default to current working directory.",
    'skip',		's', 1, "integer", 	"The number of line to skip in the matrix files.",
    'del',		'd', 1, "character",	"The comma separated list of column to delete in the matrix file."
			), byrow=TRUE, ncol=5);

opt = getopt(spec)




# if help was asked, print a friendly message
# and exit with a non-zero error code
args <- commandArgs()
if ( !is.null(opt$help) | length(args) < 5) {
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
pheno <- unlist(strsplit(opt$pheno, ","))

## columns to delete
if(!is.null(opt$del)){
	to_del <- as.integer(unlist(strsplit(opt$del, ",")))
	to_del <- to_del - 1 # row names is not counted 
	m  <- read.table(opt$input_file, sep='\t', head=TRUE, row=1, quote='', comment.char="", skip=1, )[,-to_del]
}else{
	m  <- read.table(opt$input_file, sep='\t', head=TRUE, row=1, quote='', comment.char="", skip=1, )
}


cat("Calling differentially expressed genes (DESeq2\n")
des <- DESeqDataSetFromMatrix(m, as.data.frame(pheno), design=formula(~pheno))
dds <- DESeq(des, fitType='local')
resMLE<- results(dds, addMLE=TRUE, cooksCutoff=FALSE)
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


df <- data.frame(	resMLE$baseMean, 
					resMLE$lfcMLE, 
					ifelse(is.na(res$padj), 
					FALSE, res$padj < .1)
	)


## -----------------------------------------------------------------------------
## Output
## -----------------------------------------------------------------------------


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

res <- results(dds)
resOrdered <- res[order(res$padj),]
write.table(as.data.frame(resOrdered), 
			file.path(opt$outdir,
			"DESeq2_diff_genes.txt"), sep="\t", col.names=NA, quote=F)
