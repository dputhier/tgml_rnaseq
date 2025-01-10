rule normalize_counts:  
    input: a="output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt", \
           b= "output/quantification_known_genes/gene_counts_mini.txt"
    output: a="output/normalize_counts/gene_counts_known_and_novel_mini_log2_pseudocount_norm.txt", \
            b="output/normalize_counts/gene_counts_mini_log2_pseudocount_norm.txt"
    threads: 1
    params: mem="4G"
    shell: """
    module load r/4.4.1
    echo "library(DESeq2); dir.create('output/normalize_counts', showWarnings = FALSE); \
    d <- read.table('{input.a}', head=TRUE, sep='\\t', row.names = 1); \
    pheno <- data.frame(smp=colnames(d), condition=as.factor(sample(c(1,2), ncol(d), replace = TRUE))); \
    rownames(pheno) <- colnames(d) ; \
    des <- DESeq2::DESeqDataSetFromMatrix(d, pheno, design=formula(~condition)) ; \
    dds <- DESeq(des) ; \
    log2_counts <- log2(counts(dds, normalized=TRUE) + 1) ; \
    write.table(log2_counts,file='{output.a}' , quote = FALSE, sep = '\\t', col.names = NA); \
    d <- read.table('{input.b}', head=TRUE, sep='\\t', row.names = 1); \
    pheno <- data.frame(smp=colnames(d), condition=as.factor(sample(c(1,2), ncol(d), replace = TRUE))); \
    rownames(pheno) <- colnames(d) ; \
    des <- DESeq2::DESeqDataSetFromMatrix(d, pheno, design=formula(~condition)) ; \
    dds <- DESeq(des) ; \
    log2_counts <- log2(counts(dds, normalized=TRUE) + 1) ; \
    write.table(log2_counts,file='{output.b}' , quote = FALSE, sep = '\\t', col.names = NA)" | R --slave
"""