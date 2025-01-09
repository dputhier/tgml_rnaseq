rule normalize_counts:  
    input: "output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt"
    output: "output/normalize_counts/gene_counts_known_and_novel_mini_log2_pseudocount_norm.txt"
    threads: 1
    params: mem="4G"
    shell: """
    module load r/4.4.1
    echo "library(DESeq2); dir.create('output/normalize_counts', showWarnings = FALSE); \
    d <- read.table('{input}', head=TRUE, sep='\\t', row.names = 1); \
    pheno <- data.frame(smp=colnames(d), condition=as.factor(sample(c(1,2), ncol(d), replace = TRUE))); \
    rownames(pheno) <- colnames(d) ; \
    des <- DESeq2::DESeqDataSetFromMatrix(d, pheno, design=formula(~condition)) ; \
    dds <- DESeq(des) ; \
    log2_counts <- log2(counts(dds, normalized=TRUE) + 1) ; \
    write.table(log2_counts,file='{output}' , quote = FALSE, sep = '\\t', col.names = NA)" | R --slave
"""