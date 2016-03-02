rule index_fasta:
    input:  "output/quantification_known_and_novel_genes/gene_counts_known_and_novel.txt"
    output: 
    threads: 1
    message: "--- Performing DESeq analysis."
    shell:  """
        
    """
