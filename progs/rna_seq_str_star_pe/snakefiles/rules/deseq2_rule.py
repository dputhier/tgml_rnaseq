import os
script = os.path.join(config["workingdir"], "progs", "scripts", "twoClasseComparisonWithDeseq2.R")

pheno = config["pheno"].split()
pheno=",".join(pheno)

rule deseq2:
    input:  "output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt"
    output: "output/diff_call_deseq2/DESeq2_diagnostic_MA.png",         \
            "output/diff_call_deseq2/DESeq2_diagnostic_disp.png",       \ 
            "output/diff_call_deseq2/DESeq2_norm_count_table_lin.txt",  \
            "output/diff_call_deseq2/DESeq2_norm_count_table_log2.txt", \
            "output/diff_call_deseq2/DESeq2_norm_count_table_rld.txt",  \
            "output/diff_call_deseq2/DESeq2_norm_count_table_vsd.txt",  \ 
            "output/diff_call_deseq2/DESeq2_raw_count_table.txt",       \
            "output/diff_call_deseq2/DESeq2_diff_genes.txt"
    params: pheno=pheno, script=script
    threads: 1
    message: "--- Performing DESeq analysis."
    shell:  """
        {params.script} -i {input} -p {params.pheno} -o output/diff_call_deseq2/ &> {output[0]}.log
    """
