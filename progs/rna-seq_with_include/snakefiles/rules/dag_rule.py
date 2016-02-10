rule draw_dag:
    input: "progs/rna-seq_with_include/snakefiles/Snakefile.py"
    output: "output/report/dag.png"
    threads: 1
    shell: """
    cd progs/rna-seq_with_include/
    snakemake -s snakefiles/Snakefile.py --rulegraph | dot -T png > ../../output/report/dag.png
    """