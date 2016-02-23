rule draw_dag:
    input: "progs/rna_seq_str_pe/snakefiles/Snakefile.py"
    output: "output/report/dag.png"
    threads: 1
    shell: """
    cd progs/rna_seq_str_pe
    snakemake -s snakefiles/Snakefile.py --rulegraph | dot -T png > ../../output/report/rulegraph.png
    snakemake -s snakefiles/Snakefile.py --dag | dot -T png > ../../output/report/dag.png
    """