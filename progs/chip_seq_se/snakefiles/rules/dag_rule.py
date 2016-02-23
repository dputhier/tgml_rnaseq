rule draw_dag:
    input: "progs/chip_seq_se/snakefiles/Snakefile.py"
    output: "output/report/dag.png"
    threads: 1
    shell: """
    cd progs/chip_seq_se
    snakemake -s snakefiles/Snakefile.py --rulegraph | dot -T png > ../../output/report/rulegraph.png
    snakemake -s snakefiles/Snakefile.py --dag | dot -T png > ../../output/report/dag.png
    """