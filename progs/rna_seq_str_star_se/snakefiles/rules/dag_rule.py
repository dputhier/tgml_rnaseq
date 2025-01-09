rule draw_dag:
    input: "progs/rna_seq_str_star_se/snakefiles/Snakefile.py"
    output: "output/report/dag.png"
    threads: 1
    params: mem="4G"
    shell: """
    cd progs/rna_seq_str_star_se
    mkdir -p ../../output/report/ 
    snakemake -s snakefiles/Snakefile.py --rulegraph | dot -T png > ../../output/report/rulegraph.png
    snakemake -s snakefiles/Snakefile.py --dag | dot -T png > ../../output/report/dag.png
    """