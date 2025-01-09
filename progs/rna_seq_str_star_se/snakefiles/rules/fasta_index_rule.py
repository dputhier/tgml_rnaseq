rule index_fasta:
    input:  config["fasta"]
    output: config["fasta"] + ".fai"
    threads: 1
    params: mem="4G"
    message: "--- Indexing fasta file."
    shell:  """
        samtools faidx {input} &> {output}.log
    """
