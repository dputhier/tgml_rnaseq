rule index_fasta:
    input:  config["fasta"]
    output: config["fasta"] + ".fai"
    threads: 1
    message: "--- Indexing fasta file."
    shell:  """
        samtools faidx {input} &> {output}.log
    """
