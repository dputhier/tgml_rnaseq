import os

rule index_fasta:
    input:  config["fasta"]
    output: config["fasta"] + ".fai"
    threads: 1
    conda: os.path.join(config["workingdir"], "conda", "rnaseq.yaml")
    message: "--- Indexing fasta file."
    shell:  """
        samtools faidx {input} &> {output}.log
    """
