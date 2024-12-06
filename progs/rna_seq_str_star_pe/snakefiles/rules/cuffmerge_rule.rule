import os

rule cuffmerge:
    input: config["fasta"] + ".fai", expand("output/cufflinks/{smp}/transcripts.gtf", smp=config["samples"].split())
    output: "output/cuffmerge/merged.gtf"
    params: gtf=config["gtf"], fa=config["fasta"]
    threads: config["cuffmerge"]["threads"]
    conda: os.path.join(config["workingdir"], "conda", "rnaseq.yaml")
    message: "--- Comparing transcript to the reference."
    shell:"""
    ls -1 output/cufflinks/*/transcripts.gtf > output/cuffmerge/assembly.txt
    cuffmerge -o output/cuffmerge -g {params.gtf} --keep-tmp -s {params.fa} -p 5 output/cuffmerge/assembly.txt &> {output}.log
    """