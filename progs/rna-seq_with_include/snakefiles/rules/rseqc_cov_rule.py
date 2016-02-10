rule rseqc_coverage_diag:
    input: expand("output/bam/{smp}.bam", smp=config["samples"].split())
    params: bam=",".join(config["samples"].split()), bed=config["housekeeping"]
    output: "output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.txt"
    threads: 1
    message: "--- Genebody coverage diagnosis ---"
    shell:"""
     geneBody_coverage.py -r {params.bed} -i {params.bam} -o output/rseqc_coverage_diag/rseqc_cov &> {output}.log
     touch {output}
    """