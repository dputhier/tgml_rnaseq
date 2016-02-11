rule rseqc_coverage_diag:
    input: expand("output/bam/{smp}.bam", smp=config["samples"].split())
    params: bed=config["housekeeping"]
    output: "output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.curves.pdf.png"
    threads: 1
    message: "--- Genebody coverage diagnosis ---"
    shell:"""
    ls -1 output/bam/* | grep -v "_min.bam" | grep -v "_plus.bam" |grep -v ".bam.bai" > output/rseqc_coverage_diag/bam_list.txt
     geneBody_coverage.py -r {params.bed} -i output/rseqc_coverage_diag/bam_list.txt -o output/rseqc_coverage_diag/rseqc_cov &> {output}.log
     touch {output}
     convert output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.curves.pdf {output}
    """