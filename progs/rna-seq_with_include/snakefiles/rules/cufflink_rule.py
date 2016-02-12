rule cufflinks:
    input:  bam="output/bam/{smp}.bam", bai="output/bam/{smp}.bam.bai"
    output: gtf="output/cufflinks/{smp}/transcripts.gtf"
    params: gtf=config["gtf"], args=config["cufflinks"]["args"]
    threads: config["cufflinks"]["threads"]
    message: "--- Searching novel transcript with cufflinks."
    shell:  """
        cufflinks -g {params.gtf} -p {threads}  -o output/cufflinks/{wildcards.smp} {input.bam} &> {output}.log
    """
