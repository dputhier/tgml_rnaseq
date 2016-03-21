rule bwig:
        input: "output/bam/{smp}_q30_rmDup.bam"
        output: "output/bwig/{smp}_q30_rmDup.bigwig"
        message: """--- Delete replicated reads with samtools."""
        params: gs = config["genome_size"], frag_len=config["fragment_size"]
        threads: config["bam_coverage"]["threads"]
        shell: """
            bamCoverage --bam {input} \
                --binSize 10 --extendReads {params.frag_len} \
            --numberOfProcessors {threads} --normalizeTo1x {params.gs} \
            -o  {output} 2> {output}.log
        """



