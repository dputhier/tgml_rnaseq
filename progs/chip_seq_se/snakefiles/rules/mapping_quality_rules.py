rule select_mapped:
        input: "output/bam/{smp}_srt.bam"
        output: "output/bam/{smp}_q30.bam"
        params: quality = config["samtools_view"]["quality"]
        message: """--- Selecting mapped reads with a quality over ... with samtools."""
        shell: """
                samtools view -bh -q {params.quality} {input} > {output}
                samtools index {output}
        """  
