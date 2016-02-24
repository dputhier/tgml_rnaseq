rule rm_dup:
        input: "output/bam/{smp}_q30.bam"
        output: "output/bam/{smp}_q30_rmDup.bam"
        message: """--- Delete replicated reads with samtools."""
        shell: """
            samtools_0.1.19 rmdup -s {input} {output} 2> {output}.log
            samtools index {output}    
        """