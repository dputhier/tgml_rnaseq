rule mapping_stats_R2:
        input: raw="output/fastqc_raw/{smp}/{smp}_R2.fq_fastqc/fastqc_data.txt", \
               trim="output/fastqc_trim/{smp}/{smp}_R2_t.fq_fastqc/fastqc_data.txt", \
               bam="output/bam/{smp}.bam"
        output: s="output/mapping_stats/{smp}_R2.stats", fs="output/mapping_stats/{smp}_R2.flagstat"
        message: """--- Performing some stats about mapping."""
        threads: 1
        shell:  """
                echo "RAW"                              > {output.s}
                grep "Total Sequence" {input.raw}       >> {output.s}
                echo "TRIM"                             >> {output.s}
                grep "Total Sequence" {input.trim}      >> {output.s}
                echo "MAPPED"                           >> {output.s}
                samtools flagstat {input.bam}           > {output.fs}
                grep "read1" {output.fs}                >> {output.s}
                """
