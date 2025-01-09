rule mapping_stats_R1:
        input: raw="output/fastqc_raw/{smp}/{smp}_fastqc/fastqc_data.txt", \
               trim="output/fastqc_trim/{smp}/{smp}_fastqc/fastqc_data.txt", \
               bam="output/bam/{smp}.bam"
        output: s="output/mapping_stats/{smp}.stats", fs="output/mapping_stats/{smp}.flagstat"
        message: """--- Performing some stats about mapping."""
        params: mem="5G"
        threads: 1
        shell:  """
                module load samtools/1.21
                echo "RAW"                              > {output.s}
                grep "Total Sequence" {input.raw}       >> {output.s}
                echo "TRIM"                             >> {output.s}
                grep "Total Sequence" {input.trim}      >> {output.s}
                echo "MAPPED"                           >> {output.s}
                samtools flagstat {input.bam}           > {output.fs}
                grep "mapped (" {output.fs}             >> {output.s}
                """
