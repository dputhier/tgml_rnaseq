rule mapping_stats:
        input: raw="output/fastqc_raw/{smp}/{smp}.fq_fastqc/fastqc_data.txt", \
                trim="output/fastqc_trim/{smp}/{smp}_t.fq_fastqc/fastqc_data.txt", \
                bam="output/bam/{smp}_srt.bam", \
                bamq30="output/bam/{smp}_q30.bam", \
                bamRmDup="output/bam/{smp}_q30_rmDup.bam"
        
        output: s="output/mapping_stats/{smp}.stats", \
            fs="output/mapping_stats/{smp}.flagstat", \
            fsq30="output/mapping_stats/{smp}q30.flagstat",\
            fsrm="output/mapping_stats/{smp}q30_rmDup.flagstat"

        message: """--- Performing some stats about mapping."""
        threads: 1
        shell:  """
            echo "RAW"                                            > {output.s}
            grep "Total Sequence" {input.raw}                     >> {output.s}
            echo "TRIM"                                           >> {output.s}
            grep "Total Sequence" {input.trim}                    >> {output.s}
            echo "MAPPED"                                         >> {output.s}
            samtools flagstat {input.bam}                         > {output.fs}
            grep "mapped" {output.fs} | head -1| sed 's/ .*//'    >> {output.s}
            echo "QUALITY_FILTERED"                               >> {output.s}
            samtools flagstat {input.bamq30}         > {output.fsq30}
            grep "mapped" {output.fsq30} | head -1| sed 's/ .*//' >> {output.s}
            echo "RM_DUP"                                         >> {output.s}
            samtools flagstat {input.bamRmDup}                    > {output.fsrm}
            grep "mapped" {output.fsrm} | head -1| sed 's/ .*//'   >> {output.s}
        """
