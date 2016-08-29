rule select_reads_by_strand:
    input:  "output/bam/{smp}.bam"
    output: min="output/bam/{smp}_min.bam", plus="output/bam/{smp}_plus.bam"
    threads: 1
    shell: """
        samtools view -f16 -hb {input} > {output.plus}
        samtools view -F16 -hb {input} > {output.min}
    """