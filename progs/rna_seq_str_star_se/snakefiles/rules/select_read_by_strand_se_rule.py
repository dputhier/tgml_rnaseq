rule select_reads_by_strand:
    input:  "output/bam/{smp}.bam"
    output: min="output/bam/{smp}_min.bam", plus="output/bam/{smp}_plus.bam"
    threads: 1
    params: mem="5G"
    shell: """
        module load samtools/1.21
        samtools view -f16 -hb {input} > {output.plus}
        samtools view -F16 -hb {input} > {output.min}
    """