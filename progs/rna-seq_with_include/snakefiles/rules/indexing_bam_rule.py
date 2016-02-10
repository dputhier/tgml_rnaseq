rule index_bam:
    input: full="output/bam/{smp}.bam", min="output/bam/{smp}_min.bam", plus="output/bam/{smp}_min.bam"
    output: full="output/bam/{smp}.bam.bai", min="output/bam/{smp}_min.bam.bai", plus="output/bam/{smp}_min.bam.bai"
    threads: 1
    shell: """
    samtools index {input.full}
    samtools index {input.min}
    samtools index {input.plus}
    """