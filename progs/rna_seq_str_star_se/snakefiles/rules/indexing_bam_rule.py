rule index_bam:
    input: full="output/bam/{smp}.bam", min="output/bam/{smp}_min.bam", plus="output/bam/{smp}_plus.bam"
    output: full="output/bam/{smp}.bam.bai", min="output/bam/{smp}_min.bam.bai", plus="output/bam/{smp}_plus.bam.bai"
    threads: 1
    params: mem="5G"
    shell: """
    module load samtools/1.21
    samtools index {input.full}
    samtools index {input.min}
    samtools index {input.plus}
    """