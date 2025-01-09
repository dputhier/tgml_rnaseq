rule create_bigwig:
    input: min="output/bam/{smp}_min.bam", plus="output/bam/{smp}_plus.bam"
    output: min="output/bwig/{smp}_min.bw", plus="output/bwig/{smp}_plus.bw"
    threads: 1
    params: chr=CHR, mem="5G"
    shell: """
        module load rseqc/2.6.4
        module load ucsc-wigtobigwig/377

        bam2wig.py-i {input.plus} -s {params.chr} -o output/bwig/{wildcards.smp}_plus &> {output.plus}.log
        wigToBigWig -clip output/bwig/{wildcards.smp}_plus.wig {params.chr} {output.plus} 2>&1 >> {output.plus}.log
        rm -f output/bwig/{wildcards.smp}_plus.wig
        
        bam2wig.py-i {input.min} -s {params.chr} -o output/bwig/{wildcards.smp}_min&> {output.min}.log
        wigToBigWig -clip output/bwig/{wildcards.smp}_min.wig {params.chr} {output.min} 2>&1 >> {output.min}.log
        rm -f output/bwig/{wildcards.smp}_min.wig
    """
