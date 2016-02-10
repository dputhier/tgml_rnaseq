rule select_reads_by_strand:
    input:  "output/bam/{smp}.bam"
    output: min="output/bam/{smp}_min.bam", plus="output/bam/{smp}_plus.bam"
    threads: 1
    shell: """
        samtools view -f99 -hb {input} > {output.min}_99.bam
        samtools view -f147 -hb {input} > {output.min}_147.bam
        samtools view -f83 -hb {input} > {output.plus}_83.bam
        samtools view -f163 -hb {input} > {output.plus}_163.bam
        
        samtools merge -f -h {output.min}_99.bam {output.min} {output.min}_99.bam {output.min}_147.bam
        samtools merge -f -h {output.plus}_83.bam {output.plus} {output.plus}_83.bam {output.plus}_163.bam
        
        rm -f {output.min}_99.bam {output.min}_147.bam {output.plus}_83.bam {output.plus}_163.bam
    """