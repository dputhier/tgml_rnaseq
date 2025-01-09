import os
total_bam_sum = os.path.join(config["workingdir"], "progs", "scripts", "total_bam_sum.py")

rule do_bigwig:
    input:  min="output/bam/{smp}_min.bam", plus="output/bam/{smp}_plus.bam", \
            minbai="output/bam/{smp}_min.bam.bai", plusbai="output/bam/{smp}_plus.bam.bai"
    output: min="output/bwig/{smp}_min.bw", plus="output/bwig/{smp}_plus.bw", \
            min_norm="output/bwig/{smp}_min_norm.bw", plus_norm="output/bwig/{smp}_plus_norm.bw"
    params: chr=config["chr_size"], bs=total_bam_sum, sf=config["total_cov_objective"], mem="5G"
    threads: 1
    shell: """
        module load rseqc/2.6.4
        module load ucsc-wigtobigwig/377

        bam2wig.py -t {params.sf} -i {input.plus} -s {params.chr} -o output/bwig/{wildcards.smp}_plus_norm &> {output.plus_norm}.log
        wigToBigWig -clip output/bwig/{wildcards.smp}_plus_norm.wig {params.chr} {output.plus_norm} 2>&1 >> {output.plus}.log
        rm -f output/bwig/{wildcards.smp}_plus_norm.wig


        bam2wig.py -i {input.plus} -s {params.chr} -o output/bwig/{wildcards.smp}_plus &> {output.plus}.log
        wigToBigWig -clip output/bwig/{wildcards.smp}_plus.wig {params.chr} {output.plus} 2>&1 >> {output.plus}.log
        rm -f output/bwig/{wildcards.smp}_plus.wig
        
        
        bam2wig.py -t {params.sf}  -i {input.min} -s {params.chr} -o output/bwig/{wildcards.smp}_min_norm &> {output.min_norm}.log
        wigToBigWig -clip output/bwig/{wildcards.smp}_min_norm.wig {params.chr} {output.min_norm} 2>&1 >> {output.min_norm}.log
        rm -f output/bwig/{wildcards.smp}_min_norm.wig
        

        bam2wig.py  -i {input.min} -s {params.chr} -o output/bwig/{wildcards.smp}_min &> {output.min}.log
        wigToBigWig -clip output/bwig/{wildcards.smp}_min.wig {params.chr} {output.min} 2>&1 >> {output.min}.log
        rm -f output/bwig/{wildcards.smp}_min.wig
    """