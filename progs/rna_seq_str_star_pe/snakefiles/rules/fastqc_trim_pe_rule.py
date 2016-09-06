
rule fastqc_trim_pe:
        input:  fwd="output/trimmed/{smp}_R1_t.fq.gz", rev="output/trimmed/{smp}_R2_t.fq.gz"
        threads: 1
        output: fwd="output/fastqc_trim/{smp}/done_fwd", \
                rev="output/fastqc_trim/{smp}/done_rev"
        message: """--- Quality check of trimmed data with Fastqc."""
        shell: """
            fastqc --outdir  output/fastqc_trim/{wildcards.smp} --extract  -f fastq {input.fwd} {input.rev} &> {output.fwd}.log 
            touch {output.fwd} {output.rev}
        """
