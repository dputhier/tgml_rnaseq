rule fastqc_raw:
        input:  fwd="input/fastq/{smp}_R1.fq.gz", rev="input/fastq/{smp}_R2.fq.gz"
        threads: 1
        output: fwd="output/fastqc_raw/{smp}/{smp}_R1_fastqc/fastqc_data.txt", \
                rev="output/fastqc_raw/{smp}/{smp}_R2_fastqc/fastqc_data.txt"
        message: """--- Quality check of raw data with Fastqc."""
        shell: "fastqc --outdir  output/fastqc_raw/{wildcards.smp} --extract  -f fastq {input.fwd} {input.rev} &> {output.fwd}.log "
