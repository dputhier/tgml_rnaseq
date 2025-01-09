rule fastqc_trim_pe:
        input:  fwd="output/trimmed/{smp}.fq.gz"
        threads: 1
        params: mem="5G"
        output: fwd="output/fastqc_trim/{smp}/{smp}_fastqc/fastqc_data.txt"
        message: """--- Quality check of trimmed data with Fastqc."""
        shell: """
		module load fastqc/0.12.1
		fastqc --outdir  output/fastqc_trim/{wildcards.smp} --extract  -f fastq {input.fwd} &> {output.fwd}.log """
