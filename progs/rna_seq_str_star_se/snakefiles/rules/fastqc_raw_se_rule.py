rule fastqc_raw:
        input:  fwd="input/fastq/{smp}.fq.gz"
        threads: 1
        params: mem="5G"
        output: fwd="output/fastqc_raw/{smp}/{smp}_fastqc/fastqc_data.txt"
        message: """--- Quality check of raw data with Fastqc."""
        shell: """
		module load fastqc/0.12.1
		fastqc --outdir  output/fastqc_raw/{wildcards.smp} --extract  -f fastq {input.fwd}  &> {output.fwd}.log 
		"""
