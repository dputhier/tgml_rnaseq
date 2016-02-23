"""
	Check the quality of fastq-formatted raw reads using the program fastQC (quality control).

"""

rule fastqc:
        input: "input/fastq/{smp}.fq.gz"
        output: "output/fastqc_raw/{smp}/{smp}.fq_fastqc/fastqc_data.txt"
        message: """--- Quality check of raw data with Fastqc."""
        shell: "fastqc --outdir output/fastqc_raw/{wildcards.smp} --extract  -f fastq {input} &> {output}.log "
                        
