"""
	Check the quality of fastq-formatted raw reads using the program fastQC (quality control).

"""

rule fastqc_trim:
        input: "output/trimmed/{smp}_t.fq.gz"
        output: "output/fastqc_trim/{smp}/{smp}_t.fq_fastqc/fastqc_data.txt"
        message: """--- Quality check of raw data with Fastqc."""
        shell: "fastqc --outdir output/fastqc_trim/{wildcards.smp} --extract  -f fastq {input} &> {output}.log "
                        
