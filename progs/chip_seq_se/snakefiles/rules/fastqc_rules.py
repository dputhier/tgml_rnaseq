"""
	Check the quality of fastq-formatted raw reads using the program fastQC (quality control).

"""

rule fastqc:
        input: "samples/{samples_type}/{samples}.fq.gz"
        output:"results/fastqc/{samples_type}/{samples}.fq_fastqc.zip"
        message: """--- Quality check of raw data with Fastqc."""
        shell: "fastqc --outdir results/fastqc/{wildcards.samples_type}/ fastqc/control --extract {input}"
                        
