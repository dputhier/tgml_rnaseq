rule star_pe:
  input:    fwd="output/trimmed/{smp}_R1_t.fq.gz", 
            rev="output/trimmed/{smp}_R2_t.fq.gz"
  output:   "output/bam/{smp}.bam"
  params:   index=config["star"]["index"], args=config["star"]["args"], gtf=config["gtf"] 
  threads:  config["star"]["threads"]
  message:  """--- mapping with star (pe)."""
  shell:    """
        mkdir -p output/star/{wildcards.smp}
        cd output/star/{wildcards.smp}               
        STAR  --genomeDir {params.index} --readFilesCommand gunzip -c --readFilesIn ../../../{input.fwd} ../../../{input.rev} --runThreadN {threads} --sjdbGTFfile {params.gtf}  {params.args} 
        mv Aligned.sortedByCoord.out.bam ../../bam/{wildcards.smp}.bam
        samtools index ../../bam/{wildcards.smp}.bam
        rm -f Aligned.out.bam Aligned.out.bam
            """
