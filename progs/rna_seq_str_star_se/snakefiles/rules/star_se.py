rule STAR:
  input:    fwd="output/trimmed/{smp}.fq.gz"
  output:   "output/bam/{smp}.bam"
  params:   index=config["star"]["index"], args=config["star"]["args"], gtf=config["gtf"], mem="40G"
  threads:  config["star"]["threads"]
  message:  """--- mapping with star (pe)."""
  shell:    """
        mkdir -p output/star/{wildcards.smp}
        cd output/star/{wildcards.smp}   
        module load star/2.7.5a
        module load samtools/1.21
        STAR  --genomeDir {params.index} --readFilesCommand gunzip -c --readFilesIn ../../../{input.fwd}  --runThreadN {threads} --sjdbGTFfile {params.gtf}  {params.args} 
        mv Aligned.sortedByCoord.out.bam ../../bam/{wildcards.smp}.bam
        samtools index ../../bam/{wildcards.smp}.bam
        rm -f Aligned.out.bam Aligned.out.bam
            """
