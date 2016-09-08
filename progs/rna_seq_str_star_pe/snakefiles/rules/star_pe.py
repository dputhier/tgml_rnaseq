
def getchr(wildcards):
    cnf = config["chrom_list"]
    return cnf.replace(","," ")
    

rule star_pe:
  input:    fwd="output/trimmed/{smp}_R1_t.fq.gz", 
            rev="output/trimmed/{smp}_R2_t.fq.gz"
  output:   "output/bam/{smp}.bam"
  params:   index=config["star"]["index"], args=config["star"]["args"], gtf=config["gtf"], chrom_list=getchr
  threads:  config["star"]["threads"]
  message:  """--- mapping with star (pe)."""
  shell:    """
        mkdir -p output/star/{wildcards.smp}
        cd output/star/{wildcards.smp}               
        STAR  --genomeDir {params.index} --readFilesCommand gunzip -c --readFilesIn ../../../{input.fwd} ../../../{input.rev} --runThreadN {threads} --sjdbGTFfile {params.gtf}  {params.args} 
        samtools index Aligned.sortedByCoord.out.bam
        samtools view -h -b Aligned.sortedByCoord.out.bam {params.chrom_list} > Aligned.out_chr.bam
        samtools sort -@ {threads} Aligned.out_chr.bam Aligned.out_chr_sorted
        mv Aligned.out_chr_sorted.bam ../../bam/{wildcards.smp}.bam
        samtools index ../../bam/{wildcards.smp}.bam
        rm -f Aligned.out.bam* Aligned.sortedByCoord.out.bam*

            """