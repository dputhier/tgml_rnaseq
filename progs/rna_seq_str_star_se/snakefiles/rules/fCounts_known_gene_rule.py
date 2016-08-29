rule fCounts_known_genes:
    input: gtf=config["gtf"], bam=expand("output/bam/{smp}.bam", smp=config["samples"].split()) 
    output: "output/quantification_known_genes/gene_counts.txt","output/quantification_known_genes/gene_counts_mini.txt"
    threads: 1
    shell: """
    featureCounts -s 2 -T 15 -t exon -g gene_name -a {input.gtf} -o {output[0]} {input.bam} &> {output[0]}.log
    cut -f 1,7- {output[0]}| awk 'NR > 1' | awk '{{gsub("output/bam/","",$0); print}}' | \
     awk '{{if(NR==1){{gsub(".bam","",$0); print}}else{{print}} }}' > {output[1]}> {output[1]}
"""