rule fCounts_known_novel:
    input:  gtf="output/inferred_gene_annotation/known_transcripts_and_selected_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf",\
            bam=expand("output/bam/{smp}.bam", smp=config["samples"].split())
    output: "output/quantification_known_and_novel_genes/gene_counts_known_and_novel.txt",
            "output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt"
    threads: config["featureCounts"]["threads"]
    params: mem="10G" , strand=config["featureCounts"]["strand"]
    shell: """
    module load subread/2.0.6
    featureCounts  {params.strand} -T {threads} -t exon -g gene_name -a {input.gtf} -o {output[0]} {input.bam} &> {output[0]}.log
    cut -f 1,7- {output[0]}| awk 'NR > 1' | awk '{{gsub("output/bam/","",$0); print}}' | \
    awk '{{if(NR==1){{gsub(".bam","",$0); print}}else{{print}} }}' > {output[1]} 
"""
