
rule samtools_sort:
        input: "output/bam/{smp}.bam"
        output: "output/bam/{smp}_srt.bam"
        params: index = config["bowtie2"]["index"], \
                options = config["bowtie2"]["other_options"]
                
        threads: config["samtools_sort"]["threads"]
        message: """--- Sorting bam files."""
        shell: """
        samtools sort -@ {threads} -m 7G {input} output/bam/{wildcards.smp}_srt
        samtools index output/bam/{wildcards.smp}_srt.bam
        """ 
