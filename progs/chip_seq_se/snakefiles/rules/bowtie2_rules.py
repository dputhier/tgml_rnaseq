
rule mapping:
        input: "output/trimmed/{smp}_t.fq.gz"
        output: "output/bam/{smp}.bam"
        params: index = config["bowtie2"]["index"], \
                options = config["bowtie2"]["other_options"]
                
        threads: config["bowtie2"]["threads"]
        message: """--- Mapping samples with Bowtie2."""
        shell: """
            bowtie2 {params.options} -q -x {params.index} -U {input} 2> {output}.log | samtools view -bS - > {output}
        """ 
