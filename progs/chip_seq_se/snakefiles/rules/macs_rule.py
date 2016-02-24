rule macs14:
        input: chip="output/bam/{smp}_q30_rmDup.bam", input="output/bam/" + config["input"] + "_q30_rmDup.bam"
        output: "output/macs/{smp}/{smp}_sorted_peaks.bed"
        message: """--- Searching peaks with macs."""
        params: gs=config["macs14"]["genome_size"], \
                ctrl_name=config["input"]
        shell: """
        cd output/macs/{wildcards.smp}
        macs14    -t ../../bam/{wildcards.smp}_q30_rmDup.bam \
            -c ../../bam/{params.ctrl_name}_q30_rmDup.bam -f BAM \
            -n {wildcards.smp} -g {params.gs} &> {wildcards.smp}_peaks.bed.log
        sort -g -k5,5     {wildcards.smp}_peaks.bed > {wildcards.smp}_sorted_peaks.bed
        """