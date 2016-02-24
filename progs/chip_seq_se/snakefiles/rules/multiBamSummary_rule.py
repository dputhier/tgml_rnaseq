
rule multiBamSummary:
    input:  bam= expand("output/bam/{smp}_q30_rmDup.bam", smp=config["samples"].split()), \
            bed="output/merged_peaks/merged_peaks.bed"
    output: npz="output/multiBamSummary/merge_peaks_coverage.npz", \
            tab="output/multiBamSummary/merge_peaks_coverage.tab",
    params: frag_size=config["fragment_size"], name=" ".join(list([config["samples"]]))
    shell: """
    multiBamSummary bins \
     --bamfiles {input.bam} \
     --labels  {params.name} \
     --extendReads {params.frag_size} \
     -out {output.npz} --outRawCounts {output.tab}
 """
 