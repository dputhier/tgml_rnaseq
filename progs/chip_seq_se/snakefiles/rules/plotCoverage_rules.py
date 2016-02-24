rule plotCoverage:
    input:  bam= expand("output/bam/{smp}_q30_rmDup.bam", smp=config["samples"].split())
    output: tab="output/plotCoverage/coverage.tab", \
            png="output/plotCoverage/coverage.png"
    params: labels=" ".join(list([config["samples"]]))
    shell: """
    plotCoverage \
    -b {input.bam} \
    --plotFileFormat png \
    --labels {params.labels} \
    -n 1000000 \
    --plotTitle "Coverage plot" \
    --ignoreDuplicates \
    --minMappingQuality 10 \
    --plotFile {output.png} \
    --outRawCounts {output.tab} 2> {output.tab}.log
 """
 