rule plotFingerprint:
    input:  bam= expand("output/bam/{smp}_q30_rmDup.bam", smp=config["samples"].split())
    output: tab="output/plotFingerprint/fingerprint.tab", \
            png="output/plotFingerprint/fingerprint.png"
    params: frag_size=config["fragment_size"], labels=" ".join(list([config["samples"]]))
    shell: """
    plotFingerprint \
     -b {input.bam} \
    --labels {params.labels} \
    --skipZeros \
    --numberOfSamples 50000 \
    -T "Fingerprints of the samples"  \
    --plotFile {output.png} \
    --outRawCounts {output.tab} 
 """
 