rule tophat_pe:
  input:    fwd="output/trimmed/{smp}_R1_t.fq.gz", 
            rev="output/trimmed/{smp}_R2_t.fq.gz", 
            single="output/trimmed/{smp}_R1_singletons.fq.gz"
  output:   "output/bam/{smp}.bam"
  params:   gtf=config["gtf"], args=config["tophat"]["args"], index=config["tophat"]["index"]
  threads:  config["tophat"]["threads"]
  message:  """--- mapping with tophat (pe)."""
  shell:    """
                mkdir -p output/tophat/{wildcards.smp}
                tophat2                                          \
                        -o output/tophat/{wildcards.smp}         \
                        {params.args}                            \
                        -p {threads}                             \
                        -G {params.gtf}                          \
                        {params.index}                           \
                        {input.fwd} {input.rev} &> output/tophat/{wildcards.smp}/run_tophat.log
                cd output/tophat/{wildcards.smp}
                mv accepted_hits.bam ../../bam/{wildcards.smp}.bam
            """
  