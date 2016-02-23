"""
    Definition of the sorting rule bam file
"""

rule sorting_bam:
    input: "output/bam/{smp}.bam"
    output: "output/bam/{smp}_sort.bam"
    message: """--- Sorting mapped reads by position with samtools."""
    shell: "samtools sort {input} output/bam/{wildcards.smp}_sorted"
