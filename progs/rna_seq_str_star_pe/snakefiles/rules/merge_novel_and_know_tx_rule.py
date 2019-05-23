
rule merge_novel_and_known:
    input: novel="output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + "_with_gene_name.gtf", \
    known=config["gtf"]
    output: "output/inferred_gene_annotation/known_transcripts_and_selected_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"
    threads: 1
    message: "--- Merging known and novel transcripts."
    shell:"""
    cat {input.novel} {input.known} > {output}
    """ 