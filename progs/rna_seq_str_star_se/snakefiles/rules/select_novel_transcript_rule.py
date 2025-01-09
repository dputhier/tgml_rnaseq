rule select_novel_transcripts:
    input: "output/cuffmerge/merged.gtf"
    output: "output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"
    params: selected_class_code=config["cuffmerge"]["selected_class_code"], mem="5G"
    threads: 1
    message: "--- Selecting novel_transcripts."
    shell:"""
    awk '/class_code "[{params.selected_class_code}]"/' {input} > {output}
    """
