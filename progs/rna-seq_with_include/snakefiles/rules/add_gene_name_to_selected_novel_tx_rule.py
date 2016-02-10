rule add_gene_name_to_unknown:
    input: "output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"
    output: "output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + "_with_gene_name.gtf"
    threads: 1
    message: "--- Adding gene name to novel transcript."
    run:
        import re 
        fh_in = open(input[0], "r")
        fh_out = open(output[0], "w")
        for line in fh_in:
            line = line.rstrip("\n")
            if not re.search("gene_name", line):
                gene_id = re.match('.*gene_id "(.*?)"', line).group(1)
                fh_out.write(line + ' gene_name "' + gene_id + '";\n')