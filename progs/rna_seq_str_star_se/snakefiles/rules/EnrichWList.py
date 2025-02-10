import os


script = os.path.join(config["workingdir"], "progs", "scripts", "EnrichWList.R")

rule enrichlist:
    input:  gene="output/cuffdiff/{comp}/{comp}_gene_exp.diff", liste="output/LIST_{comp}.xlsx"
    output: "output/cuffdiff/{comp}/FuncEnrichWithList/barplot.png", \
            "output/cuffdiff/{comp}/FuncEnrichWithList/dotplot.png", \
            "output/cuffdiff/{comp}/FuncEnrichWithList/upsetplot.png"
    params: script=script, mem="4G"
    threads: 4
    message: "--- Functionnal enrichment -----"
    shell:  """
        module unload r
        module load r/4.4.1
        Rscript {params.script} -i {input.gene} -l {input.liste}
    """
