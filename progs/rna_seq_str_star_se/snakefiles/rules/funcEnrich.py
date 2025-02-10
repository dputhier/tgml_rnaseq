import os

script = os.path.join(config["workingdir"], "progs", "scripts", "funcEnrich.R")

rule funcenrich:
    input:  "output/comparison/{comp}/{comp}_DESeq2_pval_and_norm_count_log2.txt"
    output: "output/comparison/{comp}/FuncEnrich/{comp}_barplot_up.png"
    params: 
        script=script, 
        mem="40G"
    threads: 30
    message: "--- Functionnal enrichment -----"
    log: "output/comparison/{comp}/FuncEnrich/{comp}.log"
    shell:  """
        module unload r
        module load r/4.4.1
        Rscript {params.script} -i {input} > {log} 2>&1
    """
