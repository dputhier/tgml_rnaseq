import os

script = os.path.join(config["workingdir"], "progs", "scripts", "corrPlot.R")

rule corrplot:
    input:  "output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt"
    output: "output/corr_plot/CoorPlot_circle.png", \
            "output/corr_plot/CoorPlot_ellipe.png", \
            "output/corr_plot/CoorPlot_pie.png",    \
            "output/corr_plot/CoorPlot_piewb.png",  \ 
            "output/corr_plot/CoorPlot_square.png", \
            "output/corr_plot/Pairs_plot.png"
    params: smp=",".join(config["samples"].split(" ")), script=script
    conda: os.path.join(config["workingdir"], "conda", "R.yaml")
    threads: 1
    conda: os.path.join(config["workingdir"], "conda", "R.yaml")
    message: "--- Producing correlation plots."
    shell:  """
        {params.script} -i {input} -o output/corr_plot  -t 1 2> {output[0]}.log
        touch {output}
    """
