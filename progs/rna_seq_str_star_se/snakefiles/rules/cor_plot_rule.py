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
    params: smp=",".join(config["samples"].split(" ")), script=script, mem="4G"
    threads: 1
    message: "--- Producing correlation plots."
    shell:  """
        module load r/4.4.1
        {params.script} -i {input} -o output/corr_plot  -w 12 -y 12 -t 1 2> {output[0]}.log
        touch {output}
    """
