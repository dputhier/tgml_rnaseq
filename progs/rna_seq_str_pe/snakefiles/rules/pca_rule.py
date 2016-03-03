from __future__ import print_function
import os

pheno_file = open(os.path.join(config["workingdir"], "output/pheno.txt"), "w")

cfg_pheno = config["pheno"].split(" ")
cfg_smp = config["samples"].split(" ")

for smp,cnd  in zip(cfg_smp, cfg_pheno):
    pheno_file.write(smp + "\t" + cnd + "\n")

pheno_file.close()

script = os.path.join(config["workingdir"], "progs", "scripts", "pca_plot.R")

pheno = config["pheno"].split()
pheno=",".join(pheno)

rule pca_mds:
    input:  "output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt"
    output: "output/pca_mds/done"
    params: pheno="output/pheno.txt", script=script
    threads: 1
    message: "--- Performing PCA and MDS analysis."
    shell:  """
        {params.script} -i {input} -f {params.pheno} -o output/pca_mds  -q -z 2> {output[0]}.log
        touch {output}
        mv {params.pheno} output/pca_mds
    """
