from __future__ import print_function
import os
dir_base = os.path.join(config["workingdir"], "output")
if not os.path.exists(dir_base):
    os.makedirs(dir_base)
pheno_file = open(os.path.join(dir_base, "pheno.txt"), "w")

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
    output: "output/pca_mds/PCA_dim_genes_2D_samples.png",\
            "output/pca_mds/PCA_dim_genes_2D_classes.png",\
            "output/pca_mds/PCA_dim_genes_2D_overlay.png",\
            "output/pca_mds/PCA_Eigen_val.png",\
            "output/pca_mds/PCA_ggplot_sampleName.png",\
            "output/pca_mds/PCA_ggplot_ClassName.png",\
            "output/pca_mds/MDS_ggplot_sampleName.png",\
            "output/pca_mds/MDS_ggplot_ClassName.png",\
            "output/pca_mds/PCA_dim_genes_2D_samples.png",\
            "output/pca_mds/PCA_dim_genes_2D_classes.png",\
            "output/pca_mds/PCA_dim_genes_2D_overlay.png",\
            "output/pca_mds/PCA_Eigen_val.png",\
            "output/pca_mds/PCA_ggplot_sampleName.png",\
            "output/pca_mds/PCA_ggplot_ClassName.png",\
            "output/pca_mds/MDS_ggplot_sampleName.png",\
            "output/pca_mds/MDS_ggplot_ClassName.png"
    params: script=script, pheno=pheno_file.name 
    threads: 1
    message: "--- Performing PCA and MDS analysis."
    shell:  """
        {params.script} -i {input} -f output/pheno.txt -o output/pca_mds  -q -z &> {output[0]}.log
        touch {output}
        mv {params.pheno} output/pca_mds
    """