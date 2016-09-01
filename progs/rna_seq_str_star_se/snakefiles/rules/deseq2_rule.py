import os
script = os.path.join(config["workingdir"], "progs", "scripts", "twoClasseComparison4Counts.R")

def getclass1(wildcards):

    comp = wildcards["comp"]
    pheno = config["comparison"][comp]
    pheno = pheno.split()

    return(pheno[0])

def getclass2(wildcards):

    comp = wildcards["comp"]
    pheno = config["comparison"][comp]
    pheno = pheno.split()

    return(pheno[1])

rule deseq2:
    input:  "output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt"
    output: "output/comparison/{comp}/DESeq2_pval_and_norm_count_log2.txt"
    params: class1=getclass1, class2=getclass2, script=script
    threads: 1
    message: "-- Performing DESeq analysis --"
    shell:  """
        {params.script} -i {input} -c {params.class1} -d  {params.class2} -o output/comparison/{wildcards.comp}/ &> {output[0]}.log
    """
