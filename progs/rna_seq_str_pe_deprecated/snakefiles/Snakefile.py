#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
import sys
import re
from snakemake.utils import report
from snakemake.utils import R

configfile: "snakefiles/config.json"

#================================================================#
#                         Includes                               #
#================================================================#

include: "rules/fastqc_raw_pe_rule.py"
include: "rules/trimming_pe_rule.py"
include: "rules/fastqc_trim_pe_rule.py"
include: "rules/tophat_pe_rule.py"
include: "rules/select_read_by_strand_pe_rule.py"
include: "rules/big_wig_rule.py"
include: "rules/indexing_bam_rule.py"
include: "rules/fCounts_known_gene_rule.py"
include: "rules/cufflink_rule.py"
include: "rules/cuffmerge_rule.py"
include: "rules/fasta_index_rule.py"
include: "rules/add_gene_name_to_selected_novel_tx_rule.py"
include: "rules/select_novel_transcript_rule.py"
include: "rules/merge_novel_and_know_tx_rule.py"
include: "rules/fCounts_known_and_novel_gene_rule.py"
include: "rules/rseqc_cov_rule.py"
include: "rules/dag_rule.py"
#include: "rules/pair_plot_rule.py"
include: "rules/mapping_stat_file_R1_rule.py"
include: "rules/mapping_stat_file_R2_rule.py"
include: "rules/mapping_stats_plot_rule.py"
include: "rules/deseq2_rule.py"
include: "rules/pca_rule.py"
include: "rules/cor_plot_rule.py"
include: "rules/stat_cuffmerge_rule.py"
include: "rules/prep_stat_cuffmerge_rule.py"

#================================================================#
#     Global variables                                           #
#================================================================#

workdir: config["workingdir"]

SAMPLES = config["samples"].split()

#================================================================#
#                         Workflow                               #
#================================================================#

FASTQC_RAW = expand("output/fastqc_raw/{smp}/{smp}_R1.fq_fastqc/fastqc_data.txt", smp=SAMPLES)

TRIMMING =  expand("output/trimmed/{smp}_R1_t.fq.gz", smp=SAMPLES)

FASTQC_TRIM = expand("output/fastqc_trim/{smp}/{smp}_R1_t.fq_fastqc/fastqc_data.txt", smp=SAMPLES)

MAPPING = expand("output/bam/{smp}.bam", smp=SAMPLES)

BAM_BY_STRAND = expand("output/bam/{smp}_min.bam", smp=SAMPLES)

BAM_INDEX = expand("output/bam/{smp}.bam.bai", smp=SAMPLES)

BIGWIG = expand("output/bwig/{smp}_min.bw", smp=SAMPLES)

QUANTIF_KNOWN_GENES = expand("output/quantification_known_genes/gene_counts.txt", smp=SAMPLES)

CUFFLINKS = expand("output/cufflinks/{smp}/transcripts.gtf", smp=SAMPLES)

FASTA_INDEX = config["fasta"] + ".fai"

CUFFMERGE = expand("output/cuffmerge/merged.gtf", smp=SAMPLES)

STAT_CUFFMERGE = "output/cuffmerge/cuffmerge_stats.png"

NOVEL_SELECTED_TX = "output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"

GENE_NAME_FOR_NOVEL_TX = "output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + "_with_gene_name.gtf"

KNOWN_AND_NOVEL_TX = "output/inferred_gene_annotation/known_transcripts_and_selected_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"

QUANTIF_KNOWN_AND_NOVEL_GENES = ["output/quantification_known_and_novel_genes/gene_counts_known_and_novel.txt", \
                                "output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt"]

RSEQC_COV = "output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.curves.pdf.png"
 
DAG_PNG = "output/report/dag.png"

#DIAGNOSIS_PLOT = "output/diagnostic_plot/diagnostic.pdf"

MAPPING_STATS_R1 = expand("output/mapping_stats/{smp}_R1.stats", smp=SAMPLES)

MAPPING_STATS_R2 = expand("output/mapping_stats/{smp}_R2.stats", smp=SAMPLES)

MAPPING_STAT_PLOT = expand("output/mapping_stats/{smp}.stats.png", smp=SAMPLES)

DESEQ2 = [  "output/diff_call_deseq2/DESeq2_diagnostic_MA.png",         \
            "output/diff_call_deseq2/DESeq2_diagnostic_disp.png",       \ 
            "output/diff_call_deseq2/DESeq2_norm_count_table_lin.txt",  \
            "output/diff_call_deseq2/DESeq2_norm_count_table_log2.txt", \
            "output/diff_call_deseq2/DESeq2_norm_count_table_rld.txt",  \
            "output/diff_call_deseq2/DESeq2_norm_count_table_vsd.txt",  \ 
            "output/diff_call_deseq2/DESeq2_raw_count_table.txt"]
DESEQ2_F = "output/diff_call_deseq2/DESeq2_diff_genes.txt"

PCA_MDS =   ["output/pca_mds/PCA_dim_genes_2D_samples.png",\
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
            "output/pca_mds/MDS_ggplot_ClassName.png"]

CORR_PLOT= [ "output/corr_plot/CoorPlot_circle.png", \
            "output/corr_plot/CoorPlot_ellipe.png", \
            "output/corr_plot/CoorPlot_pie.png",    \
            "output/corr_plot/CoorPlot_piewb.png",  \ 
            "output/corr_plot/CoorPlot_square.png", \
            "output/corr_plot/Pairs_plot.png"]

ruleorder:

rule all:
    input: "output/report/report.html"

rule final:
    input:  FASTQC_RAW, FASTQC_TRIM,        \  
            BAM_BY_STRAND, BIGWIG,          \
            QUANTIF_KNOWN_GENES,            \
            QUANTIF_KNOWN_AND_NOVEL_GENES,  \
            KNOWN_AND_NOVEL_TX,             \
            DAG_PNG,                        \
            RSEQC_COV,                      \
#            DIAGNOSIS_PLOT,                 \
            MAPPING_STATS_R1,               \
            MAPPING_STATS_R2,               \
            MAPPING_STAT_PLOT,              \
            DESEQ2,                         \
            PCA_MDS,                        \
            CORR_PLOT,                      \
            STAT_CUFFMERGE

    output: "output/code/Snakefile.py"
    params: wdir = config["workingdir"] + "/progs/rna_seq_str_pe/snakefiles/*nake*"
    shell: """
    cp {params.wdir} {output}
    """

#================================================================#
#            functions for Report                                #
#================================================================#

def report_numbered_list(list):
    result = ""
    n = 0
    for element in list:
        n = n + 1
        result = result + "   " + str(n) + ". " + element + "\n"
        
    return(result)


def report_link_list(list):
    if isinstance(list, str):
        list = [list]
    result = ""
    n = 0
    for element in list:
        n+=1
        result += "- `" + element + " <../../" + element + ">`_ \n\n" 
    
    return(result + "\n") 

def report_bullet_list(alist):

    return "\n".join(["- "+ x + "\n" for x in  alist])

#================================================================#
#           Table of images                                      #
#================================================================#


def image_fastq(alist, prefix="a_"):

    if isinstance(alist, str):
        alist = [alist]
        
    table ="""
+---------+-------------------+-------------------+-------------------+      
+ Sample  +Per base Quality   +Duplication levels +      K-mers       +
+---------+-------------------+-------------------+-------------------+"""

    row = """   
+    {s}  + |{p}pbq{n}|     +   |{p}dup{n}|   +     |{p}kmr{n}| +
+---------+-------------------+-------------------+-------------------+"""

    alist = [x.replace("output/","") for x in alist]

    pbq = [x.replace("fastqc_data.txt", "Images/per_base_quality.png") for x in alist]
    dup = [x.replace("fastqc_data.txt", "Images/duplication_levels.png") for x in alist]
    kmr = [x.replace("fastqc_data.txt", "Images/kmer_profiles.png") for x in alist]
      
    pbq = "\n".join([" .. |" + prefix + "pbq"  + str(p).zfill(3) + "| image:: ../" + x + "\n" for p,x in  enumerate(pbq)])
    dup = "\n".join([" .. |" + prefix + "dup"  + str(p).zfill(3) + "| image:: ../" + x + "\n" for p,x in  enumerate(dup)])
    kmr = "\n".join([" .. |" + prefix + "kmr"  + str(p).zfill(3) + "| image:: ../" + x + "\n" for p,x in  enumerate(kmr)])
     
    for i in range(len(alist)):
        table += row.format(n=str(i).zfill(3), p=prefix, s=str(i + 1).zfill(3),)
    
    result = "\n\n" + pbq + "\n" + dup + "\n" + kmr + "\n" + table + "\n\n"

    return result

def image_other(alist, name="mapstat", addheader=True):

    if isinstance(alist, str):
        alist = [alist]

    alist = [x.replace("output/","") for x in alist]

    header="""
+---------+----------------------+      
+ Sample  + Mapping statistics   +"""

    table ="""
+---------+----------------------+"""

    row = """   
+    {s}  + |{name}{n}|         +
+---------+----------------------+"""

    if addheader:
        table = header + table

    mapstat = "\n".join([" .. |" + name  + str(p).zfill(3) + "| image:: ../" + x + "\n" for p,x in  enumerate(alist)])

    for i in range(len(alist)):
        table += row.format(n=str(i).zfill(3), s=str(i + 1).zfill(3),name=name)

    result = "\n\n" + mapstat + "\n" + table + "\n\n"

    return result
    
#================================================================#
#           Variables  for report                                #
#================================================================#

## info
SAMPLES_L = report_numbered_list(SAMPLES)

## dag
DAG_PNG_L = os.path.basename(DAG_PNG)

SINGLE_DAG_PNG_L = os.path.basename(DAG_PNG.replace("dag","rulegraph"))
 
## BAM

BAM_L = report_link_list(MAPPING)
BWIG_L = report_link_list(BIGWIG)

## FASTQC
# link
FASTQC_RAW_R1_L = report_link_list(FASTQC_RAW)
FASTQC_RAW_R2_L = FASTQC_RAW_R1_L.replace("R1.fq_fastqc","R2.fq_fastqc")
# images
FASTQC_TRIM_R1_L = report_link_list(FASTQC_TRIM)
FASTQC_TRIM_R2_L = FASTQC_TRIM_R1_L.replace("R1_t.fq_fastqc","R2_t.fq_fastqc")

FASTQC_RAW_SEQ_R1_QUAL_IT = image_fastq(FASTQC_RAW, prefix="a____")
FASTQC_RAW_SEQ_R2_QUAL_IT = FASTQC_RAW_SEQ_R1_QUAL_IT.replace("R1.fq_fastqc","R2.fq_fastqc")
FASTQC_RAW_SEQ_R2_QUAL_IT = FASTQC_RAW_SEQ_R2_QUAL_IT.replace("a____","b____")

FASTQC_TRIM_SEQ_R1_QUAL_IT = image_fastq(FASTQC_TRIM, prefix="c____")
FASTQC_TRIM_SEQ_R2_QUAL_IT = FASTQC_TRIM_SEQ_R1_QUAL_IT.replace("R1_t.fq_fastqc","R2_t.fq_fastqc")
FASTQC_TRIM_SEQ_R2_QUAL_IT = FASTQC_TRIM_SEQ_R2_QUAL_IT.replace("c____","d____")

## Mapping statistics
MAPPING_STAT_PLOT_I = image_other(MAPPING_STAT_PLOT)
MAPPING_STAT_PLOT_L = report_link_list([x.replace(".png","") for x in MAPPING_STAT_PLOT])

## RSEQC
RSEQC_COV_L = report_link_list(RSEQC_COV.replace(".curves.pdf.png",".txt"))
RSEQC_COV_I =  ".. image:: " + RSEQC_COV.replace("output/","../") + "\n\n"


## DESEQ2
DESEQ2_I = image_other([x for x in DESEQ2 if x.endswith("png")], name="deseq2_", addheader=False)
DESEQ2_F = report_link_list(DESEQ2_F)

## CORR_PLOT
CORR_PLOT_I = image_other(CORR_PLOT, name="corrplo", addheader=False)

## PCA and MDS
PCA_MDS_I = image_other(PCA_MDS, name="pcaplot", addheader=False)

## Statistics cuffmerge
STAT_CUFFMERGE_I = image_other(STAT_CUFFMERGE, name="cuffmer", addheader=False)

#================================================================#
#           Report                                               #
#================================================================#



rule report:
    """
    Generate a report with the list of datasets + summary of the results.
    """
    input:  code="output/code/Snakefile.py"
            

    params: wdir=config["workingdir"], \
            user=config["user"], \
            dag_png=os.path.basename(DAG_PNG), \
            nb_smp=str(len(config["samples"].split()))

    output: html="output/report/report.html"

    run:
        report("""
        =========================
        RNA-Seq analysis summary
        =========================
        
        About
        ======
        
        - User name : {params.user}
        - Directory path : {params.wdir}

        Contents
        ========
        
        - `Flowcharts`_
        - `Datasets description`_
        - `FastQC raw reads (R1)`_
        - `FastQC trimmed reads (R1)`_
        - `FastQC raw reads (R2)`_
        - `FastQC trimmed reads (R2)`_
        - `RSeQC genebody coverage`_
        - `Mapping statistics`_
        - `BAM files`_
        - `Bigwig files`_

    
        -----------------------------------------------------

        Flowcharts
        ==========

        - Layout

        .. image:: {SINGLE_DAG_PNG_L}
        
        - Sample-wise worflow

        .. image:: {DAG_PNG_L}


        -----------------------------------------------------

        
        Datasets description
        =====================
        
        - Number of samples: {params.nb_smp}

        - Sample names

        {SAMPLES_L}


        -----------------------------------------------------
        
        FastQC raw reads (R1)
        ======================
        
        {FASTQC_RAW_SEQ_R1_QUAL_IT}
           
        {FASTQC_RAW_R1_L}


        FastQC trimmed reads (R1)
        ==========================
        
        {FASTQC_TRIM_SEQ_R1_QUAL_IT}
             
        {FASTQC_TRIM_R1_L}
        

        
        FastQC raw reads (R2)
        ======================

        {FASTQC_RAW_SEQ_R2_QUAL_IT}
             
        {FASTQC_RAW_R2_L}


        FastQC trimmed reads (R2)
        ==========================
                
        {FASTQC_TRIM_SEQ_R2_QUAL_IT}
           
        {FASTQC_TRIM_R2_L}


        -----------------------------------------------------

        RSeQC genebody coverage
        =========================
        
        {RSEQC_COV_I}
        
        {RSEQC_COV_L}
        
        -----------------------------------------------------

        Mapping statistics
        =========================
        
        {MAPPING_STAT_PLOT_I}
        
        {MAPPING_STAT_PLOT_L}
        
        -----------------------------------------------------
                
        BAM files
        ==============
        
        {BAM_L}
        -----------------------------------------------------
                
        Bigwig files
        ==============
        {BWIG_L}
        -----------------------------------------------------

        Correlation plot
        ================

        {CORR_PLOT_I}

        -----------------------------------------------------

        Cuffmerge statistics
        =====================
        
        {STAT_CUFFMERGE_I}

        DESeq2 output
        ==============

        {DESEQ2_I}
        {DESEQ2_F}

        PCA and MDS: various representations
        ====================================
        
        {PCA_MDS_I}
        
        """, output.html, metadata="D. Puthier", **input)
