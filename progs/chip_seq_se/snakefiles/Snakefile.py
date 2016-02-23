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

include: "rules/fastqc_rules.py"
include: "rules/trimming_se_rules.py"
include: "rules/fastqc_trim_se_rules.py"
include: "rules/bowtie2_rules.py"
include: "rules/mapping_stat_chipSeq_rule.py"
include: "rules/delete_duplicate_read_rules.py"
include: "rules/mapping_quality_rules.py"
include: "rules/samtools_sort_rules.py"
include: "rules/dag_rule.py"
include: "rules/mapping_stats_plot_chip_seq_rule.py"

#================================================================#
#     Global variables                                           #
#================================================================#

workdir: config["workingdir"]

SAMPLES = config["samples"].split()


#================================================================#
#                         TARGETS                                #
#================================================================#

FASTQC_RAW = expand("output/fastqc_raw/{smp}/{smp}.fq_fastqc/fastqc_data.txt", smp=SAMPLES)

TRIMMING =  expand("output/trimmed/{smp}_t.fq.gz", smp=SAMPLES)

FASTQC_TRIM = expand("output/fastqc_trim/{smp}/{smp}_t.fq_fastqc/fastqc_data.txt", smp=SAMPLES)

MAPPING = expand("output/bam/{smp}.bam", smp=SAMPLES)

MAPPING_SORT = expand("output/bam/{smp}_srt.bam", smp=SAMPLES)

MAPPING_QUAL =  expand("output/bam/{smp}_q30.bam", smp=SAMPLES)

MAPPING_DUP = expand("output/bam/{smp}_q30_rmDup.bam", smp=SAMPLES)

MAPPING_STATS=expand("output/mapping_stats/{smp}q30_rmDup.flagstat", smp=SAMPLES)

MAPPING_STAT_PLOT = expand( "output/mapping_stats/{smp}.stats.png", smp=SAMPLES)


DAG_PNG = "output/report/dag.png"

#================================================================#
#                         LAST RULES                             #
#================================================================#


rule all:
    input: "output/report/report.html"

rule final:
    input:  FASTQC_RAW, FASTQC_TRIM, MAPPING_STATS, MAPPING_STAT_PLOT, DAG_PNG
    output: "output/code/Snakefile.py"
    params: wdir = config["workingdir"] + "progs/chip_seq_se/snakefiles/Snakefile.py"
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

def image_mapp_stats(alist):

    if isinstance(alist, str):
        alist = [alist]

    alist = [x.replace("output/","") for x in alist]
        
    table ="""
+---------+----------------------+      
+ Sample  + Mapping statistics   +
+---------+----------------------+"""

    row = """   
+    {s}  + |mapstat{n}|         +
+---------+----------------------+"""

    mapstat = "\n".join([" .. |" + "mapstat"  + str(p).zfill(3) + "| image:: ../" + x + "\n" for p,x in  enumerate(alist)])

    for i in range(len(alist)):
        table += row.format(n=str(i).zfill(3), s=str(i + 1).zfill(3),)

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

BAM_L = report_link_list(MAPPING_DUP)


## FASTQC
# link
FASTQC_RAW_L = report_link_list(FASTQC_RAW)
FASTQC_TRIM_L = report_link_list(FASTQC_TRIM)

# images

FASTQC_RAW_SEQ_QUAL_IT = image_fastq(FASTQC_RAW, prefix="a____")

FASTQC_TRIM_SEQ_QUAL_IT = image_fastq(FASTQC_TRIM, prefix="c____")

## Mapping statistics
MAPPING_STAT_PLOT_I = image_mapp_stats(MAPPING_STAT_PLOT)
MAPPING_STAT_PLOT_L = report_link_list([x.replace(".png","") for x in MAPPING_STAT_PLOT])


#================================================================#
#                         REPORT                                 #
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
        ChIP-Seq analysis summary
        =========================
        
        About
        ======
        
        - User name : {params.user}
        - Directory path : {params.wdir}

        Contents
        ========
        
        - `Flowcharts`_
        - `Datasets description`_
        - `FastQC raw reads`_
        - `FastQC trimmed reads`_
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
        
        FastQC raw reads
        ======================
        
        {FASTQC_RAW_SEQ_QUAL_IT}
           
        {FASTQC_RAW_L}


        FastQC trimmed reads (R2)
        ==========================
                
        {FASTQC_TRIM_SEQ_QUAL_IT}
           
        {FASTQC_TRIM_L}

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


        """, output.html, metadata="D. Puthier", **input)
        