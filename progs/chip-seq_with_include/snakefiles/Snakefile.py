#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
configfile: "snakefiles/config.json"

#================================================================#
#                         Includes                               #
#================================================================#

include: "rules/fastqc_rules.py"
include: "rules/trimming.rules"
include: "rules/mapping.rules"
include: "rules/sorting_bam.rules"
include: "rules/indexing.rules"
include: "rules/flagstat.rules"
include: "rules/select_mapped.rules"
include: "rules/delete_replicats.rules"
include: "rules/rulegraph.rules"
include: "rules/dag.rules"

#================================================================#
#     Global variables                                           #
#================================================================#

workdir: config["workingdir"]

SAMPLES = config["samples"].split()
SAMPLES_TYPE = config["samples_type"].split() 

#================================================================#
#     Data should be downloaded from SRA (or not)                                           #
#================================================================#

if config["sra"] ==  "yes":

#================================================================#
#                         Workflow                               #
#================================================================#

FASTQC = expand("results/fastqc/{samples_type}/{samples}.fq_fastqc.zip", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES)

TRIMMING =  expand("results/trimmed/{samples_type}/{samples}_trimmed.fq.gz", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES)

MAPPING = expand("results/mapped/{samples_type}/{samples}_mapped.bam", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES)

MAPPING_Q30 =  expand("results/mapped/{samples_type}/{samples}_mapped_q30.bam", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES)

MAPPING_UNIQUE = expand("results/mapped/{samples_type}/{samples}_mapped_unique.bam", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES)

STAT = expand("results/stat/{samples_type}/{samples}_{stat_type}.txt", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES, stat_type="mapped mapped".split())

STAT_Q30 = expand("results/stat/{samples_type}/{samples}_{stat_type}.txt", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES, stat_type="mapped_q30 mapped_q30".split())

STAT_UNIQUE = expand("results/stat/{samples_type}/{samples}_{stat_type}.txt", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES, stat_type="mapped_unique mapped_unique".split())

INDEX = expand("results/mapped/{samples_type}/{samples}_mapped_sorted.bam.bai", zip, samples_type=SAMPLES_TYPE, samples=SAMPLES)


rule all:
    input: FASTQC , \
      STAT, \
      STAT_Q30, \
      STAT_UNIQUE,\
      INDEX,\
      "reports/rulegraph.png",\
      "reports/dag.png",\
      "reports/report.html"
    



#================================================================#
#                          Report                                #
#================================================================#

from snakemake.utils import report

WD = config["workingdir"]

def report_numbered_list(list):
    result = ""
    n = 0
    for element in list:
        n = n + 1
        result = result + str(n) + ". " + element + "\n"
        
    return(result)


def report_link_list(list):
    result = ""
    n = 0
    for element in list:
        n+=1
        result += "`" + element + " <../" + element + ">`_ \n" 

    return(result) 

SAMPLES_TYPE_L = report_numbered_list(SAMPLES_TYPE)
SAMPLES_L = report_numbered_list(SAMPLES)
TRIMMING_L = report_numbered_list(TRIMMING)
INDEX_L = report_numbered_list(INDEX)
STAT_L = report_link_list(STAT)
STAT_Q30_L = report_link_list(STAT_Q30)
STAT_UNIQUE_L = report_link_list(STAT_UNIQUE)

rule report:
    """
    Generate a report with the list of datasets + summary of the results.
    """
    input:  dag_png="reports/dag.png", \
            rulegraph_png="reports/rulegraph.png"
    output: html="reports/report.html"
    run:
        report("""
        ==================
        Chip-seq analysis
        ==================
        
        :Analysis workflow:    Justine Long, Jeanne Chèneby
        
        Contents
        ========
        
        - `Flowcharts`_
        - `Datasets`_
             - `Sample types`_
             - `Sample names`_

        -----------------------------------------------------

        Flowcharts
        ==========

        - Sample treatment: dag_png_
        - Workflow: rulegraph_png_

        .. image:: rulegraph.png

        -----------------------------------------------------

        Datasets
        ========
        
        Sample types
        -------------------
        {SAMPLES_TYPE_L}

        Sample names
        ------------------
        {SAMPLES_L}

        Result files
        ============

        Quality control (raw reads)
        ---------------------------

        .. |logo| image:: ../results/fastqc/chip_seq/siNT_ER_E2_r3_SRX176860_chr21_0.6_Noise.fq_fastqc/Images/per_base_quality.png



        TRIMMING
        -------------
        {TRIMMING_L}

        INDEX
        --------
        Files to visualize with a genome browser \n

        {INDEX_L}

        STAT
        --------
        {STAT_L}

        STAT_Q30
        ----------
        {STAT_Q30_L}

        STAT_UNIQUE
        -------------
        {STAT_UNIQUE_L}



        """, output.html, metadata="Jeanne Chèneby and Justine Long", **input)
