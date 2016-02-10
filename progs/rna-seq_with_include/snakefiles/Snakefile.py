#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
configfile: "snakefiles/config.json"

#================================================================#
#                         Includes                               #
#================================================================#

include: "rules/fastqc_raw_pe_rule.py"
include: "rules/trimming_pe_rule.py"
include: "rules/fastqc_trim_pe_rule.py"
include: "rules/tophat_pe_rule.py"
include: "rules/select_read_by_strand_pe_rule.py"
include: "rules/indexing_bam_rule.py"
include: "rules/fCounts_known_gene_rule.py"
include: "rules/cufflink_rule.py"

#================================================================#
#     Global variables                                           #
#================================================================#

workdir: config["workingdir"]

SAMPLES = config["samples"].split()
SAMPLES_TYPE = config["samples_type"].split()


#================================================================#
#                         Workflow                               #
#================================================================#

FASTQC_RAW = expand("output/fastqc_raw/{smp}/{smp}_R1.fq_fastqc.zip", zip,  smp=SAMPLES)

TRIMMING =  expand("output/trimmed/{smp}_R1_t.fq.gz", zip,  smp=SAMPLES)

FASTQC_TRIM = expand("output/fastqc_trim/{smp}/{smp}_R1_t.fq_fastqc.zip", zip,  smp=SAMPLES)

MAPPING = expand("output/bam/{smp}.bam", zip,  smp=SAMPLES)

BAM_BY_STRAND= expand("output/bam/{smp}_min.bam", zip,  smp=SAMPLES)

BAM_INDEX= expand("output/bam/{smp}.bam.bai", zip,  smp=SAMPLES)

QUANTIF_KNOWN_GENES= expand("output/quantification_known_genes/gene_counts.txt", zip,  smp=SAMPLES)

CUFFLINKS = expand("output/cufflinks/{smp}/transcripts.gtf", zip,  smp=SAMPLES)



MAPPING_Q30 =  expand("results/mapped/{smp}/{smp}_mapped_q30.bam", zip,  smp=SAMPLES)

MAPPING_UNIQUE = expand("results/mapped/{smp}/{smp}_mapped_unique.bam", zip,  smp=SAMPLES)

STAT = expand("results/stat/{smp}/{smp}_{stat_type}.txt", zip,  smp=SAMPLES, stat_type="mapped mapped".split())

STAT_Q30 = expand("results/stat/{smp}/{smp}_{stat_type}.txt", zip,  smp=SAMPLES, stat_type="mapped_q30 mapped_q30".split())

STAT_UNIQUE = expand("results/stat/{smp}/{smp}_{stat_type}.txt", zip,  smp=SAMPLES, stat_type="mapped_unique mapped_unique".split())

INDEX = expand("results/mapped/{smp}/{smp}_mapped_sorted.bam.bai", zip,  smp=SAMPLES)


rule final:
    input:  FASTQC_RAW,     \
            TRIMMING,       \
            FASTQC_TRIM,    \
            MAPPING,        \
            BAM_BY_STRAND,  \
            BAM_INDEX,      \
            QUANTIF_KNOWN_GENES, \
            CUFFLINKS
            
    



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
