#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
from snakemake.utils import report

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

#================================================================#
#     Global variables                                           #
#================================================================#

workdir: config["workingdir"]

SAMPLES = config["samples"].split()
SAMPLES_TYPE = config["samples_type"].split()


#================================================================#
#                         Workflow                               #
#================================================================#

FASTQC_RAW = expand("output/fastqc_raw/{smp}/{smp}_R1.fq_fastqc.html", zip,  smp=SAMPLES)

TRIMMING =  expand("output/trimmed/{smp}_R1_t.fq.gz", zip,  smp=SAMPLES)

FASTQC_TRIM = expand("output/fastqc_trim/{smp}/{smp}_R1_t.fq_fastqc.html", zip,  smp=SAMPLES)

MAPPING = expand("output/bam/{smp}.bam", zip,  smp=SAMPLES)

BAM_BY_STRAND = expand("output/bam/{smp}_min.bam", zip,  smp=SAMPLES)

BAM_INDEX = expand("output/bam/{smp}.bam.bai", zip,  smp=SAMPLES)

BIGWIG = expand("output/bwig/{smp}_min.bw", zip,  smp=SAMPLES)

QUANTIF_KNOWN_GENES = expand("output/quantification_known_genes/gene_counts.txt", zip,  smp=SAMPLES)

CUFFLINKS = expand("output/cufflinks/{smp}/transcripts.gtf", zip,  smp=SAMPLES)

FASTA_INDEX = config["fasta"] + ".fai"

CUFFMERGE = expand("output/cuffmerge/merged.gtf", zip,  smp=SAMPLES)

NOVEL_SELECTED_TX = "output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"

GENE_NAME_FOR_NOVEL_TX = "output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + "_with_gene_name.gtf"

KNOWN_AND_NOVEL_TX = "output/inferred_gene_annotation/known_transcripts_and_selected_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"

QUANTIF_KNOWN_AND_NOVEL_GENES = "output/quantification_known_and_novel_genes/gene_counts_known_and_novel.txt"

RSEQC_COV = "output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.txt"
 
DAG_PNG = "output/report/dag.png"

MAPPING_Q30 =  expand("results/mapped/{smp}/{smp}_mapped_q30.bam", zip,  smp=SAMPLES)

MAPPING_UNIQUE = expand("results/mapped/{smp}/{smp}_mapped_unique.bam", zip,  smp=SAMPLES)

STAT = expand("results/stat/{smp}/{smp}_{stat_type}.txt", zip,  smp=SAMPLES, stat_type="mapped mapped".split())

STAT_Q30 = expand("results/stat/{smp}/{smp}_{stat_type}.txt", zip,  smp=SAMPLES, stat_type="mapped_q30 mapped_q30".split())

STAT_UNIQUE = expand("results/stat/{smp}/{smp}_{stat_type}.txt", zip,  smp=SAMPLES, stat_type="mapped_unique mapped_unique".split())

INDEX = expand("results/mapped/{smp}/{smp}_mapped_sorted.bam.bai", zip,  smp=SAMPLES)


rule final:
    input:  FASTQC_RAW,     \
            FASTQC_TRIM,    \   
            MAPPING,        \
            BAM_BY_STRAND,  \
            QUANTIF_KNOWN_GENES, \
            QUANTIF_KNOWN_AND_NOVEL_GENES,\
            DAG_PNG,        \
            BIGWIG,         \
            RSEQC_COV,      \
            "output/report/report.html"

#================================================================#
#            functions for Report                                #
#================================================================#

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
        result += "    - `" + element + " <../../" + element + ">`_ \n\n" 

    return(result) 

#================================================================#
#           Variables  for report                                #
#================================================================#


SAMPLES_L = report_numbered_list(SAMPLES)
TRIMMING_L = report_numbered_list(TRIMMING)

FASTQC_RAW_L = report_link_list(FASTQC_RAW)
FASTQC_TRIM_L = report_link_list(FASTQC_TRIM)


#================================================================#
#           Report                                               #
#================================================================#


rule report:
    """
    Generate a report with the list of datasets + summary of the results.
    """
    input:  dag_png=DAG_PNG, \
            rseqc=RSEQC_COV
            

    params: wdir=config["workingdir"], \
            user=config["user"], \
            dag_png=os.path.basename(DAG_PNG), \
            nb_smp=str(len(config["samples"].split())),\
            bwig_path="\n".join(["- "+ x + "\n" for x in  BIGWIG]), \
            bam_path="\n".join(["- "+ x + "\n" for x in  MAPPING])

    output: html="output/report/report.html"

    run:
        report("""
        =========================
        RNA-Seq analysis summary
        =========================
        
        :Analysis workflow:
        
        - User name : {params.user}
        - Directory path : {params.wdir}

        Contents
        ========
        
        - `Flowcharts`_
        - `Datasets description`_
            - `Number of samples`_
            - `Sample names`_
        - `FastQC reports`_
            - `FastQC raw reads`_
            - `FastQC trimmed reads`_
        - `RSeQC genebody coverage`_
        - `BAM files`_
        - `Bigwig files`_
        

        -----------------------------------------------------

        Flowcharts
        ==========

        - Sample processing file: {params.dag_png}

        .. image:: ../../{input.rseqc}

        -----------------------------------------------------

        
        Datasets description
        =====================
        
        Number of samples
        ------------------
 
        - {params.nb_smp}

        Sample names
        -------------

        {SAMPLES_L}

        -----------------------------------------------------
        
        FastQC reports
        ===============
        
        FastQC raw reads
        ------------------
        
        {FASTQC_RAW_L}
        
        FastQC trimmed reads
        ---------------------
        
        {FASTQC_TRIM_L}        

        
        -----------------------------------------------------

        RSeQC genebody coverage
        =========================
        
        .. image:: {params.dag_png}
        
        -----------------------------------------------------
                
        BAM files
        ==============
        
        {params.bam_path}
        -----------------------------------------------------
                
        Bigwig files
        ==============
        
        {params.bwig_path}

        """, output.html, metadata="D. Puthier", **input)
