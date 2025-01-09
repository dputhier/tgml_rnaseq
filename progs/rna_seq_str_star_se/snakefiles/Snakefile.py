#================================================================#
# D.Puthier
# Last modification: Mon Jul 18 16:30:04 CEST 2016
#================================================================#



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

include: "rules/fastqc_raw_se_rule.py"
include: "rules/trimming_se_rule.py"
include: "rules/fastqc_trim_se_rule.py"
include: "rules/star_se.py"
include: "rules/select_read_by_strand_se_rule.py"
include: "rules/big_wig_rule.py"
include: "rules/indexing_bam_rule.py"
include: "rules/cufflink_rule.py"
include: "rules/cuffmerge_rule.py"
include: "rules/fCounts_known_gene_rule.py"
include: "rules/add_gene_name_to_selected_novel_tx_rule.py"
include: "rules/select_novel_transcript_rule.py"
include: "rules/merge_novel_and_know_tx_rule.py"
include: "rules/fCounts_known_and_novel_gene_rule.py"
include: "rules/mapping_stat_file.py"
include: "rules/mapping_stats_plot_rule.py"
include: "rules/rseqc_cov_rule.py"
include: "rules/pca_rule.py"
include: "rules/cor_plot_rule.py"
include: "rules/deseq2_rule.py"
include: "rules/cuffdiff.py"
"""



include: "rules/fasta_index_rule.py"

include: "rules/dag_rule.py"
#include: "rules/pair_plot_rule.py"




include: "rules/reportRmd.py"
"""

#================================================================#
#     Global variables                                           #
#================================================================#

workdir: config["workingdir"]

SAMPLES = config["samples"].split()

for curr_key in config:
    if "_INFO_" in curr_key:
        sys.stderr.write(config[curr_key] + "\n")

COMPARISON = config["comparison"].keys()
FASTA_INDEX = config["fasta"] + ".fai"
print(COMPARISON)

#================================================================#
#                         Workflow                               #
#================================================================#

FASTQC_RAW = expand("output/fastqc_raw/{smp}/{smp}_fastqc/fastqc_data.txt", smp=SAMPLES)
TRIMMING =  expand("output/trimmed/{smp}.fq.gz", smp=SAMPLES)
FASTQC_TRIM = expand("output/fastqc_trim/{smp}/{smp}_fastqc/fastqc_data.txt", smp=SAMPLES)
MAPPING = expand("output/bam/{smp}.bam", smp=SAMPLES)
BAM_BY_STRAND = expand("output/bam/{smp}_min.bam", smp=SAMPLES)
BAM_INDEX = expand("output/bam/{smp}.bam.bai", smp=SAMPLES)
BIGWIG = expand("output/bwig/{smp}_min.bw", smp=SAMPLES)
CUFFLINKS = expand("output/cufflinks/{smp}/transcripts.gtf", smp=SAMPLES)
CUFFMERGE = expand("output/cuffmerge/merged.gtf", smp=SAMPLES)
STAT_CUFFMERGE = "output/cuffmerge/cuffmerge_stats.png"
QUANTIF_KNOWN_GENES = expand("output/quantification_known_genes/gene_counts.txt", smp=SAMPLES)
NOVEL_SELECTED_TX = "output/cuffmerge/selected_novel_transcript_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"
GENE_NAME_FOR_NOVEL_TX = "output/cuffmerge/selected_novel_transcript_class_code_" + \
                          config["cuffmerge"]["selected_class_code"] + "_with_gene_name.gtf"
KNOWN_AND_NOVEL_TX = "output/inferred_gene_annotation/known_transcripts_and_selected_class_code_" + config["cuffmerge"]["selected_class_code"] + ".gtf"
QUANTIF_KNOWN_AND_NOVEL_GENES = ["output/quantification_known_and_novel_genes/gene_counts_known_and_novel.txt",\
                          "output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt"]
RSEQC_COV = "output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.curves.pdf.png"
MAPPING_STATS = expand("output/mapping_stats/{smp}.stats", smp=SAMPLES)
MAPPING_STAT_PLOT = expand("output/mapping_stats/{smp}.stats.png", smp=SAMPLES)
PCA = ["output/pca_mds/MDS_ggplot_ClassName.png"]
CORR_PLOT= [ "output/corr_plot/CoorPlot_circle.png", \
             "output/corr_plot/CoorPlot_ellipe.png", \
             "output/corr_plot/CoorPlot_pie.png",    \
             "output/corr_plot/CoorPlot_piewb.png",  \ 
             "output/corr_plot/CoorPlot_square.png", \
             "output/corr_plot/Pairs_plot.png"]

DESEQ2 = expand("output/comparison/{comp}/DESeq2_pval_and_norm_count_log2.txt", comp=COMPARISON)

DAG_PNG = "output/report/dag.png"

CUFFDIFF = "output/cuffdiff/RipmOVA_OTII.done"
"""


 

#DIAGNOSIS_PLOT = "output/diagnostic_plot/diagnostic.pdf"







REPORT = "output/report/report.html"


rule all:
    input: REPORT



, FASTQC_TRIM,        \  
            BAM_BY_STRAND, BIGWIG,          \
            QUANTIF_KNOWN_GENES,            \
            QUANTIF_KNOWN_AND_NOVEL_GENES,  \
            KNOWN_AND_NOVEL_TX,             \
            DAG_PNG,                        \
            RSEQC_COV,                      \
            MAPPING_STATS,               \
            MAPPING_STAT_PLOT,              \
            DESEQ2,                         \
            CORR_PLOT,                      \
            STAT_CUFFMERGE
"""

rule final:
    input:  FASTQC_RAW, \
            FASTQC_TRIM, \
            MAPPING, \
            BIGWIG, \
            CUFFLINKS, \
            CUFFMERGE, \
            QUANTIF_KNOWN_AND_NOVEL_GENES, \
            RSEQC_COV, \
            MAPPING_STATS, \
            MAPPING_STAT_PLOT, \
            CORR_PLOT, \
			DESEQ2 , \
			CUFFDIFF
    output: "output/code/Snakefile.py"
    params: wdir = config["workingdir"] + "/progs/rna_seq_str_star_se/snakefiles/*nake*", mem="2G"
    shell: """
    cp {params.wdir} {output}
    """

