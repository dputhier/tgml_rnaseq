#!/usr/bin/env python
"""
Description: Compute the number of reads per chromosome.
Version: {v}
"""
from __future__ import print_function
from collections import defaultdict
import sys
import argparse
import os
import pysam
import re
import pandas as pd
import seaborn as sns


__version__ = 0.1

def make_parser():

    parser = argparse.ArgumentParser()

    parser_grp_main = parser.add_argument_group('Arguments')


    parser_grp_main.add_argument('-i', '--bam',
                          help="Path to the BAM file. Default to STDIN",
                          default=sys.stdin,
                          metavar="BAM",
                          type=argparse.FileType("r"))

    parser_grp_main.add_argument('-o', '--outputfile',
                          help="output file",
                          default=sys.stdout,
                          metavar="GTF",
                          type=argparse.FileType('w'))

    parser_grp_main.add_argument("-q",
                                  "--mapq",
                                  type=int,
                                  help=" Only keep read with mapping quality equal or greater than q.",
                                  default=0,
                                  metavar="MAPQ",
                                  required=False)

    parser_grp_main.add_argument("-d",
                                  "--diagram-path",
                                  metavar="MAPQ",
                                  default=None,
                                  type=str,
                                  help="A png file name.",
                                  required=False)

    parser_grp_main.add_argument("-s",
                                  "--show-these-chr",
                                  type=str,
                                  default=None,
                                  help="Comma separated list of selected chromosomes to analyse. Data will be ordered accordingly.",
                                  required=False)

    return parser

def close_properly(*args):
    """Close a set of file if they are not None."""
    for afile in args:
        if afile is not None:
            if afile != sys.stdout:
                afile.close()

def write_properly(string, afile):
    """Write a string to a file. If file is None, write string to stdout."""
    if afile is not None:
        afile.write(string + "\n")
    else:
        sys.stdout.write(string + "\n")




def reads_per_chrom(bam=None,
                outputfile=None,
                mapq=None,
                mode=None,
                diagram_path=False,
                show_these_chr=None):

    chr_nb_read = defaultdict(int)


    if show_these_chr is not None:
        show_these_chr = show_these_chr.split(",")

    if not os.path.isfile(bam.name + ".bai"):
        bam_fn_test = re.sub(".bam", "", bam.name, re.IGNORECASE)
        if not os.path.isfile(bam_fn_test + ".bai"):
                sys.stderr.write("Can't find any bam index "
                                "for file : " + bam.name + "\n")
                sys.exit()

    samfile = pysam.Samfile(bam.name, "rb")

    sam_records = samfile.fetch()

    for read in sam_records:
        if read.mapq >= mapq:
            # read.rlen is the number of bases in the aligned read
            # excluding soft clipped bases.
            try:
                chr_nb_read[read.reference_name] += 1
            except:
                chr_nb_read[read.rname] += 1

    if show_these_chr is None:

        chr_nb_read_out = []
        chr_name_out = []

        for i in chr_nb_read.keys():
            out = i + "\t" + str(chr_nb_read[i])
            chr_nb_read_out += [chr_nb_read[i]]
            chr_name_out += [i]
            write_properly(out, outputfile)

    else:

        chr_nb_read_out = []
        chr_name_out = []

        for i in show_these_chr:
            try:
                out = i + "\t" + str(chr_nb_read[i])
                chr_nb_read_out += [chr_nb_read[i]]
                chr_name_out += [i]
                write_properly(out, outputfile)
            except:
                pass

    if diagram_path is not None:
        sns_plot = sns.barplot(x=chr_name_out,
                               y=chr_nb_read_out,
                               palette="Greys_d")

        sns_plot.figure.savefig(diagram_path)
    close_properly(outputfile, bam)

if __name__ == '__main__':
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    reads_per_chrom(**args)
else:
    cmd = CmdObject(name="bam_sum",
                    message="Compute the number of read per chromosome.",
                    parser=make_parser(),
                    fun=reads_per_chrom,
                    desc=__doc__.format(v=__version__))



