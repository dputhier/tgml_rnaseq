#!/usr/bin/env python
"""
Description: Compute the coverage sum. Eg. If 1000 reads of size 10 are encountered this will be 1000*10.
Version: {v}
"""
from __future__ import print_function
import sys
import argparse
import os
import pysam
import re

__version__ = 0.1

def make_parser():
    parser = argparse.ArgumentParser(add_help=True)

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
                                  metavar="\b",
                                  required=False)
    parser_grp_main.add_argument("-m",
                                  "--mode",
                                  choices=["se", "pe"],
                                  help=" Only keep read with mapping quality equal or greater than q.",
                                  default=0,
                                  metavar="\b",
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

def compute_cov(bam=None,
                outputfile=None,
                mapq=None,
                mode=None):

    cur_len = 0

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
            cur_len += abs(int(read.qlen))

    write_properly(str(cur_len), outputfile)

    close_properly(outputfile, bam)

if __name__ == '__main__':
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    compute_cov(**args)



