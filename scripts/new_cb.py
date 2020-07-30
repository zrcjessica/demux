#!/usr/bin/env python

import csv
import pysam
import argparse
from pathlib import Path


def parse_format(fmt, f):
    """ what input or output format should we use? """
    if fmt is None:
        fmt = 'b' if Path(f).suffix == '.bam' else 's'
    return fmt

def read_barcodes(tsv):
    """ retrieve the tsv file containing the barcodes """
    with open(tsv) as infile:
        print(infile)
        barcodes = {
            rows[0]: rows[1]
            for rows in csv.reader(infile, delimiter="\t")
        }
    return barcodes

def main(barcodes, reads, in_format=None):
    """
        read the reads and alter their barcodes, returning each line as a
        generator function; the first line will be the pysam object
    """
    if in_format != "b":
        in_format = ''
    reads = pysam.AlignmentFile(reads, "r"+in_format)
    yield reads
    for read in reads:
        yield read

def write_reads(out, reads, out_format=None):
    """
        read each read from the iterator 'reads' and write it to out
    """
    if out_format != 'b':
        out_format = ''
    out = pysam.AlignmentFile(out, "w"+out_format, template=next(reads))
    for read in reads:
        out.write(read)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert CB tags in a SAM or BAM file")
    parser.add_argument(
        "-o", "--out", default="-", help="the filename to save the sam (or bam) to; defaults to stdout"
    )
    parser.add_argument(
        "barcodes", type=Path, help="a TSV file containing two columns: 1) barcodes that should be changed, and 2) the new barcodes they should be changed to"
    )
    parser.add_argument(
        "reads", nargs="?", default="-", help="the filename from which to read the scRNA-seq reads from; defaults to stdin"
    )
    parser.add_argument(
        "-i", "--in-format", choices=['s', 'b'], default=None, help="whether the reads are in the sam or bam file format; defaults to inferring from the file extension or if stdin, defaults to sam"
    )
    parser.add_argument(
        "-f", "--out-format", choices=['s', 'b'], default=None, help="whether the output reads should be written in the sam or bam file format; defaults to inferring from the file extension or if stdout, defaults to sam"
    )
    args = parser.parse_args()

    # parse the optional format args correctly
    args.in_format, args.out_format = parse_format(args.in_format, args.reads), parse_format(args.out_format, args.out)

    new_reads = write_reads(args.out, main(read_barcodes(args.barcodes), args.reads, args.in_format), args.out_format)

