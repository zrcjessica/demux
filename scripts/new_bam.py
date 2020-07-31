#!/usr/bin/env python

import pysam
import argparse
from pathlib import Path


def parse_format(fmt, f):
    """ what input or output format should we use? """
    if fmt is None:
        fmt = 'b' if Path(f).suffix == '.bam' else 's'
    return fmt

def tsv_reader(f):
    """ read lines from a file """
    for line in f:
        yield line.split("\t")

def read_barcodes(tsv):
    """ retrieve the tsv file containing the barcodes """
    with open(tsv) as infile:
        barcodes = {
            row[0]: (row[1] if len(row)-1 and len(row[1]) else row[0])
            for row in tsv_reader(infile)
        }
    return barcodes

def get_reads(reads):
    """ get the reads but first get the tags from the first line"""
    read_iter = iter(reads)
    first = read_iter.next()
    yield first.tags
    yield first
    for read in reads:
        yield read

def main(barcodes, reads, in_format=None, filter_also=False):
    """
        read the reads and alter their barcodes, returning each line as a
        generator function; the first line will be the pysam object
    """
    if in_format != "b":
        in_format = ''
    reads = pysam.AlignmentFile(reads, "r"+in_format)
    yield reads
    reads = get_reads(reads)
    # get the indices of each tag for fast lookup later
    tag_idxs = {}
    tags = next(reads)
    for i in range(len(tags)):
        if tags[i][0] == 'CB':
            tag_idxs['CB'] = i
    # iterate through each read
    for read in reads:
        # initialize a pointer to the tags
        tags = read.tags
        # check to see whether the CB tag needs to be changed
        if tags[tag_idxs['CB']][1] in barcodes:
            # get the new CB tag
            tags[tag_idxs['CB']] = (tags[tag_idxs['CB']][0], barcodes[tags[tag_idxs['CB']][1]])
        elif filter_also:
            continue
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
        "barcodes", type=Path, help="a TSV file containing two columns: 1) barcodes in the input bam file, and 2) the new barcodes if they should be changed, otherwise leave the 2nd column out or leave it empty"
    )
    parser.add_argument(
        "reads", nargs="?", default="-", help="the filename from which to read the scRNA-seq reads from; defaults to stdin"
    )
    parser.add_argument(
        "--in-format", choices=['s', 'b'], default=None, help="whether the reads are in the sam or bam file format; defaults to inferring from the file extension or if stdin, defaults to sam"
    )
    parser.add_argument(
        "--out-format", choices=['s', 'b'], default=None, help="whether the output reads should be written in the sam or bam file format; defaults to inferring from the file extension or if stdout, defaults to sam"
    )
    parser.add_argument(
        "-f", "--filter-also", action='store_true', help="whether to discard reads whose barcodes don't appear in the barcodes file"
    )
    args = parser.parse_args()

    # parse the optional format args correctly
    args.in_format, args.out_format = parse_format(args.in_format, args.reads, args.filter_also), parse_format(args.out_format, args.out)

    new_reads = write_reads(args.out, main(read_barcodes(args.barcodes), args.reads, args.in_format), args.out_format)

