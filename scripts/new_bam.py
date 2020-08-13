#!/usr/bin/env python

import gzip
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
    # first, open a reader for the file
    tsv = gzip.open(tsv) if tsv.suffix == '.gz' else open(tsv)
    with tsv as infile:
        barcodes = {
            row[0]: (row[1] if len(row)-1 and len(row[1]) else row[0])
            for row in tsv_reader(infile)
        }
    return barcodes

def get_reads(reads):
    """ get the reads but first get the tags from the first line"""
    read_iter = iter(reads)
    first = next(read_iter)
    yield first.tags
    yield first
    for read in reads:
        yield read

def main(barcodes, reads, in_format=None, no_filter=False, keep_tags=False):
    """
        read the reads and alter their barcodes, returning each line as a
        generator function; the first line will be the pysam object
    """

    if in_format != "b":
        in_format = ''
    reads = pysam.AlignmentFile(reads, "r"+in_format)

    # parse and output the header
    head = reads.header.to_dict()
    if not keep_tags:
        # delete sample-specific tags
        for tag in ['PG', 'CO']:
            del head[tag]
        # change the RG tag too, so that it is consistent across every sample
        RG_ID = 'Rat:0:1:HFYJTDRXX:1'
        head['RG'] = [{
            'ID': RG_ID,
            'SM': 'Rat',
            'LB': '0.1',
            'PU': 'Rat:0:1:HFYJTDRXX:1',
            'PL': 'ILLUMINA'
        }]
    yield head

    reads = get_reads(reads)
    # get the indices of each tag for fast lookup later
    tag_idxs = {}
    tags = next(reads)
    for i in range(len(tags)):
        if tags[i][0] in ('CB', 'RG', 'PG'):
            tag_idxs[tags[i][0]] = i

    # iterate through each read
    for read in reads:
        tags = read.tags
        # check to see whether the CB tag needs to be changed
        if tags[tag_idxs['CB']][1] in barcodes:
            # get the new CB tag
            tags[tag_idxs['CB']] = ('CB', barcodes[tags[tag_idxs['CB']][1]])
        elif not no_filter:
            continue
        if not keep_tags:
            # also change the RG and PG tags so they are consistent across every sample
            tags[tag_idxs['RG']] = ('RG', RG_ID)
            if 'PG' in tag_idxs:
                tags.pop(tag_idxs['PG'])
        # apply the tags back to the read
        read.tags = tags
        yield read

def write_reads(out, reads, out_format=None):
    """
        read each read from the iterator 'reads' and write it to out
    """
    if out_format != 'b':
        out_format = ''
    out = pysam.AlignmentFile(out, "w"+out_format, header=next(reads))
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
        "-t", "--keep-tags", action='store_true', help="keep PG and CO tags in the header and don't change the RG tags"
    )
    parser.add_argument(
        "-f", "--no-filter", action='store_true', help="keep reads even if their barcodes don't appear in the barcodes file"
    )
    args = parser.parse_args()

    # parse the optional format args correctly
    args.in_format, args.out_format = parse_format(args.in_format, args.reads), parse_format(args.out_format, args.out)

    new_reads = write_reads(args.out, main(read_barcodes(args.barcodes), args.reads, args.in_format, args.no_filter, args.keep_tags), args.out_format)

