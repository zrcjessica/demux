#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path


# OUR QUESTIONS
# 1) how many droplets were correctly classified as doublets vs singlets?
# 2) were doublets assigned to the correct samples? (cohen's kappa)
# 3) were singlets assigned to the correct samples? (also cohen's kappa)


def get_predicts(demux):
    """ parse the .best file into a pd dataframe """
    barcodes = pd.read_csv(demux, sep="\t", usecols=['BARCODE', 'BEST'], index_col='BARCODE')
    barcodes = barcodes['BEST'].str.split('-', n=1, expand=True)
    barcodes.columns = ['type', 'sample']
    barcodes['sample'] = barcodes['sample'].str.split('-')
    # remove the prob suffix from each doublet
    barcodes['sample'][barcodes['type'] == 'DBL'] = barcodes['sample'][barcodes['type'] == 'DBL'].apply(lambda row: row[:2])
    barcodes
    return barcodes

def get_truth(truth):
    """ parse the true barcode assignments into an equivalent pd df """
    truth = pd.read_csv(truth, sep="\t", header=None, names=['samp', 'old', 'BARCODE'])
    truth = truth.groupby('BARCODE').apply(lambda grp: list(grp['samp'])).to_frame('sample')
    truth.insert(0, 'type', truth['sample'].map(lambda x: ['SNG', 'DBL'][len(x)-1]), True)
    return truth

def main(demux, truth, out):
    # retrieve the predicted samples from demuxlet
    predicts = get_predicts(demux)
    # retrieve the true samples from the simulation script
    truth = get_truth(truth)
    return truth, predicts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize results from a demultiplexing simulation.")
    parser.add_argument(
        "demux", type=Path, help="demuxlet's .best file"
    )
    parser.add_argument(
        "truth", type=Path, help="the true labels"
    )
    parser.add_argument(
        "out", type=Path, help="a directory for our output"
    )
    args = parser.parse_args()

    barcodes = main(args.demux, args.truth, args.out)
