#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn import metrics
from sklearn.preprocessing import MultiLabelBinarizer
from IPython import embed


# OUR QUESTIONS
# 1) how many droplets were correctly classified as doublets vs singlets?
# 2) were doublets assigned to the correct samples? (hamming loss)
# 3) were singlets assigned to the correct samples? (also cohen's kappa)


def get_predicts(demux):
    """ parse the .best file into a pd dataframe """
    barcodes = pd.read_csv(demux, sep="\t", usecols=['BARCODE', 'BEST'], index_col='BARCODE')
    # TODO: include probabilities so that you can look at curves?
    barcodes = barcodes['BEST'].str.split('-', n=1, expand=True)
    barcodes.columns = ['type', 'sample']
    barcodes['sample'] = barcodes['sample'].str.split('-')
    # remove the prob suffix from each doublet
    barcodes['sample'][barcodes['type'] == 'DBL'] = barcodes['sample'][barcodes['type'] == 'DBL'].apply(lambda row: row[:2])
    # convert lists to tuples b/c they're immutable
    barcodes['sample'] = barcodes['sample'].apply(tuple)
    return barcodes

def get_truth(truth):
    """ parse the true barcode assignments into an equivalent pd df """
    barcodes = pd.read_csv(truth, sep="\t", header=None, names=['samp', 'old', 'BARCODE'])
    barcodes = barcodes.groupby('BARCODE').apply(lambda grp: list(grp['samp'])).to_frame('sample')
    barcodes.insert(0, 'type', barcodes['sample'].map(lambda x: ['SNG', 'DBL'][len(x)-1]), True)
    # convert lists to tuples b/c they're immutable
    barcodes['sample'] = barcodes['sample'].apply(tuple)
    return barcodes

def type_metrics(predicts, truth):
    """ were droplots correctly classified by their type? """
    labels = ['DBL', 'SNG']
    scores = metrics.precision_recall_fscore_support(
        truth['type'],
        predicts['type'],
        labels=labels
    )
    scores = pd.DataFrame(
        scores, columns=labels,
        index=['precision', 'recall', 'fscore', 'support']
    )
    return scores

def hamming(predicts, truth):
    """ return the hamming loss for these droplets """
    preds = predicts.copy()
    # convert ambiguous droplets to None
    preds['sample'][preds['type'] == 'AMB'] = (['N/A'],)
    samples = preds['sample'].explode().unique()
    # convert to matrix format
    prds = MultiLabelBinarizer(classes = samples)
    prds = prds.fit_transform(preds['sample'])
    trth = MultiLabelBinarizer(classes = samples)
    trth = trth.fit_transform(truth['sample'])
    return metrics.hamming_loss(trth, prds)

def hammings(predicts, truth):
    labels = {'DBL':0, 'SNG':0, 'BOTH':0}
    for label in labels:
        lab = list(label)
        if lab == ['BOTH']:
            lab = ['DBL', 'SNP']
        trth = truth[truth['type'].isin(lab)]
        prds = predicts.loc[trth.index]
        labels[label] = hamming(prds, trth)
    embed()
    return labels

def main(demux, truth, out):
    # retrieve the predicted samples from demuxlet
    predicts = get_predicts(demux)
    # retrieve the true samples from the simulation script
    truth = get_truth(truth)
    type_scores = type_metrics(predicts, truth)
    dbl_score = hammings(predicts, truth)
    return type_scores, dbl_score

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

    results = main(args.demux, args.truth, args.out)
