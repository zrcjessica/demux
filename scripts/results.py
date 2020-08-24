#!/usr/bin/env python

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn import metrics
from sklearn.preprocessing import MultiLabelBinarizer


# OUR QUESTIONS
# 1) how many droplets were correctly classified as doublets vs singlets?
# 2) were doublets assigned to the correct samples? (hamming loss)
# 3) were singlets assigned to the correct samples? (also hamming loss)


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
    barcodes['sample'] = barcodes['sample'].apply(lambda samp: tuple(str(s) for s in samp))
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

def hamming_exactmatch(predicts, truth, samples, exact=False):
    """ return the hamming loss (or exact match accuracy) for these droplets """
    # convert to matrix format
    preds_matrix = MultiLabelBinarizer(classes = samples)
    preds_matrix = preds_matrix.fit_transform(predicts['sample'].copy())
    trth_matrix = MultiLabelBinarizer(classes = samples)
    trth_matrix = trth_matrix.fit_transform(truth['sample'].copy())
    if exact:
        return metrics.accuracy_score(trth_matrix, preds_matrix)
    return metrics.hamming_loss(trth_matrix, preds_matrix)

def hammings_accuracy(predicts, truth, exact=False):
    """ return the hamming loss (or exact match accuracy) for doublets, singlets, and both """
    # preprocess the predicts dataframe by converting ambiguous droplets to None
    preds = predicts.copy()
    preds.loc[preds['type'] == 'AMB', 'sample'] = [('NA',)] * sum(preds['type'] == 'AMB')
    # and get the union of the samples from the predicts and the truth
    samples = set(preds['sample'].explode().unique())
    samples |= set(truth['sample'].explode().unique())
    samples = sorted(tuple(samples))
    # get the hamming losses
    labels = {'DBL':0, 'SNG':0, 'BOTH':0}
    for label in labels:
        lab = [label]
        if lab == ['BOTH']:
            lab = ['DBL', 'SNG']
        trth = truth[truth['type'].isin(lab)]
        prds = preds.loc[trth.index].copy()
        labels[label] = hamming_exactmatch(prds, trth, samples, exact)
    return labels

def cohen_kappa(predicts, truth):
    """ calculate cohen's kappa for the singlets """
    preds = predicts.loc[preds['type'] == 'SNG', 'sample']
    trth = truth.loc[preds['type'] == 'SNG', 'sample']
    # agh this won't work because they won't share the same droplets!
    # potential solution: just calculate cohen's kappa among predicted and simulated singlets
    return metrics.cohen_kappa_score(preds, trth)

def main(demux, truth):
    # retrieve the predicted samples from demuxlet
    predicts = get_predicts(demux)
    # retrieve the true samples from the simulation script
    truth = get_truth(truth)
    type_scores = type_metrics(predicts, truth)
    ham_score = hammings_accuracy(predicts, truth)
    accuracy_score = hammings_accuracy(predicts, truth, exact=True)
    return type_scores, ham_score, accuracy_score

def write_out(out, type_scores, ham_score, accuracy_score):
    print("precision/recall:", file=out)
    print(type_scores, file=out)
    print("\nhamming loss:", file=out)
    print(ham_score, file=out)
    print("\nsubset accuracy:", file=out)
    print(accuracy_score, file=out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize results from a demultiplexing simulation.")
    parser.add_argument(
        "demux", type=Path, help="demuxlet's .best file"
    )
    parser.add_argument(
        "truth", type=Path, help="the true labels"
    )
    parser.add_argument(
        "out", type=Path, nargs='?', default=sys.stdout, help="a directory for our output"
    )
    args = parser.parse_args()

    results = main(args.demux, args.truth)
    write_out(sys.stdout, *results)
