#!/usr/bin/env python3

import numpy as np
import gzip
import glob
import random
from scipy import stats
import argparse

def main():
	parser = argparse.ArgumentParser(description = 'process inputs')
	parser.add_argument("-b","--barcodes",
		action = "store", 
		dest = "barcodes_dir",
		help = "full path to directory with all unique filtered feature barcodes for all samples to multiplex")
	parser.add_argument("-d", "--doublets",
		action = "store", type = "float",
		dest = "d",
		help = "doublet rate")
	parser.add_argument("-o", "--out",
		action = "store",
		dest = "out_dir",
		help = "full path to directory to save new barcode files in")
	parser.add_argument("-r", "--ref",
		action = "store",
		dest = "ref_dir",
		help = "full path to directory to save reference tables in")

	args = parser.parse_args()
	print('args parsed')
	print(args)
	
	# load barcodes and sample names
	barcode_files = glob.glob('%s/*' % args.barcodes_dir)
	sample_names = []
	barcodes = []

	for file in barcode_files:
		sample = file.split('.')[0].split('/')[-1]
		sample_names.append(sample)
		with gzip.open(file) as fh:
			sample_barcodes = fh.read().splitlines()
	#         sample_barcodes = [bc.decode('utf-8') for bc in sample_barcodes]
			barcodes.append(sample_barcodes)
	
	print(barcode_files)
	print(sample_names)

	# compute metrics
	x = len([x for l in barcodes for x in l])
	N = len(sample_names)

	print('total cells = %d\ntotal samples = %d' % (x, N))

	# compute expected number of doublets for simulation
	d = opt.d

	totDoublets = int(round(d*x/2))
	a = int(round(totDoublets*(1-1/N)))
	b = totDoublets-a

	print('expected number of doublets = %d\n' % totDoublets)
	print('expected number of doublets containing cells from different individuals = %d\n' % a)
	print('expeced nmber of doublets containing cels from the same individual = %d\n' % b)

	# shuffle list of barcodes for each sample
	shuffled_barcodes = [random.sample(bc, len(bc)) for bc in barcodes]

	# simulate doublets containing cells from the same sample 
	weights = [len(bc)/x for bc in barcodes]
	doublets_same = np.unique(np.random.choice(len(sample_names), b, p = weights), return_counts = True)
	simulated_doublets_same = []
	for sample, counts in zip(doublets_same[0], doublets_same[1]):
		for i in range(counts):
			bcA = shuffled_barcodes[sample].pop()
			bcB = shuffled_barcodes[sample].pop()
			doublet = [bcA, bcB]
			simulated_doublets_same.append(doublet)

	# simulate doublets containing cells from different samples
	doublets_diff = []
	for i in range(a):
		dbl = np.random.choice(len(sample_names),2,p = weights,replace = False)
		doublets_diff.append(dbl)

	doublets_diff = np.array(doublets_diff)

	simulated_doublets_diff = []
	for pair in doublets_diff:
		bcA = shuffled_barcodes[pair[0]].pop()
		bcB = shuffled_barcodes[pair[1]].pop()
		doublet = [bcA, bcB]
		simulated_doublets_diff.append(doublet)

	# change barcodes for doublets
	# make dict mapping second barcode in each doublet to first barcode
	new_barcodes_map = {}
	for doublet in simulated_doublets_same + simulated_doublets_diff:
		new_barcodes_map[doublet[1]] = doublet[0]

	# generate new filtered.tsv.gz files
	for sample,bcs in zip(sample_names, barcodes):
		print('generating new barcodes for %s\n' % sample)
		new_bcs = []
		mismatches = 0
		for i in range(len(bcs)):
			if bcs[i] in new_barcodes_map.keys():
				new_bcs.append(new_barcodes_map[bcs[i]])
				mismatches += 1
			else:
				new_bcs.append(bcs[i])
		print('%d mismatches observed\n' % mismatches)
		new_bcs = sorted(new_bcs)
		with gzip.open('%s/%s.tsv.gz' % (opt.out_dir, sample), 'wb') as fh:
			for bc in new_bcs:
				fh.write(bc)
				fh.write('\n'.encode('utf-8'))
			fh.close()

	# check that number of mismatches is as expected
	expected_mismatches = doublets_same[1] + np.unique(doublets_diff[:,1], return_counts = True)[1]
	for sample, expected in zip(sample_names, expected_mismatches):
		print('%d mismatches expected for %s\n' % (expected, sample))

	# write reference files for simulated doublets
	# for doublets with cells from same samples
	doublets_same_reference = np.concatenate((np.repeat(sample_names, doublets_same[1]).reshape(-1,1),
				   np.array(simulated_doublets_same)), axis = 1)

	with open(opt.ref_dir + '/same_sample_doublets_reference.tsv', 'w') as fh:
		fh.write('sample\tbarcodeA\tbarcodeB\n')
		for l in doublets_same_reference:
			fh.write('\t'.join(l) + '\n')

	# for doublets with cells from different samples
	map_samples = dict(zip(range(len(sample_names)), sample_names))
	doublets_diff_with_sample_names = []
	for i in range(doublets_diff.shape[0]):
		sampleA = map_samples[doublets_diff[i,0]]
		sampleB = map_samples[doublets_diff[i,1]]
		doublets_diff_with_sample_names.append([sampleA, sampleB])

	doublets_diff_reference = np.concatenate((doublets_diff_with_sample_names, simulated_doublets_diff), axis = 1)

	with open(opt.ref_dir + '/diff_sample_doublets_reference.tsv', 'w') as fh:
		fh.write('sampleA\tsampleB\tbarcodeA\tbarcodeB\n')
		for l in doublets_diff_reference:
			fh.write('\t'.join(l) + '\n')

main()
