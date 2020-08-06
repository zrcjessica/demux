#!/usr/bin/env python3

import pandas as pd
import gzip
import glob
from collections import Counter
import argparse

def main():
	parser = OptionParser()
	parser.add_option("-b","--barcodes",
		action = "store", type = "string", nargs = "+",
		dest = "filtered_barcodes",
		help = "list of filtered barcode files")
	parser.add_option("-s", "--samples",
		action = "store", type = "string", nargs = "+",
		dest = "samples",
		help = "list of sample names corresponding to barcode files")
	parser.add_option("-o","--out",
		action = "store", type = "string",
		dest = "out_dir",
		help = "full path to directory where output files will be written")
	args = parser.parse_args()

	# load barcodes
	all_barcodes = []
	for file in args.filtered_barcodes:
	    with gzip.open(file) as f:
	        sample_barcodes = f.read().splitlines()
	        all_barcodes.append(sample_barcodes)
	all_barcodes_flat = [x for l in all_barcodes for x in l]

	print('%d barcodes across all filtered matrices' % len(all_barcodes_flat))
	print('%d unique barcodes across all filtered matrices' % len(set(all_barcodes_flat)))

	# identify redundant barcodes and which sample they come from
	barcodes_counted = Counter(all_barcodes_flat).items()
	duplicated_barcodes = [bc for bc, count in barcodes_counted if count > 1]
	
	# remove duplicated barcodes from barcodes.tsv.gz files and write to new files

	all_filtered_barcodes = []
	for sample, barcodes in zip(args.samples, all_barcodes):
	    filtered_barcodes = sorted(set(barcodes) - set(duplicated_barcodes))
	    all_filtered_barcodes.append(filtered_barcodes)
	    
	    with gzip.open('%s/%s.tsv.gz' % (args.out_dir, sample), 'wb') as fh:
	        for barcode in filtered_barcodes:
	            fh.write(barcode)
	            fh.write('\n'.encode('utf-8'))
	        fh.close()
main()
