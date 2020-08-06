# scripts
This directory contains various scripts used by the pipeline.
However, you can use most of these scripts on their own, too. Some may even be helpful in day-to-day use.

All python scripts implement the `--help` argument.

### [get_unique_filtered_barcodes.py](get_unique_filtered_barcodes.py)
A python script that filters out shared cellular barcodes among the samples.

### [new_bam.py](new_bam.py)
A python script that changes cellular barcodes in a BAM file.

### [simulate_doublets.py](simulate_doublets.py)
A python script that simulates doublets by deciding how to rename barcodes.
