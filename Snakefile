import warnings
from pathlib import Path
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.yml"
container: "docker://continuumio/miniconda3:4.8.2"


def check_config(value, place=config, default=False):
    """ return true if config value exists and is true """
    return place[value] if (value in place and place[value]) else default
# parse samples
if check_config('samples'):
    # double check that the user isn't asking for samples they haven't provided
    if not set(config['samples']).issubset(config['data'].keys()):
        warnings.warn("Not all of the samples requested have provided input. Proceeding with as many samples as is possible...")
        config['samples'] = set(config['samples']).intersection(config['data'].keys())
else:
    config['samples'] = list(config['data'].keys())


rule all:
    input:
        expand(
            config['out']+"/new_reads/{samp}.bam", samp=config['samples']
        )

# @zrcjessica: You can add your code here. I hope you don't mind that I've
# already named the rule 'simulate' ;)
# By the way, I've already created the output directive. It's a directory, right?
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#directories-as-outputs
# Each tsv in the directory has to be named sample_name + '.tsv.gz'
rule simulate:
    output: directory(config['out']+"/barcodes")

rule new_bam:
    input:
        barcode_dir = rules.simulate.output[0],
        reads = lambda wildcards: config['data'][wildcards.samp]['reads']
    params:
        barcodes = lambda wildcards, input: str(Path(input.barcode_dir))+"/"+wildcards.samp+".tsv.gz"
    output: config['out']+"/new_reads/{samp}.bam"
    conda: "env.yml"
    shell:
        "scripts/new_bam.py -o {output} {params} {input.reads}"

rule merge:
    input:
        bams = expand(rules.new_bam.output, samp=config['samples'])
    output:
        bam = config['out']+"/new_reads/{samp}.bam"
    conda: "env.yml"
    shell:
        "samtools merge {output} {input}"
