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
# also set default value for output dir
if not check_config('out'):
    config['out'] = 'out'


rule all:
    input:
        expand(
            config['out']+"/merged.bam", samp=config['samples']
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
        bam = config['out']+"/merged.bam"
    conda: "env.yml"
    shell:
        "samtools merge -h {input[0]} -cf {output} {input}"

# TODO: create rules to sort and index the merged.bam file before passing it to demuxlet
# TODO: create a rule to handle conflicting UMIs?

rule demux:
    input:
        bam = rules.merge.output.bam,
        vcf = config['vcf']
    params:
        out = lambda wildcards, output: output.best[:-len('.best')]
    output:
        best = config['out']+"/demuxlet/out.best",
        sing = config['out']+"/demuxlet/out.sing",
        sing2 = config['out']+"/demuxlet/out.sing2",
        pair = config['out']+"/demuxlet/out.pair"
    conda: "env.yml"
    shell:
        "demuxlet --sam {input.bam} --vcf {input.vcf} --out {params.out}"

# TODO: create a rule to summarize the simulation results
