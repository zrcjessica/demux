import warnings
from pathlib import Path
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.yml"
container: "docker://continuumio/miniconda3:4.8.2"


def check_config(value, place=config, default=False):
    """ return true if config value exists and is true """
    return place[value] if (value in place and place[value]) else default

# also parse samples
if check_config('samples'):
    # double check that the user isn't asking for samples they haven't provided
    if not set(config['samples']).issubset(config['data'].keys()):
        warnings.warn("Not all of the samples requested have provided input. Proceeding with as many samples as is possible...")
        config['samples'] = set(config['samples']).intersection(config['data'].keys())
else:
    config['samples'] = list(config['data'].keys())

# and set a default value for the output dir
config['out'] = check_config('out', default='out')


rule all:
    input:
        expand(
            config['out']+"/{rate}/demuxlet/out.best", rate=check_config('rate', default=0.3)
        )

rule unique_barcodes:
    input:
        barcodes = [config['data'][samp]['barcodes'] for samp in config['samples']]
    params:
        samples = config['samples']
    output: directory(config['out'] + "/unique_filtered_barcodes")
    conda: "envs/default.yml"
    shell:
        "mkdir -p {output} && "
        "scripts/get_unique_filtered_barcodes.py -b {input.barcodes} "
        "-s {params.samples} -o {output}"

rule simulate:
    input: 
        barcodes_dir = rules.unique_barcodes.output
    output:
        new_barcodes_dir = directory(config['out'] + "/{rate}/renamed_filtered_barcodes"),
        reference_dir = directory(config['out'] + "/{rate}/reference_tables")
    conda: "envs/default.yml"
    shell:
        "mkdir -p {output} && "
        "scripts/simulate_doublets.py -b {input.barcodes_dir} -d {wildcards.rate} "
        "-o {output.new_barcodes_dir} -r {output.reference_dir}"

rule new_bam:
    input:
        barcode_dir = rules.simulate.output[0],
        reads = lambda wildcards: config['data'][wildcards.samp]['reads']
    params:
        barcodes = lambda wildcards, input: str(Path(input.barcode_dir))+"/"+wildcards.samp+".tsv.gz"
    output: config['out']+"/{rate}/new_reads/{samp}.bam"
    conda: "envs/default.yml"
    shell:
        "scripts/new_bam.py -o {output} {params} {input.reads}"

rule merge:
    input:
        bams = expand(rules.new_bam.output, samp=config['samples'], rate='{rate}')
    output:
        bam = temp(config['out']+"/{rate}/merge.bam")
    conda: "envs/default.yml"
    shell:
        "samtools merge -h {input[0]} -cf {output} {input}"

rule sort:
    input: rules.merge.output.bam
    output:
        bam = config['out']+"/{rate}/merge.sort.bam",
        idx = config['out']+"/{rate}/merge.sort.bam.bai"
    conda: "envs/default.yml"
    shell:
        "samtools sort -o {output.bam} {input} && samtools index {output.bam}"

# TODO: create a rule to handle conflicting UMIs?

rule demux:
    input:
        bam = rules.sort.output.bam,
        vcf = config['vcf']
    params:
        samps = "\\n".join((str(config['data'][samp]['vcf_id']) for samp in config['samples'])),
        out = lambda wildcards, output: output.best[:-len('.best')]
    output:
        best = config['out']+"/{rate}/demuxlet/out.best",
        sing = config['out']+"/{rate}/demuxlet/out.sing",
        sing2 = config['out']+"/{rate}/demuxlet/out.sing2"
    conda: "envs/demuxlet.yml"
    shell:
        "demuxlet --sm-list <(echo -e \"{params.samps}\") --sam {input.bam} --vcf {input.vcf} --out {params.out}"

# TODO: create a rule to summarize the simulation results
