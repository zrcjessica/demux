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
            config['out']+"/{rate}/results", rate=check_config('rate', default=0.3)
        )

rule unique_barcodes:
    input:
        barcodes = [config['data'][samp]['barcodes'] for samp in config['samples']]
    params:
        samples = config['samples'],
        output_dir = lambda wildcards, output: Path(output[0]).parent
    output:
        expand(config['out'] + "/unique_filtered_barcodes/{samp}.tsv.gz", samp=config['samples'])
    conda: "envs/default.yml"
    shell:
        "mkdir -p {output} && "
        "scripts/get_unique_filtered_barcodes.py -b {input.barcodes} "
        "-s {params.samples} -o {params.output_dir}"

rule simulate:
    input: rules.unique_barcodes.output
    params:
        barcodes_dir = lambda wildcards, input: Path(input[0]).parent,
        reference_dir = lambda wildcards, output: Path(output.diff_reference).parent
    output:
        new_barcodes_dir = directory(config['out'] + "/{rate}/renamed_filtered_barcodes"),
        diff_reference = config['out'] + "/{rate}/reference_tables/diff_sample_doublets_reference.tsv",
        same_reference = config['out'] + "/{rate}/reference_tables/same_sample_doublets_reference.tsv"
    conda: "envs/default.yml"
    shell:
        "mkdir -p {output.new_barcodes_dir} {params.reference_dir} && "
        "scripts/simulate_doublets.py -b {params.barcodes_dir} -d {wildcards.rate} "
        "-o {output.new_barcodes_dir} -r {params.reference_dir}"

rule table:
    """ create a table containing all of the barcodes before and after """
    input:
        old = rules.unique_barcodes.output,
        new = rules.simulate.output.new_barcodes_dir,
    params:
        new = lambda wildcards, input: expand(input.new+"/{samp}.tsv.gz", samp=config['samples'])
    output:
        config['out'] + "/{rate}/barcodes_table.tsv"
    shell:
        "cat "+' '.join([
            "<(paste <(zcat {input.old["+i+"]:q}) <(zcat {params.new["+i+"]:q})"+\
            " | sed 's/^/{config[data]["+config['samples'][int(i)]+"][vcf_id]}\\t/g')"
            for i in map(lambda x: str(x), range(len(config['samples'])))
        ])+" >{output}"

rule new_bam:
    input:
        barcodes = rules.table.output,
        reads = lambda wildcards: config['data'][wildcards.samp]['reads']
    params:
        vcf_id = lambda wildcards: config['data'][wildcards.samp]['vcf_id']
    output: config['out']+"/{rate}/new_reads/{samp}.bam"
    conda: "envs/default.yml"
    shell:
        "scripts/new_bam.py -o {output} "
        "<(grep -P '^{params.vcf_id}\\t' {input.barcodes} | cut -f 2-) {input.reads}"

rule merge:
    input:
        bams = expand(rules.new_bam.output, samp=config['samples'], rate='{rate}')
    output:
        bam = temp(config['out']+"/{rate}/merge.bam")
    conda: "envs/default.yml"
    threads: 12
    shell:
        "samtools merge -@ {threads} -h {input[0]} -cf {output} {input}"

rule sort:
    input: rules.merge.output.bam
    output:
        bam = config['out']+"/{rate}/merge.sort.bam",
        idx = config['out']+"/{rate}/merge.sort.bam.bai"
    conda: "envs/default.yml"
    threads: 12
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input} && samtools index {output.bam}"

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
        sing = config['out']+"/{rate}/demuxlet/out.single",
        sing2 = config['out']+"/{rate}/demuxlet/out.sing2"
    conda: "envs/demuxlet.yml"
    shell:
        "demuxlet --sm-list <(echo -e \"{params.samps}\") --sam {input.bam} --vcf {input.vcf} --out {params.out}"

rule results:
    input:
        rules.demux.output.best,
        rules.table.output
    output: directory(config['out']+"/{rate}/results")
    conda: "envs/default.yml"
    shell:
        "scripts/results.py {input} {output}"
