# demux
A pipeline for running single-cell demultiplexing simulations.

# download
Execute the following command.
```
git clone https://github.com/zrcjessica/demux.git
```

# setup
## dependencies
The pipeline is written as a Snakefile which can be executed via [Snakemake](https://snakemake.readthedocs.io). We recommend installing version 5.18.0:
```
conda create -n snakemake -c bioconda -c conda-forge 'snakemake==5.18.0' --no-channel-priority
```
We highly recommend you install [Snakemake via conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda) like this so that you can use the `--use-conda` flag when calling `snakemake` to let it [automatically handle all dependencies](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) of the pipeline. Otherwise, you must manually install the dependencies listed in the [env files](envs).

## input
Symlink your data into the gitignored `data/` folder:
```
ln -s /iblm/netapp/data1/jezhou/Telese_Rat_Amygdala data
```
If you ever need to switch the input to a different dataset, you can just change the symlink path.

# execution
Locally:
```
./run &
```
__or__ on an SGE cluster:
```
qsub run
```

#### Executing the pipeline on your own data
You must modify [the config.yml file](config.yml) to specify paths to your data. The config file is currently configured to run the pipeline on our data (in the git-ignored `data/` folder).

# files and directories
### [Snakefile](Snakefile)
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for running the demultiplexing simulation.

### [config.yml](config.yml)
Config file that defines options and input for the pipeline.

### [scripts/](scripts)
Various scripts used by the pipeline. See the [script README](scripts/README.md) for more information.

### [envs/](envs)
The dependencies of our pipeline, specified as `conda` [environment files](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually). These are used by Snakemake to automatically install our dependencies at runtime.

### [run](run)
An example bash script for executing the pipeline using `snakemake` and `conda`. Any arguments to this script are passed directly to `snakemake`.
