## Setup

There are 4 steps required to get up and running with bam2bakR

1. [Install conda (or mamba) on your system](#conda). This is the package manager that bam2bakR uses to make setting up the necessary dependencies a breeze.
1. [Deploy workflow](#deploy) with [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html)
1. [Edit the config file](#config) (located in config/ directory of deployed/cloned repo) to your liking
1. [Run it!](#run)

The remaining documentation on this page will describe each of these steps in greater detail and point you to additional documentation that might be useful.

### Install conda (or mamba)<a name="conda"></a>
[Conda](https://docs.conda.io/projects/conda/en/latest/index.html) is a package/environment management system. [Mamba](https://mamba.readthedocs.io/en/latest/) is a newer, faster, C++ reimplementation of conda. While often associated with Python package management, lots of software, including all of the TimeLapse pipeline dependencies, can be installed with these package managers. They have pretty much the same syntax and can do the same things, so I highly suggest using Mamba in place of Conda whenever possible. 

One way to install Mamba is to first install Conda following the instructions at [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Then you can call:

``` bash
conda install -n base -c conda-forge mamba
```
to install Mamba.

A second strategy would be to install Mambaforge, which is similar to something called Miniconda but uses Mamba instead of Conda. I will reproduce the instructions to install Mambaforge below, as this is probably the easiest way to get started with the necessary installation of Mamba. These instructions come from the [Snakemake Getting Started tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html), so go to that link if you'd like to see the full original details:

* For Linux users with a 64-bit system, run these two lines of code from the terminal:

``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
* For Mac users with x86_64 architecture: 
``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
bash Mambaforge-MacOSX-x86_64.sh
```
* And for Mac users with ARM/M1 architecture:
``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh -o Mambaforge-MacOSX-arm64.sh
bash Mambaforge-MacOSX-arm64.sh
```

When asked this question:
``` bash
Do you wish the installer to preprend the install location to PATH ...? [yes|no]
```
answer with `yes`. Prepending to PATH means that after closing your current terminal and opening a new one, you can call the `mamba` (or `conda`) command to install software packages and create isolated environments. We'll be using this in the next step.

### Deploy workflow<a name="deploy"></a>

Version 1.0.1 of bam2bakR is now compatible with deployment using the tool [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html). To get started with Snakedeploy, you first need to create a simple conda environment with Snakemake and Snakedeploy:


``` bash
mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy
```

Next, create a directory that you want to run bam2bakR in (I'll refer to it as `workdir`) and move into it:
``` bash
mkdir workdir
cd workdir
```

Now, activate the `deploy_snakemake` environment and deploy the workflow as follows:

``` bash
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/simonlabcode/bam2bakR.git . --branch main
```

`snakedeploy deploy-workflow https://github.com/simonlabcode/bam2bakR.git` copies the content of the `config` directory in the bam2bakR Github repo into the directoy specified (`.`, which means current directory, i.e., `workdir` in this example). It also creates a directory called `workflow` that contains a singular Snakefile that instructs Snakemake to use the workflow hosted on the main branch (that is what `--branch main` determines) of the bam2bakR Github repo. `--branch main` can also be replaced with `--tag 1.0.1` to ensure that you are consistently using the same version of bam2bakR (version 1.0.1 release).

### Edit the config file<a name="config"></a>
In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. The first parameter that you have to set is at the top of the file:

``` yaml
samples:
  WT_1: data/samples/WT_replicate_1.bam
  WT_2: data/samples/WT_replicate_2.bam
  WT_ctl: data/samples/WT_nos4U.bam
  KO_1: data/samples/KO_replicate_1.bam
  KO_2: data/samples/KO_replicate_2.bam
  KO_ctl: data/samples/KO_nos4U.bam
```
`samples` is the list of sample IDs and paths to .bam files that you want to process. Delete the existing sample names and paths and add yours. The sample names in this example are `WT_1`, `WT_2`, `WT_ctl`, `KO_1`, `KO_2`, and `KO_ctl`. These are the sample names that will show up in the `sample` column of the output cB.csv file. The `:` is necessary to distinguish the sample name from what follows, the path to the relevant bam file. Note, the path is NOT an absolute path, it is relative to the directory that you deployed to (i.e., `workdir` in this example). Thus, in this example, the bam files are located in a directory called `samples` that is inside of a directory called `data` located in `workdir`. Your data can be wherever you want it to be, but it might be easiest if you put it in a `data` directory inside the bam2bakR directory as in this example. 

As another example, imagine that the `data` directory was in the directory that contains `workdir`, and that there was no `samples` subdirectory inside of `data`. In that case, the paths would look something like this:

``` yaml
samples:
  WT_1: ../data/WT_replicate_1.bam
  WT_2: ../data/WT_replicate_2.bam
  WT_ctl: ../data/WT_nos4U.bam
  KO_1: ../data/KO_replicate_1.bam
  KO_2: ../data/KO_replicate_2.bam
  KO_ctl: ../data/KO_nos4U.bam
```
where `../` means navigate up one directory. 

The next parameter you have to set denotes the sample names of any -s4U control samples (i.e., samples that were not fed s4U or a similar metabolic label):

``` yaml
control_samples: ["WT_ctl", "KO_ctl"]
```

In this case, the samples named WT_ctl and KO_ctl are the -s4U control samples. -s4U controls will be used to call any single nucleotide polymorphisms (SNPs) in your cell line so as to avoid conflating them with T-to-C mutations induced by the nucleotide recoding chemistry. 

The third crucial parmaeter immediately follows:

``` yaml
annotation: data/annotation/GRCh38.gtf
```
This is the path to the GTF annotation file for the genome that reads were mapped to. The same rules apply when it comes to specifying this path.

Finally, the path to the genome fasta file that you used for alignment must also be specified:

``` yaml
genome_fasta: data/genome/GRCh38.fasta
```

The other parameters that can be altered are:

* `strandedness`: whether the first read in a pair (or the only read if single-end) represents the original sequence of the RNA (F), or its reverse complement (R). For example, set this parameter to "F" if your library is an FR paired-end library, and "R" if it is an RF paired-end library.
* `FORMAT`: whether the reads are paired-end (PE) or single-end (SE).
* `mut_tracks`: the type of mutation (e.g., T-to-C mutations) that sequencing tracks will be colored by. If you are most interested in the T-to-C mutational content, then `mut_tracks` should be TC. If G-to-A, then `mut_tracks` should be GA. If both, then `mut_tracks` should be "TC,GA".
* `minqual`: Minimum base quality to call it a mutation.
* `keepcols`: Names of columns to keep in cB.csv output file. See [Output](output.md) for details of columns you can keep.
* `spikename`: If spike-ins are present, this should be a string that is common to all gene_ids for spike-in transcripts in annotation gtf. For example, in Ensembl annotations for Drosophila melanogaster, all gene_ids start with "FBgn". Therefore, if you have Drosophila spike-ins, `spikename` should be "FBgn".
* `normalize`: If True, then scale factor calculated with edgeR is used to normalize sequencing tracks.
* `WSL`: whether you are running this on the Windows subsystem for linux (0 = yes; 1= no)

 
Edit the values in the config file as necessary and move on to the last step.

### Run it!<a name="run"></a>

Once steps 1-3 are complete, bam2bakR can be run from the directory you deployed the workflow to as follows:

``` bash
snakemake --cores all --use-conda
```
There are **A LOT** of adjustable parameters that you can play with when running a Snakemake pipeline. I would point you to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) 
for the details on everything you can change when running the pipeline.

