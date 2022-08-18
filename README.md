# TL-Snakemake
This is a Snakemake implementation of a portion of the [TimeLapse pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by the [Simon lab](https://simonlab.yale.edu/) at Yale. The contributors to the original pipeline are Matthew Simon, Jeremy Schofield, Martin Machyna, Lea Kiefer, and Joshua Zimmer. The original TimeLapse pipeline was developed to process fastq files from TimeLapse-seq (or similar methods, e.g., SLAM-seq, TUC-seq, etc.), map them to a reference genome, and U-to-C or G-to-A mutations in mapped reads. TL-Snakemake currently includes the basic TimeLapse pipeline functionality downstream of alignment. Thus, the input to TL-Snakemake is a set of .bam files and the output is a cB.csv file, which as described on the TimeLapse pipeline bitbucket, contains the following columns:
* XF - Mature feature: ENSEMBL ID if the read mapped solely to exonic parts of a feature
* GF - Gene feature: ENSEMBL ID when read is aligned to any part of feature (intronic or exonic)
* rname - Chromosome name
* nT - Number of uridines in read (additional columns are added if considering other types of mutations, like G-to-A mutations)
* TC - Number of called T-to-C mutations in read (additional columns are added for other types of mutations)
* sj - Logical: TRUE if read contains exon-exon spliced junction
* ai - Logical: TRUE if read aligns to any intronic region
* io - Logical: TRUE if read aligns to only intronic regions
* ei - Logical: TRUE if read aligns to intronic and exonic regions
* sample - Sample name
* n - Number of reads which have the identical set of values described above

## Requirements
TL-Snakemake (as the name implies) uses the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/). The minimal version of Snakemake is techncially compatible with Windows, Mac, and Linux OS, but several of the software dependencies (e.g., HTSeq) are only Mac and Linux compatible. In addition, TL-Snakemake has so far been exclusively tested on a Linux OS, so Mac users be warned. If you are a Windows user like me, don't sweat it, I would suggest looking to the Windows subsystem for linux which can be easily installed (assuming you are running Windows 10 version 2004 or higher).

In addition, you will need Git installed on your system so that you can clone this repository. Head to [this link](https://git-scm.com/downloads) for installation instructions if you don't already have Git.

## Setup
There are 5 steps required to get up and running with TL-Snakemake
1. [Install conda (or mamba) on your system](#conda). This is the package manager that TL-Snakemake uses to make setting up the necessary dependencies a breeze.
1. [Clone this repository](#clone) to wherever you want to run it from
1. [Install dependencies](#depend) with conda/mamba
1. [Edit the config file](#config) (located in config/ directory of cloned repo) to your liking
1. [Run it!](#run)

The remaining documentation will describe each of these steps in greater detail and point you to additional documentation that might be useful.

### Install conda (or mamba)<a name="conda"></a>
[Conda](https://docs.conda.io/projects/conda/en/latest/index.html) is a package/environment management system. [Mamba](https://mamba.readthedocs.io/en/latest/) is a newer, faster, C++ reimplementation of conda. While often associated with Python package management, lots of software, including all of the TimeLapse pipeline dependencies, can be installed with these package managers. They have pretty much the same syntax and can do the same things, so I highly suggest using Mamba in place of conda whenever possible. 

One way to install Mamba is to first install Conda following the instructions at [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Then you can call:
```
$ conda install -n base -c conda-forge mamba
```
to install Mamba.

A second strategy would be to install Mambaforge, which is similar to something called Miniconda but uses Mamba instead of Conda. I will reproduce the instructions to install Mambaforge below, as this is probably the easiest way to get started with the necessary installation of Mamba. These instructions come from the [Snakemake Getting Started tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html), so go to that link if you'd like to see the full original details:

* For Linux users with a 64-bit system, run these two lines of code from the terminal:
``` 
$ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
$ bash Mambaforge-Linux-x86_64.sh
```
* For Mac users with x86_64 architecture: 
``` 
$ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
$ bash Mambaforge-MacOSX-x86_64.sh
```
* And for Mac users with ARM/M1 architecture:
```
$ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh -o Mambaforge-MacOSX-arm64.sh
$ bash Mambaforge-MacOSX-arm64.sh
```

When asked this question:
```
Do you wish the installer to preprend the install location to PATH ...? [yes|no]
```
answer with `yes`. Prepending to PATH means that after closing your current terminal and opening a new one, you can call the `conda` command to install software packages and create isolated environments. We'll be using this in the two steps.

### Clone TL-Snakemake<a name="clone"></a>
Clone the TL-Snakemake repository to wherever you would like on your system. You will eventually be navigating to this repo directory in the terminal and running Snakemake from inside the directory, so make sure your chosen location is conducive to this. Navigate to the directory in the terminal and run:
```
$ git clone https://github.com/isaacvock/TL-Snakemake.git
$ cd TL-Snakemake
```
You should be in the TL-Snakemake repo directory now!

### Install dependencies<a name="depend"></a>
Inside the repo, you will find a file named `current_env.yaml`. This is a YAML file with the list of exact dependencies that need to be installed to run the Snakemake workflow. Luckily, installation of everything can be completed automatically thanks to Mamba/Conda! If you followed the earlier instructions for installing Mamba, just run:
```
$ conda activate base
$ mamba env create --file current_env.yaml
```
The first line activates the so-called `base` conda environment, which is what allows you to call `mamba` in the next line. This will create a new environment called pipeline-env-tmp2. If you don't like that environment name, you can change it by editing the first line of `current_env.yaml`, which looks like:
```
name: pipeline-env-tmp2
```
by replacing `pipeline-env-tmp2` with the name of your choice! 

To make sure we are all on the same page, an environment is a collection of installed software tucked away in its own directory, separate from everything else installed on your computer. When installing software with Mamba/Conda, I highly suggest installing all of the software you need for a particular task into an isolated environment like this. What's great about it is that you can have multiple versions of any given software installed at the same time! 

Now, whenever you want to run the TL-Snakemake pipeline, you can just run:

```
$ conda activate pipeline-env-tmp2
```
replacing `pipeline-env-tmp2` with whatever name you ended up giving the environment, and voila, all of the dependencies are ready to be called upon.

### Edit the config file<a name="config"></a>
In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. The first and arguably most important parameter is at the top of the file:

```
samples:
  WT_1: data/samples/DCP2_subsample1Aligned.out.bam
  WT_2: data/samples/DCP2_subsample2Aligned.out.bam
```
`samples` is the list of sample IDs and paths to .bam files that you want to process. Delete the existing sample names and paths and add yours. The sample names in this example are `WT_1` and `WT_1` and will be what shows up in the `sample` column of the output cB.csv file. The `:` is necessary to distinguish the sample name from what follows, the path to the relevant bam file. Note, the path is NOT an absolute path, it is relative to the top of the TL-Snakemake repo directory. Thus, in this example, the bam files are located in a directory called `samples` that is inside of a directory called `data` located at the top of the TL-Snakemake directory. Your data can be wherever you want it to be, but it might be easiest if you put it in a `data` directory inside the TL-Snakemake directory as in this example. 

As another example, imagine that the `data` directory was in the directory that you cloned the TL-Snakemake repo to, and that there was no `samples` subdirectory inside of `data`. In that case, the paths would look something like this:

```
samples:
  Sample_1: ../data/Sample1Aligned.out.bam
  Sample_2: ../data/Sample2Aligned.out.bam
```
where `../` means navigate up one directory. 

A second almost equally important parmaeter immediately follows:

```
annotation: data/annotation/Homo_sapiens.GRCh38.104.chr_chr.gtf
```
This is the path to the GTF file for the genome that reads were mapped to. The same rules apply when it comes to specifying this path.

The other parameters that can be altered are:
* `cpus`: the number of cores (i.e., cpus) you want to use.
* `FORMAT`: whether the reads are paired-end (PE) or single-end (SE).
* `fragment_size`: For parallel processing, bam files will be split up into temporary files with `fragment_size` reads in each. The default value for this is what was used for testing with a downsampled .bam file, and thus is likely a bit small for typicaly analyses. A value of around (total # of reads)/(cpus) should work more generally. For example, if you have on average 20 million mapped reads in your bam files and will be running this with 20 cores, `fragment_size` should be around 1 million.
* `mut_tracks`: the type of mutation (e.g., T-to-C mutations) that you are most interested in. If T-to-C, then `mut_tracks` should be TC. If G-to-A, then `mut_tracks` should be GA. If both, then `mut_tracks` should be "TC,GA".
* `minqual`: Minimum base quality to call it a mutation. I wouldn't usually worry about editing this.
* `keepcols`: Names of columns to keep in cB.csv output file. I wouldn't usually worry about editing this.

Edit the values in the config file as necessary and move on to the last step.

### Run it!<a name="run"></a>
Technically, all you need to do to run TL-Snakemake is activate the pipeline environment and call snakemake from the top of the TL-Snakemake directory as follows (replacing `pipeline-env-tmp2` with whatever you named the environment and `8` replaced with whatever you entered for `cpus` in the config). Also, TL-Snakemake doesn't currently support use of -s4U control samples to call SNPs, though this will hopefully change soon, so you have to create an empty snp.txt file before running. You should still process these control samples though as they are helpful for downstream analysis and data quality assessment:
```
$ touch snp.txt
$ conda activate pipeline-env-tmp2
$ snakemake --cores 8
```
There are **A LOT** of adjustable parameters that you can play with when running a Snakemake pipeline. I would point you to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) 
for the details on everything you can change when running the pipeline.

## Output
All output files will be placed in a directory named `results` that will be created the first time you run TL-Snakemake. The output of greatest interest, the gzipped cB.csv file, will be in `results/cB/`. 

## Questions?
If you have any questions or run into any problems, feel free to reach out to me (Isaac Vock) at isaac.vock@gmail.com
