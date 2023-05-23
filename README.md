# bam2bakR version 2.0.0

## Check out [bam2bakR's website](https://tl-snakemake.readthedocs.io/en/latest/) for revamped and extended documentation!

## bam2bakR now includes optional fastq2bakR functionality. That is, you can now provide fastq files as input and adapter trimming + alignment will be performed

This is a Snakemake implementation of the [TimeLapse pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by the [Simon lab](https://simonlab.yale.edu/) at Yale. The contributors to the original pipeline are Matthew Simon, Jeremy Schofield, Martin Machyna, Lea Kiefer, and Joshua Zimmer. The original TimeLapse pipeline was developed to process fastq files from TimeLapse-seq (or similar methods, e.g., SLAM-seq, TUC-seq, etc.), map them to a reference genome, and U-to-C or G-to-A mutations in mapped reads. bam2bakR currently includes the basic TimeLapse pipeline functionality downstream of alignment. Thus, the input to bam2bakR is a set of .bam or .fastq files and the output is a cB.csv file, which as described on the TimeLapse pipeline bitbucket, contains the following columns by default:
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
bam2bakR uses the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/). The minimal version of Snakemake is techncially compatible with Windows, Mac, and Linux OS, but several of the software dependencies (e.g., HTSeq) are only Mac and Linux compatible. If you are a Windows user like me, don't sweat it, I would suggest looking to the Windows subsystem for linux which can be easily installed (assuming you are running Windows 10 version 2004 or higher).

In addition, you will need Git installed on your system so that you can clone this repository. Head to [this link](https://git-scm.com/downloads) for installation instructions if you don't already have Git.

Finally, bam2bakR requires bam or fastq files from stranded library preps (i.e., information about whether a given read represents the original sequence of the RNA or its reverse complement must be retained). In addition, make sure that the aligner you used was configured to output bam files with an MD tag, which will keep track of the reference nucleotide at any mismatches between the aligned read and the reference genome. While this tag is not always included by default, most aligners have an option to add this tag to the list of default tags.

## Setup for bam2bakR
There are 4 steps required to get up and running with bam2bakR

1. [Install conda (or mamba) on your system](#conda). This is the package manager that bam2bakR uses to make setting up the necessary dependencies a breeze.
1. [Deploy workflow](#deploy) with [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html)
1. [Edit the config file](#config) (located in config/ directory of deployed/cloned repo) to your liking
1. [Run it!](#run)

Each of these steps will be described in greater detail below.

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

`snakedeploy deploy-workflow https://github.com/simonlabcode/bam2bakR.git` copies the content of the `config` directory in the bam2bakR Github repo into the directoy specified (`.`, which means current directory, i.e., `workdir` in this example). It also creates a directory called `workflow` that contains a singular Snakefile that instructs Snakemake to use the workflow hosted on the main branch (that is what `--branch main` determines) of the bam2bakR Github repo. `--branch main` can also be replaced with `--tag 1.0.2` to ensure that you are consistently using the same version of bam2bakR (version 1.0.2 release).

### Edit the config file<a name="config"></a>
In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. The first parameter that you have to set is at the top of the file:

``` yaml
samples:
  WT_1: data/bam/WT_replicate_1.bam
  WT_2: data/bam/WT_replicate_2.bam
  WT_ctl: data/bam/WT_nos4U.bam
  KO_1: data/bam/KO_replicate_1.bam
  KO_2: data/bam/KO_replicate_2.bam
  KO_ctl: data/bam/KO_nos4U.bam
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
* `keepcols`: Names of columns to keep in cB.csv output file. 
* `spikename`: If spike-ins are present, this should be a string that is common to all gene_ids for spike-in transcripts in annotation gtf. For example, in Ensembl annotations for Drosophila melanogaster, all gene_ids start with "FBgn". Therefore, if you have Drosophila spike-ins, `spikename` should be "FBgn".
* `normalize`: If True, then scale factor calculated with edgeR is used to normalize sequencing tracks.
* `WSL`: whether you are running this on the Windows subsystem for linux (0 = yes; 1= no)

 
Edit the values in the config file as necessary and move on to the last step.

### Run it!<a name="run"></a>

Once steps 1-3 are complete, bam2bakR can be run from the directory you deployed the workflow to as follows:

``` bash
snakemake --cores all --use-conda
```
There are **A LOT** of adjustable parameters that you can play with when running a Snakemake pipeline. I will point you to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) 
for the details on everything you can change when running the pipeline.

## Setup for fastq2bakR

fastq2bakR is an optional add-on to the bam2bakR workflow that includes adapter trimming and alignment of fastq files. Adapter trimming is performed with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and alignment can be performed with one of three aligners:

1. [HISAT2](http://daehwankimlab.github.io/hisat2/) (default)
2. [STAR](https://github.com/alexdobin/STAR)
3. [HISAT-3N](http://daehwankimlab.github.io/hisat2/hisat-3n/) (requires manual installation of HISAT-3N prior to running)

Running fastq2bakR requires all of the same setup as running bam2bakR. The difference is that there are additional config file parameters that you may need to change to run fastq2bakR.

Most importantly, you will need to change the value of `bam2bakr` from `True` (it's default value) to `False`.

Secondly, the input (`samples`) is a little different. There is an example of what the input will look like for fastq2bakR in a commented out code block just below the bam2bakR example input. The major difference in that the instead of providing paths to individual files, you now have to provide paths to directories. Each fastq file or pair of fastq files must be in its own directory. 

Finally, there are a number of parameters below the line that reads `##### Parameters that are only relevant if bam2bakr is False #####`. These parameters are:

1. `HISAT2`: If using hisat2 for alignment, this should represent the path to the directory containing pre-built hisat2 alignment indices. As with other directory paths specified in the config file, this path is relative to the directory containing the config and workflow directories (the directory you deployed the workflow inside).
2. `HISAT_3N`: If using HISAT-3N for alignment, this should represent the path to HISAT-3N alignment indices. NOTE: unlike for `HISAT2`, this should be the path to the directory containing the indices, **plus** (as described in the [HISAT-3N documentation](http://daehwankimlab.github.io/hisat2/hisat-3n/)) the name of the index files up to but not including the suffix `.3n.*.*.ht2`. 
3. `STAR_index`: If using STAR for alignment, this should represent the path to STAR alignment indices. fastq2bakR also allows you to build STAR indices as part of the workflow, in which case this will be the directory at which the indices are created.
4. `use_hisat3n`: If True, HISAT-3N will be used for alignment. **NOTE:** Unlike all other dependencies, HISAT-3N cannot be installed automatically via conda. Therefore, you will have to install HISAT-3N yourself to use it as an aligner.
5. `use_star`: If True, STAR will be used for alignment.
6. `build_star`: If True, STAR alignment indices will be created as part of the workflow
7. `hisat3n_path`: As described in the description of the `use_hisat3n` parameter, you need to install HISAT-3N yourself if you want to use it as a part of fastq2bakR. `hisat3n_path` is thus the path to the executable to run hisat-3n (or just `hisat-3n` if HISAT-3N has been added to the PATH environment variable).
8. `chr_tag`: If True, "chr" will be appended to chromosome names in bam files output by HISAT-2 or HISAT-3N. If your annotation/genome files have "chr" in the chromosome names (e.g., chr1, chr2, etc.), then you will need to set this to TRUE.
9. `Yale`: If True, then HISAT-3N will be loaded as a module via the Lmod module system. The parameter is called `Yale` because Yale HPC uses an Lmod system, and the name of the module is hard coded in the workflow as `HISAT-3N`. In other words, HISAT-3N is loaded with `module load HISAT-3N` if this parameter is True.
10. `flattened`: If TRUE, fastq2bakR will expect that the provided annotation file (`annotation` is flattened by DEXSeq).
11. `adapter`: Code passed to cutadapt specifying adapters to be trimmed
12. `cutadapt_extra`: Additional parameters to pass to cutadapt
13. `star_extra`: Additional parameters to pass to STAR
14. `hisat2_extra`: Additional parameters to pass to HISAT2

## Common problems

If there are very few T-to-C mutations in the final cB.csv file (e.g., if sample-wide mutation rates in +s4U samples are < 0.003), then you may have used the incorrect value for the `strandedness` parameter in the config. One way to tell if this is the case is by looking at one of the +s4U sample counts.csv files in `results/counts/` and checking for an abundance of A-to-G mutations. If this is the case, flip the value of `strandedness` to the opposite of whatever you used.

Related to the first point, a good sanity check after running the pipeline is going into R and checking the raw mutation rates as such:

```
library(data.table)

# To unzip and read cB, also need to have R.utils package installed
cB <- fread("path/to/cB.csv.gz")

# Assess sample-wide T-to-C mutation rate in each sample
cB[,.(mutrate = sum(TC*n)/sum(nT*n), by = sample]
  # Want to see that +s4U samples has higher mutation rate than -s4U samples
```

Similarly, checking a counts.csv file for an abundance of A-to-G mutations can be done as follows:

```
library(data.table)

counts <- fread("path/to/+s4U/counts.csv.gz")

## Check if A-to-G mutation rate is higher than T-to-C mutation rate:

# A-to-G mutation rate
sum(counts$AG)/sum(counts$nA)

# T-to-C mutation rate
sum(counts$TC)/sum(counts$nT)
```

## Questions?
If you have any questions or run into any problems, feel free to post them to Issues.
