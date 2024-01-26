# Welcome to bam2bakR version 3.0.0!

## Running bam2bakR

There is now a [single website](https://pipelinedocs.readthedocs.io/en/latest/) that hosts documentation for all of the Snakemake pipelines that I have and will develop, including bam2bakR. Information relevant to bam2bakR on this website includes:

1. [Instructions](https://pipelinedocs.readthedocs.io/en/latest/deploy/) on how to deploy and run any of these pipelines. [Instructions](https://pipelinedocs.readthedocs.io/en/latest/simon/) specific to Yale HPC (and somewhat more generally systems that use a SLURM scheduler) are also available.
2. [Details](https://pipelinedocs.readthedocs.io/en/latest/bam2bakR/configuration/) about configuring bam2bakR.
3. [Description](https://pipelinedocs.readthedocs.io/en/latest/bam2bakR/output/) of all output produced by bam2bakR.
4. [FAQs](https://pipelinedocs.readthedocs.io/en/latest/bam2bakR/faqs/).

## Pipeline summary

bam2bakR is a Snakemake implementation of the [TimeLapse pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by the [Simon lab](https://simonlab.yale.edu/) at Yale. The contributors to the original pipeline are Matthew Simon, Jeremy Schofield, Martin Machyna, Lea Kiefer, and Joshua Zimmer. bam2bakR takes either bam files (preferred) or fastq files as input, and produces an easy to work with table that is compatible with an R package developed by Isaac Vock (also the developer of bam2bakR) called [bakR](https://github.com/simonlabcode/bakR).  bakR analyzes and performs comparative analysis of NR-seq data (TimeLapse-seq, SLAM-seq, etc.); see the bakR repository for more details.

Version 3.0.0 of bam2bakR comes with a major change meant to significantly cut down on pipeline runtime. In previous versions, bam2bakR used HTSeq to assign sequencing reads to annotated genomic features (e.g., genes). In version 3.0.0, HTseq was replaced with featureCounts. While a HTSeq might take over an hour to assign all reads in a bam file with 25 million reads to features in a standard human annotation, featureCounts performs the same assignment in a matter of minutes. featureCounts, unlike HTseq, is multi-threaded, and thus can be further sped up by providing it with additional cores. 


## Quick Start

Below is an efficient walk-through of all of the steps necessary to run bam2bakR. A more detailed description of all of these steps can be found [here](https://pipelinedocs.readthedocs.io/en/latest/deploy/). A description of all of the tunable parameters in the bam2bakR config file can be found [here](https://pipelinedocs.readthedocs.io/en/latest/bam2bakR/configuration/).

```
### 
# PREREQUISITES: INSTALL MAMBA AND GIT (only need to do ONCE)
###

# CREATE ENVIRONMENT (only need to do ONCE)
mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy

# CREATE AND NAVIGATE TO WORKING DIRECTORY (only need to do ONCE)
mkdir path/to/working/directory
cd path/to/working/directory

# DEPLOY PIPELINE TO YOUR WORKING DIRECTORY (only need to do ONCE)
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/simonlabcode/bam2bakR.git . --branch main

###
# EDIT CONFIG FILE (need to do ONCE PER NEW DATASET)
###

# RUN PIPELINE

# See [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for details on all of the configurable parameters
snakemake --cores all --use-conda --rerun-triggers mtime

```

## Citations

Please cite the following if you end up using bam2bakR in published work:

[TimeLapse-seq paper](https://www.nature.com/articles/nmeth.4582), where initial pipeline was introduced:

- Schofield JA, Duffy EE, Kiefer L, Sullivan MC, and Simon MD. 2018. TimeLapse-seq: adding a temporal dimension to RNA sequencing through nucleoside recoding. *Nature Methods*. **15**:221-225. doi:10.1038/nmeth.4582.

[bakR paper](https://rnajournal.cshlp.org/content/29/7/958.abstract), where Snakemake implementation was introduced:

- Vock IW and Simon MD. 2023. bakR: uncovering differential RNA synthesis and degradation kinetics transcriptome-wide with Bayesian hierarchical modeling. *RNA*:*rna.079451.122*. doi:10.1261/rna.079451.122.


## Questions?
If you have any questions or run into any problems, feel free to post them to Issues.
