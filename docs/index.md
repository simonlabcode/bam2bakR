# Welcome to bam2bakR

bam2bakR is a Snakemake implementation of a portion of the [TimeLapse pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by the [Simon lab](https://simonlab.yale.edu/) at Yale. The contributors to the original pipeline are Matthew Simon, Jeremy Schofield, Martin Machyna, Lea Kiefer, and Joshua Zimmer.

## What bam2bakR does

bam2bakR currently includes the basic TimeLapse pipeline functionality downstream of alignment. Thus, the input to bam2bakR is a set of .bam files and the output is a cB.csv file, which as described on the TimeLapse pipeline bitbucket, contains by default the following columns:

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

Additional columns can be included and any of these columns can be omitted per your choice.

## Requirements for bam2bakR

bam2bakR uses the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/). The minimal version of Snakemake is techncially compatible with Windows, Mac, and Linux OS, but several of the software dependencies are only Mac and Linux compatible. If you are a Windows user like me, don't sweat it, I would suggest looking to the Windows subsystem for linux which can be easily installed (assuming you are running Windows 10 version 2004 or higher).

In addition, you will need Git installed on your system so that you can clone this repository. Head to [this link](https://git-scm.com/downloads) for installation instructions if you don't already have Git.

Finally, bam2bakR requires bam files from stranded library preps (i.e., information about whether a given read represents the original sequence of the RNA or its reverse complement must be retained). In addition, make sure that the aligner you used was configured to output bam files with an MD tag, which will keep track of the reference nucleotide at any mismatches between the aligned read and the reference genome. While this tag is not always included by default, most aligners have an option to add this tag to the list of default tags.


## Getting Started

There are a number of ways to use bam2bakR:

1. [Deploying with Snakedeploy](deploy.md) (recommended)
1. [Cloning the repo locally](alt.md); there are several possible variants of this option:
    * Running with --use-conda option
    * Running within full environment
    * Running inside a Docker container