import glob

### Sample details
SAMP_NAMES = list(config["samples"].keys())

CTL_NAMES = list(config["control_samples"])

s4U_SAMPS = list(set(SAMP_NAMES) - set(CTL_NAMES))

nctl = len(CTL_NAMES)
nsamps = len(SAMP_NAMES)


### Anticipate output of genome indexing
def get_index_name():
    genome = config["genome_fasta"]
    index = str(genome) + ".fai"
    return index


### Config/other parameters used in multiple places

FORMAT = config["FORMAT"]

ALIGNER = config["use_hisat3n"]

STAR = config["use_star"]

NORMALIZE = config["normalize"]

PAIRS = [1, 2]

# Get mutation types to track
MutTypes = config["mut_tracks"]
Mutation_Types = MutTypes.split(",")


### Functions for returning input to first rules


# Relevant if providing bam files as input
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


# Relevant if providing fastq files as input
def get_input_fastqs(wildcards):
    if config["samples"][wildcards.sample].endswith("/"):
        fastq_path = str(config["samples"][wildcards.sample])
        fastq_path = fastq_path[:-1]
    else:
        fastq_path = str(config["samples"][wildcards.sample])

    fastq_path = config["samples"][wildcards.sample]
    fastq_files = sorted(glob.glob("{}/*.fastq*".format(fastq_path)))
    return fastq_files


### FeatureCounts parameters


## All of the files to merge
def get_merge_input(wildcards):
    MERGE_INPUT = []

    MERGE_INPUT.extend(
        expand("results/counts/{SID}_counts.csv.gz", SID=wildcards.sample)
    )

    if config["features"]["genes"]:
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_genes/{SID}.featureCounts", SID=wildcards.sample
            )
        )

    if config["features"]["exons"]:
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_exons/{SID}.featureCounts", SID=wildcards.sample
            )
        )

    return MERGE_INPUT


# Columns in final cB
keepcols = ["sample", "sj", "rname"]


if config["features"]["genes"]:
    keepcols.append("GF")

if config["features"]["exons"]:
    keepcols.append("XF")


keepcols = ",".join(keepcols)


# Get strandedness parameter
if config["strandedness"] == "R":
    FC_STRAND = 2

elif config["strandedness"] == "F":
    FC_STRAND = 1

else:
    FC_STRAND = 0


## Get extra parameters for gene calling

if FORMAT == "PE":
    FC_GENES_PARAMS = " -R CORE -g gene_id -t transcript -p --countReadPairs"

else:
    FC_GENES_PARAMS = " -R CORE -g gene_id -t transcript"


## Get extra parameters for exon calling

if FORMAT == "PE":
    FC_EXONS_PARAMS = " -R CORE -g gene_id -J -p --countReadPairs"

else:
    FC_EXONS_PARAMS = " -R CORE -g gene_id -J"


### TRACK NORMALIZATION HELPERS

# What is input for normalization?
if config.get("use_exons_only", True):
    NORMALIZATION_INPUT = expand(
        "results/featurecounts_exons/{sample}.featureCounts", sample=SAMP_NAMES
    )

else:
    NORMALIZATION_INPUT = expand(
        "results/featurecounts_genes/{sample}.featureCounts", sample=SAMP_NAMES
    )

# Are exons being used?
if config.get("use_exons_only", True):
    NORMALIZATION_EXTRA = "--use_exons_only"

else:
    NORMALIZATION_EXTRA = ""
