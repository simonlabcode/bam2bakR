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


### Functions for returning input to first rules


# Relevant if providing bam files as input
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


# Relevant if providing fastq files as input
def get_input_fastqs(wildcards):
    fastq_path = config["samples"][wildcards.sample]
    fastq_files = sorted(glob.glob(f"{fastq_path}/*.fastq*"))
    return fastq_files


### FeatureCounts parameters


## All of the files to merge
def get_merge_input(wildcards):
    MERGE_INPUT = []

    MERGE_INPUT.extend(
        expand("results/counts/{SID}_counts.csv.gz", SID=wildcards.sample)
    )

    MERGE_INPUT.extend(
        expand("results/featurecounts_genes/{SID}.featureCounts", SID=wildcards.sample)
    )

    MERGE_INPUT.extend(
        expand("results/featurecounts_exons/{SID}.featureCounts", SID=wildcards.sample)
    )

    return MERGE_INPUT


# Columns in final cB
keepcols = ["sample", "sj", "rname"]


keepcols.append("GF")


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
    FC_GENES_PARAMS = " -R CORE -g gene_id -t transcript -p --countReadPairs  --nonOverlap 0 --primary"

else:
    FC_GENES_PARAMS = " -R CORE -g gene_id -t transcript  --nonOverlap 0 --primary"


## Get extra parameters for exon calling

if FORMAT == "PE":
    FC_EXONS_PARAMS = (
        " -R CORE -g gene_id -J -p --countReadPairs --nonOverlap 0 --primary"
    )

else:
    FC_EXONS_PARAMS = " -R CORE -g gene_id -J  --nonOverlap 0 --primary"
