### THESE RULES PERTAIN TO THE ASSIGNMENT OF READS TO FEATURES WITH FEATURECOUNTS


# Assign reads to genes
rule featurecounts_genes:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_genes/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=FC_GENES_PARAMS,
    log:
        "logs/featurecounts_gene/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


# Assign reads to exons
rule featurecounts_exons:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_exons/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=FC_EXONS_PARAMS,
    log:
        "logs/featurecounts_exons/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"
