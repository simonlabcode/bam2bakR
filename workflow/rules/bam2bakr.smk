rule sort_filter:
    input:
        get_input_bams
    output:
        "results/sf_reads/{sample}.s.sam",
        "results/sf_reads/{sample}_fixed_mate.bam",
        "results/sf_reads/{sample}.f.sam",
    log:
        "logs/sort_filter/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/sort.yaml"
    shell:
        "workflow/scripts/sort_filter.sh {threads} {wildcards.sample} {input} {output} {config[FORMAT]} 1> {log} 2>&1"

rule htseq_cnt:
    input:
        "results/sf_reads/{sample}.s.sam"
    output:
        "results/htseq/{sample}_tl.bam",
        temp("results/htseq/{sample}_check.txt")
    log:
        "logs/htseq_cnt/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/htseq.yaml"
    shell:
        "workflow/scripts/htseq.sh {threads} {wildcards.sample} {input} {output} {config[annotation]} {config[mutcnt]} 1> {log} 2>&1"

rule normalize:
    input:
        expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
    output:
        "results/normalization/scale"
    log:
        "logs/normalize/normalize.log"
    threads: 1
    conda:
        "../envs/normalize.yaml"
    shell:
        r"""
       	./workflow/scripts/normalize.R --dirs ./results/htseq/ --spikename {config[spikename]}
	mv scale {output}
        """

rule index:
    input:
        str(config["genome_fasta"])
    output:
        get_index_name()
    log:
        "logs/genome-faidx.log",
    wrapper:
        "v1.21.4/bio/samtools/faidx"


rule call_snps:
    input:
        get_index_name(),
        expand("results/htseq/{ctl}_tl.bam", ctl = CTL_NAMES)
    params:
        nsamps = nctl
    output:
        "results/snps/snp.txt",
        "results/snps/snp.vcf",
        temp("results/snps/mkdir.txt")
    log:
        "logs/call_snps/ctl_samps.log"
    threads: workflow.cores
    conda:
        "../envs/snps.yaml"
    shell:
        "workflow/scripts/call_snps.sh {threads} {params.nsamps} {output} {config[genome_fasta]} {input} 1> {log} 2>&1"

rule cnt_muts:
    input:
        "results/htseq/{sample}_tl.bam",
        "results/snps/snp.txt"
    params:
        format = config["FORMAT"],
        fragment_size = config["fragment_size"],
        minqual = config["minqual"],
        mut_tracks = config["mut_tracks"],
        strand = config["strandedness"]
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt")
    log:
        "logs/cnt_muts/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "./workflow/scripts/mut_call.sh {threads} {wildcards.sample} {input} {output} {params.fragment_size} {params.minqual} {params.mut_tracks} {params.format} {params.strand} 1> {log} 2>&1"

rule maketdf:
    input:
        "results/counts/{sample}_counts.csv.gz",
        "results/htseq/{sample}_tl.bam",
	    "results/normalization/scale"
    output:
        temp("results/tracks/{sample}_success.txt"),
        expand("results/tracks/{{sample}}.{mut}.{id}.{strand}.tdf", mut=config["mut_tracks"], id=[0,1,2,3,4,5], strand = ['pos', 'min'])
    log:
        "logs/maketdf/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/tracks.yaml"
    shell:
        "workflow/scripts/tracks.sh {threads} {wildcards.sample} {input} {config[mut_tracks]} {config[genome_fasta]} {config[WSL]} {config[normalize]} {output} 1> {log} 2>&1"

rule makecB:
    input:
        expand("results/counts/{sample}_counts.csv.gz", sample=config["samples"])
    output:
        "results/cB/cB.csv.gz"
    log:
        "logs/makecB/master.log"
    threads: workflow.cores
    conda:
        "../envs/cB.yaml"
    shell:
        "workflow/scripts/master.sh {threads} {output} {config[keepcols]} {config[mut_tracks]} 1> {log} 2>&1"