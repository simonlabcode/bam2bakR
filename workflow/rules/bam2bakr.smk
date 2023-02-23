CPU_num = config['cpus']

SAMP_NAMES = list(config['samples'].keys())

CTL_NAMES = list(config['control_samples'])

nctl = len(CTL_NAMES)

def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


rule sort_filter:
    input:
        get_input_bams
    output:
        "results/sf_reads/{sample}.s.sam",
        "results/sf_reads/{sample}_fixed_mate.bam",
        "results/sf_reads/{sample}.f.sam",
    log:
        "logs/sort_filter/{sample}.log"
    threads: CPU_num
    conda:
        "../envs/sort.yaml"
    shell:
        "workflow/scripts/sort_filter.sh {config[cpus]} {wildcards.sample} {input} {output} {config[FORMAT]}"

rule htseq_cnt:
    input:
        "results/sf_reads/{sample}.s.sam"
    output:
        "results/htseq/{sample}_tl.bam",
        temp("results/htseq/{sample}_check.txt")
    log:
        "logs/htseq_cnt/{sample}.log"
    threads: CPU_num
    conda:
        "../envs/htseq.yaml"
    shell:
        "workflow/scripts/htseq.sh {config[cpus]} {wildcards.sample} {input} {output} {config[annotation]} {config[mutcnt]}"

rule normalize:
    input:
        expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
    output:
        "results/normalization/scale"
    log:
        "logs/normalize/normalize.log"
    threads: CPU_num
    conda:
        "../envs/normalize.yaml"
    shell:
        r"""
       	./workflow/scripts/normalize.R --dirs ./results/htseq/ --spikename {config[spikename]}
	mv scale {output}
        """


rule call_snps:
    input:
        expand("results/htseq/{ctl}_tl.bam", ctl = CTL_NAMES)
    params:
        nsamps = nctl
    output:
        "results/snps/snp.txt",
        temp("results/snps/mkdir.txt")
    log:
        "logs/call_snps/ctl_samps.log"
    threads: CPU_num
    conda:
        "../envs/snps.yaml"
    shell:
        "workflow/scripts/call_snps.sh {config[cpus]} {params.nsamps} {output} {config[genome_fasta]} {input}"

rule cnt_muts:
    input:
        "results/htseq/{sample}_tl.bam",
        "results/snps/snp.txt"
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt")
    log:
        "logs/cnt_muts/{sample}.log"
    threads: CPU_num
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "workflow/scripts/mut_call.sh {config[cpus]} {wildcards.sample} {input} {output} {config[awkscript]} {config[fragment_size]} {config[minqual]} {config[mut_tracks]} {config[mutcall]} {config[FORMAT]}"

rule maketdf:
    input:
        "results/counts/{sample}_counts.csv.gz",
        "results/htseq/{sample}_tl.bam",
	    "results/normalization/scale"
    output:
        temp("results/tracks/{sample}_success.txt"),
        expand("results/tracks/{{sample}}.{mut}.{id}.{strand}.tdf", mut=config["mut_tracks"], id=[0,1,2,3,4,5], strand = ['pos', 'min'])
    threads: CPU_num
    conda:
        "../envs/tracks.yaml"
    shell:
        "workflow/scripts/tracks.sh {config[cpus]} {wildcards.sample} {input} {config[mut_tracks]} {config[genome_fasta]} {config[WSL]} {config[normalize]} {output}"

rule makecB:
    input:
        expand("results/counts/{sample}_counts.csv.gz", sample=config["samples"])
    output:
        "results/cB/cB.csv.gz"
    log:
        "logs/makecB/master.log"
    threads: CPU_num
    conda:
        "../envs/cB.yaml"
    shell:
        "workflow/scripts/master.sh {config[cpus]} {output} {config[keepcols]} {config[mut_tracks]}"
