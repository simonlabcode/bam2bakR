rule sort_filter:
    input:
        get_input_bams
    output:
        "results/sf_reads/{sample}.s.sam",
        "results/sf_reads/{sample}_fixed_mate.bam",
        "results/sf_reads/{sample}.f.sam",
    log:
        "logs/sort_filter/{sample}.log"
    params: 
        shellscript=workflow.source_path("../scripts/sort_filter.sh")
    threads: workflow.cores
    conda:
        "../envs/sort.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {config[FORMAT]} 1> {log} 2>&1
        """

rule htseq_cnt:
    input:
        "results/sf_reads/{sample}.s.sam"
    output:
        "results/htseq/{sample}_tl.bam",
        temp("results/htseq/{sample}_check.txt")
    params: 
        shellscript=workflow.source_path("../scripts/htseq.sh"),
        pythonscript=workflow.source_path("../scripts/count_triple.py"),
        strand=lambda wildcards: "yes" if config["strandedness"] == 'F' else: "reverse"
    log:
        "logs/htseq_cnt/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/htseq.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {config[annotation]} {params.strand} {params.pythonscript} 1> {log} 2>&1
        """

rule normalize:
    input:
        expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
    output:
        "results/normalization/scale"
    log:
        "logs/normalize/normalize.log"
    params:
        rscript=workflow.source_path("../scripts/normalize.R")
    threads: 1
    conda:
        "../envs/normalize.yaml"
    shell:
        r"""
        chmod +x {params.rscript}
       	{params.rscript} --dirs ./results/htseq/ --spikename {config[spikename]}
	    mv scale {output}
        """

rule index:
    input:
        str(config["genome_fasta"])
    output:
        get_index_name()
    log:
        "logs/genome-faidx.log",
    conda:
        "../envs/faidx.yaml"
    script:
        "../scripts/genome-faidx.py"


rule call_snps:
    input:
        get_index_name(),
        expand("results/htseq/{ctl}_tl.bam", ctl = CTL_NAMES)
    params:
        nsamps = nctl,
        shellscript = workflow.source_path("../scripts/call_snps.sh")
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
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {params.nsamps} {output} {config[genome_fasta]} {input} 1> {log} 2>&1
        """

rule cnt_muts:
    input:
        "results/htseq/{sample}_tl.bam",
        "results/snps/snp.txt"
    params:
        format = config["FORMAT"],
        minqual = config["minqual"],
        mut_tracks = config["mut_tracks"],
        strand = config["strandedness"],
        shellscript = workflow.source_path("../scripts/mut_call.sh"),
        pythonscript = workflow.source_path("../scripts/mut_call.py"),
        awkscript = workflow.source_path("../scripts/fragment_sam.awk")
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt")
    log:
        "logs/cnt_muts/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.minqual} {params.mut_tracks} {params.format} {params.strand} {params.pythonscript} {params.awkscript} 1> {log} 2>&1
        """

rule maketdf:
    input:
        "results/counts/{sample}_counts.csv.gz",
        "results/htseq/{sample}_tl.bam",
	    "results/normalization/scale"
    output:
        temp("results/tracks/{sample}_success.txt"),
        expand("results/tracks/{{sample}}.{mut}.{id}.{strand}.tdf", mut=config["mut_tracks"], id=[0,1,2,3,4,5], strand = ['pos', 'min'])
    params:
        shellscript = workflow.source_path("../scripts/tracks.sh"),
        pythonscript = workflow.source_path("../scripts/count_to_tracks.py")
    log:
        "logs/maketdf/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/tracks.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {config[mut_tracks]} {config[genome_fasta]} {config[WSL]} {config[normalize]} {params.pythonscript} {output} 1> {log} 2>&1
        """

rule makecB:
    input:
        expand("results/counts/{sample}_counts.csv.gz", sample=config["samples"])
    output:
        "results/cB/cB.csv.gz"
    params:
        shellscript = workflow.source_path("../scripts/master.sh")
    log:
        "logs/makecB/master.log"
    threads: workflow.cores
    conda:
        "../envs/cB.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {output} {config[keepcols]} {config[mut_tracks]} 1> {log} 2>&1
        """