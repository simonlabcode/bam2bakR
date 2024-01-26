
if config["bam2bakr"]:
    if config["remove_tags"]:
        rule remove_tags:
            input:
                input_bam=get_input_bams,
            output:
                output_bam="results/remove_tags/{sample}_no_jI_jM.bam",
            log:
                "logs/remove_tags/{sample}.log"
            conda:
                "../envs/full.yaml"
            script:
                "../scripts/remove_tags.py"

        # Filter out multi-mappers and sort reads
        rule sort_filter:
            input:
                "results/remove_tags/{sample}_no_jI_jM.bam",
            output:
                "results/sf_reads/{sample}.s.sam",
                "results/sf_reads/{sample}_fixed_mate.bam",
                "results/sf_reads/{sample}.f.sam",
            log:
                "logs/sort_filter/{sample}.log"
            params: 
                shellscript=workflow.source_path("../scripts/sort_filter.sh")
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} {config[FORMAT]} 1> {log} 2>&1
                """

    else:
        # Filter out multi-mappers and sort reads
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
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} {config[FORMAT]} 1> {log} 2>&1
                """


else:
    # Filter out multi-mappers and sort reads
    rule sort_filter:
        input:
            "results/bams/{sample}Aligned.out.bam"
        output:
            "results/sf_reads/{sample}.s.sam",
            "results/sf_reads/{sample}_fixed_mate.bam",
            "results/sf_reads/{sample}.f.sam",
        log:
            "logs/sort_filter/{sample}.log"
        params: 
            shellscript=workflow.source_path("../scripts/sort_filter.sh"),
            read_format=config["FORMAT"],
        threads: 8
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.read_format} 1> {log} 2>&1
            """



# Use custom htseq script to quanity features 
# Also creates bam files with tag designating feature that each read was mapped to; useful during mutation counting
rule htseq_cnt:
    input:
        "results/sf_reads/{sample}.s.sam"
    output:
        "results/htseq/{sample}_tl.bam",
        temp("results/htseq/{sample}_check.txt")
    params: 
        shellscript=workflow.source_path("../scripts/htseq.sh"),
        pythonscript=workflow.source_path("../scripts/count_triple.py"),
        strand=config["strandedness"],
        flattened=config["flattened"],
        annotation=config["annotation"],
    log:
        "logs/htseq_cnt/{sample}.log"
    threads: 3
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.annotation} {params.strand} {params.pythonscript} {params.flattened} 1> {log} 2>&1
        """

# Calculate normalization scale factor to be applied to tracks        
if NORMALIZE:
    rule normalize:
        input:
            expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
        output:
            "results/normalization/scale"
        log:
            "logs/normalize/normalize.log"
        params:
            rscript=workflow.source_path("../scripts/normalize.R"),
            spikename=config["spikename"],
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            r"""
            chmod +x {params.rscript}
            {params.rscript} --dirs ./results/htseq/ --spikename {params.spikename}
            mv scale {output}
            """
else:
    rule normalize:
        input:
            expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
        output:
            "results/normalization/scale"
        log:
            "logs/normalize/normalize.log"
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            """
            touch {output}
            """

# Index genome fasta file for snp calling
rule index:
    input:
        str(config["genome_fasta"])
    output:
        get_index_name()
    log:
        "logs/genome-faidx.log",
    threads: 1
    conda:
        "../envs/index.yaml"
    script:
        "../scripts/genome-faidx.py"

# Identify SNPs to be accounted for when counting mutations
rule call_snps:
    input:
        str(config["genome_fasta"]),
        get_index_name(),
        expand("results/htseq/{ctl}_tl.bam", ctl = CTL_NAMES)
    params:
        nctl = nctl,
        shellscript = workflow.source_path("../scripts/call_snps.sh"),
    output:
        "results/snps/snp.txt",
        "results/snps/snp.vcf",
        temp("results/snps/mkdir.txt")
    log:
        "logs/call_snps/ctl_samps.log"
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {params.nctl} {output} {input} 1> {log} 2>&1
        """

# TO-DO:
# 1) Add mutation position optimizations and functionality
# Count mutations 
rule cnt_muts:
    input:
        "results/sf_reads/{sample}.s.bam",
        "results/snps/snp.txt"
    params:
        format = FORMAT,
        minqual = config["minqual"],
        mut_tracks = config["mut_tracks"],
        strand = config["strandedness"],
        shellscript = workflow.source_path("../scripts/mut_call.sh"),
        pythonscript = workflow.source_path("../scripts/mut_call.py"),
        awkscript = workflow.source_path("../scripts/fragment_sam.awk"),
        mutpos = False
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt")
    log:
        "logs/cnt_muts/{sample}.log"
    threads: 32
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.minqual} {params.mut_tracks} {params.format} {params.strand} {params.pythonscript} {params.awkscript} {params.mutpos} 1> {log} 2>&1
        """


# Merge mutation counts with feature assignment
rule merge_features_and_muts:
    input:
        get_merge_input,
    output:
        "results/merge_features_and_muts/{sample}_counts.csv.gz"
    params:
        genes_included = True,
        exons_included = True,
        exonbins_included = False,
        transcripts_included = False,
        rscript = workflow.source_path("../scripts/merge_features_and_muts.R")
    log:
        "logs/merge_features_and_muts/{sample}.log"
    threads: 20
    conda: 
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.rscript}

        {params.rscript} -g {params.genes_included} -e {params.exons_included} -b {params.exonbins_included} \
        -t {params.transcripts_included} -o {output} -s {wildcards.sample} 1> {log} 2>&1
        """


# Make color-coded tracks
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
        pythonscript = workflow.source_path("../scripts/count_to_tracks.py"),
        mut_tracks=config["mut_tracks"],
        genome_fasta=config["genome_fasta"],
        wsl=config["WSL"],
        normalize=config["normalize"],
    log:
        "logs/maketdf/{sample}.log"
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {params.mut_tracks} {params.genome_fasta} {params.wsl} {params.normalize} {params.pythonscript} {output} 1> {log} 2>&1
        """

# Make cB file that will be input to bakR
rule makecB:
    input:
        expand("results/merge_features_and_muts/{sample}_counts.csv.gz", sample=config["samples"])
    output:
        "results/cB/cB.csv.gz"
    params:
        shellscript = workflow.source_path("../scripts/master.sh"),
        keepcols = config["keepcols"],
        mut_tracks = config["mut_tracks"],
        mut_pos = False,
    log:
        "logs/makecB/master.log"
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {output.cB} {params.keepcols} {params.mut_tracks} \
        ./results/merge_features_and_muts/ {params.mut_pos} 1> {log} 2>&1
        """