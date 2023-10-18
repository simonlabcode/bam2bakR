## Create STAR index if not available
if config['build_star']:
    rule star_index:
        input:
            fasta=config["genome_fasta"],
            annotation=config["annotation"],
        output:
            directory(config['STAR_index']),
        threads: 4
        params:
            extra="--sjdbGTFfile {}".format(str(config["annotation"])),
        log:
            "logs/star_index_genome.log",
        wrapper:
            "v1.25.0/bio/star/index"


## Paired-end reads
if FORMAT == 'PE':

    # Run cutadapt
    rule cutadapt:
        input:
            get_input_fastqs
        output:
            fastq1="results/fastq_cut/{sample}.t.r1.fastq",
            fastq2="results/fastq_cut/{sample}.t.r2.fastq",
            qc = "results/fastq_cut/{sample}.qc.txt",
        params:
            adapters=config["adapter"],
            extra=config["cutadapt_extra"], 
        log:
            "logs/cutadapt/{sample}.log",
        threads: 8  # set desired number of threads here
        wrapper:
            "v1.25.0/bio/cutadapt/pe"       

    # Run hisat-3n
    if ALIGNER:
        rule align:
            input:
                "results/fastq_cut/{sample}.t.r1.fastq",
                "results/fastq_cut/{sample}.t.r2.fastq"
            output:
                "results/bams/{sample}Aligned.out.bam",
            log:
                "logs/align/{sample}.log"
            params:
                shellscript = workflow.source_path("../scripts/hisat_3n.sh"),
                format = config["FORMAT"],
                strand = config["strandedness"],
                chr = config["chr_tag"],
                h3n = config["HISAT_3N"],
                h3n_path = config["hisat3n_path"],
                muts = config["mut_tracks"],
                yale = config["Yale"]
            threads: 24
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h3n} {params.h3n_path} {params.muts} {params.yale} {input} {output} 1> {log} 2>&1
                """

    # Run STAR
    elif STAR:

        rule align:
            input:
                fq1 = "results/fastq_cut/{sample}.t.r1.fastq",
                fq2 = "results/fastq_cut/{sample}.t.r2.fastq",
                index = config['STAR_index'],
            output:
                aln="results/bams/{sample}Aligned.out.bam",
                reads_per_gene="results/bams/{sample}-ReadsPerGene.out.tab",
                aln_tx="results/bams/{sample}-Aligned.toTranscriptome.out.bam",
            log:
                "logs/bams/{sample}.log",
            params:
                idx=lambda wc, input: input.index,
                extra="--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile {} {}".format(
                    str(config["annotation"]), config["star_extra"]
                ),
            conda:
                "../envs/star.yaml"
            threads: 24
            script: 
                "../scripts/star-align.py"

        rule RSEM_index:
            input:
                reference_genome=config['genome_fasta'],
            output:
                seq="index/reference.seq",
                grp="index/reference.grp",
                ti="index/reference.ti",
                tfa="index/reference.transcripts.fa",
                idxfa="index/reference.idx.fa",
                n2g="index/reference.n2g.idx.fa",
            params:
                extra="--gtf {}".format(str(config["annotation"])),
            log:
                "logs/rsem/prepare-reference.log",
            threads: 20
            wrapper:
                "v2.3.1/bio/rsem/prepare-reference"
                
        rule RSEM:
            input:
                bam="results/bams/{sample}-Aligned.toTranscriptome.out.bam",
                reference=multiext("index/reference", ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa"),
            output:
                genes_results="results/rsem/{sample}.genes.results",
                isoforms_results="results/rsem/{sample}.isoforms.results",
                bam="results/rsem/{sample}.transcript.bam"
            params:
                # optional, specify if sequencing is paired-end
                paired_end= True,
                # additional optional parameters to pass to rsem, for example,
            log:
                "logs/rsem/calculate_expression/{sample}.log",
            conda:
                "../envs/rsem.yaml"
            threads: 20
            script:
                "../scripts/rsem-calc.py"

        rule rsem_to_csv:
            input:
                "results/rsem/{sample}.transcript.bam",
            output:
                "results/rsem_csv/{sample}_rsem.csv.gz",
                temp("results/rsem_csv/{sample}_check.txt"),
            params:
                shellscript = workflow.source_path("../scripts/rsem_to_csv.sh"),
                pythonscript = workflow.source_path("../scripts/rsem_csv.py"),
                awkscript = workflow.source_path("../scripts/fragment_sam_rsem.awk")
            log:
                "logs/rsem_to_csv/{sample}.log"
            threads: 20
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                chmod +x {params.pythonscript}
                chmod +x {params.awkscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.pythonscript} {params.awkscript} 1> {log} 2>&1
                """
    # Run hisat2
    else:
        ## Old way of running hisat2 with custom script
       # rule align:
           # input:
           #     "results/fastq_cut/{sample}.t.r1.fastq",
           #     "results/fastq_cut/{sample}.t.r2.fastq"
           # output:
           #     "results/bams/{sample}Aligned.out.bam",
           # log:
           #     "logs/align/{sample}.log"
           # params:
           #     shellscript = workflow.source_path("../scripts/hisat2.sh"),
           #     format = config["FORMAT"],
           #     strand = config["strandedness"],
           #     chr = config["chr_tag"],
           #     h2 = config["HISAT2"]
           # threads: 20
           # conda:
           #     "../envs/full.yaml"
           # shell:
           #     """
           #     chmod +x {params.shellscript}
           #     {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h2} {input} {output} 1> {log} 2>&1
           #     """
        
        # Running hisat2 with Snakemake wrapper
        rule align:
            input:
                reads=["results/fastq_cut/{sample}.t.r1.fastq","results/fastq_cut/{sample}.t.r2.fastq"],
                idx=config["HISAT2"],
            output:
                "results/bams/{sample}Aligned.out.bam",
            log:
                "logs/hisat2_align/{sample}.log",
            params:
                extra=config["hisat2_extra"],
            threads: 20
            wrapper:
                "v1.25.0/bio/hisat2/align"

## Single end data
else:

    # Run cutadapt
    rule cutadapt:
        input:
            get_input_fastqs
        output:
            fastq="results/fastq_cut/{sample}.t.fastq",
            qc = "results/fastq_cut/{sample}.qc.txt",
        params:
            adapters=config["adapter"],
            extra=config["cutadapt_extra"], 
        log:
            "logs/cutadapt/{sample}.log",
        threads: 8  # set desired number of threads here
        wrapper:
            "v1.25.0/bio/cutadapt/se" 

    # Run hisat-3n
    if ALIGNER:
        rule align:
            input:
                "results/fastq_cut/{sample}.t.fastq",
            output:
                "results/bams/{sample}Aligned.out.bam",
            log:
                "logs/align/{sample}.log"
            params:
                shellscript = workflow.source_path("../scripts/hisat_3n.sh"),
                format = config["FORMAT"],
                strand = config["strandedness"],
                chr = config["chr_tag"],
                h3n = config["HISAT_3N"],
                h3n_path = config["hisat3n_path"],
                muts = config["mut_tracks"],
                yale = config["Yale"]
            threads: 20
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h3n} {params.h3n_path} {params.muts} {params.yale} {input} {output} 1> {log} 2>&1
                """

    # Run STAR
    elif STAR:
        rule align:
            input:
                fq1 = "results/fastq_cut/{sample}.t.fastq",
                index = config['STAR_index'],
            output:
                aln="results/bams/{sample}Aligned.out.bam",
                reads_per_gene="results/bams/{sample}-ReadsPerGene.out.tab",
                aln_tx="results/bams/{sample}-Aligned.toTranscriptome.out.bam",
            log:
                "logs/bams/{sample}.log",
            params:
                idx=lambda wc, input: input.index,
                extra="--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile {} {}".format(
                    str(config["annotation"]), config["star_extra"]
                ),
            conda:
                "../envs/star.yaml"
            threads: 24
            script: 
                "../scripts/star-align.py"

        rule RSEM_index:
            input:
                reference_genome=config['genome_fasta'],
            output:
                seq="index/reference.seq",
                grp="index/reference.grp",
                ti="index/reference.ti",
                tfa="index/reference.transcripts.fa",
                idxfa="index/reference.idx.fa",
                n2g="index/reference.n2g.idx.fa",
            params:
                extra="--gtf {}".format(str(config["annotation"])),
            log:
                "logs/rsem/prepare-reference.log",
            threads: 20
            wrapper:
                "v1.23.4/bio/rsem/prepare-reference"

        rule RSEM:
            input:
                bam="results/bams/{sample}-Aligned.toTranscriptome.out.bam",
                reference=multiext("index/reference", ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa"),
            output:
                genes_results="results/rsem/{sample}.genes.results",
                isoforms_results="results/rsem/{sample}.isoforms.results",
            params:
                # optional, specify if sequencing is paired-end
                paired_end= False,
                # additional optional parameters to pass to rsem, for example,
            log:
                "logs/rsem/calculate_expression/{sample}.log",
            conda:
                "../envs/rsem.yaml"
            threads: 20
            script:
                "../scripts/rsem-calc.py"

        rule rsem_to_csv:
            input:
                "results/rsem/{sample}.transcript.bam",
            output:
                "results/rsem_csv/{sample}_rsem.csv.gz",
                temp("results/rsem_csv/{sample}_check.txt"),
            params:
                shellscript = workflow.source_path("../scripts/rsem_to_csv.sh"),
                pythonscript = workflow.source_path("../scripts/rsem_csv.py"),
                awkscript = workflow.source_path("../scripts/fragment_sam.awk")
            log:
                "logs/rsem_to_csv/{sample}.log"
            threads: 20
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                chmod +x {params.pythonscript}
                chmod +x {params.awkscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.pythonscript} {params.awkscript} 1> {log} 2>&1
                """

    # Run hisat2
    else:
        ## Old way of running hisat2 with custom script
        #rule align:
            #input:
            #    "results/fastq_cut/{sample}.t.fastq",
            #output:
            #    "results/bams/{sample}Aligned.out.bam",
            #log:
            #    "logs/align/{sample}.log"
            #params:
            #    shellscript = workflow.source_path("../scripts/hisat2.sh"),
            #    format = config["FORMAT"],
            #    strand = config["strandedness"],
            #    chr = config["chr_tag"],
            #    h2 = config["HISAT2"]
            #threads: 20
            #conda:
            #    "../envs/full.yaml"
            #shell:
            #    """
            #    chmod +x {params.shellscript}
            #    {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h2} {input} {output} 1> {log} 2>&1
            #    """
        
        # Run hisat2 with Snakemake wrapper
        rule align:
            input:
                reads=["results/fastq_cut/{sample}.t.fastq"],
                idx=config["HISAT2"],
            output:
                "results/bams/{sample}Aligned.out.bam",
            log:
                "logs/hisat2_align/{sample}.log",
            params:
                extra=config["hisat2_extra"],
            threads: 20
            wrapper:
                "v1.25.0/bio/hisat2/align"
