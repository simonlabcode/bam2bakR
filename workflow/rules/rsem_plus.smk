if STAR:
    rule transcript_fn:
        input:
            rsem="results/rsem_csv/{s4U_sample}_rsem.csv.gz",
            counts="results/counts/{s4U_sample}_counts.csv.gz",
        output:
            outfile="results/transcript_fn/{s4U_sample}_RSEM_plus.csv",
        params:
            rscript = workflow.source_path("../scripts/RSEM_plus.R"),
        log:
            "logs/transcript_fn/{s4U_sample}.log"
        threads: workflow.cores
        conda:
            "../envs/full.yaml"
        shell:
            r"""
            chmod +x {params.rscript}
            {params.rscript} -o {output.outfile} -c {input.counts} -r {input.rsem} -s {wildcards.s4U_sample} 1> {log} 2>&1
            """

    rule combine_fn:
        input:
            expand("results/transcript_fn/{samps}_RSEM_plus.csv", samps = s4U_SAMPS),
        output:
            "results/transcript_fn/RSEM_plus.csv"
        log:
            "logs/transcript_fn/combine_fn.log"
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            """
            cat {input} >> {output}
            """