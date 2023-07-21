if STAR:
    rule transcript_fn:
        input:
            rsem="results/rsem_csv/{sample}_rsem.csv.gz",
            counts="results/counts/{sample}_counts.csv.gz",
        output:
            outfile="results/transcript_fn/{sample}_RSEM_plus.csv",
        params:
            rscript = workflow.source_path("../scripts/RSEM_plus.R"),
            pnew = get_pnew,
            pold = get_pold,
        log:
            "logs/transcript_fn/{sample}.log"
        threads: 20
        conda:
            "../envs/full.yaml"
        shell:
            r"""
            chmod +x {params.rscript}
            {params.rscript} -o {output.outfile} -c {input.counts} -r {input.rsem} -s {wildcards.s4U_sample} -n {params.pnew} -b {params.pold} 1> {log} 2>&1
            """

    rule combine_fn:
        input:
            expand("results/transcript_fn/{samps}_RSEM_plus.csv", samps = SAMP_NAMES),
        output:
            "results/transcript_fn/RSEM_plus.csv"
        log:
            "logs/combine_fn/combine_fn.log"
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            """
            cat {input} >> {output}
            """