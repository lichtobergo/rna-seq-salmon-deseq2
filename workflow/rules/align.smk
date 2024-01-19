rule salmon_quant:
    input:
        unpack(get_fq),
        index = config["index"]["location"]
    output:
        "results/pseudoaligment/quants/{sample}/quant.sf"
    params:
        dir = "results/pseudoaligment/quants/{sample}"
    benchmark:
        "workflow/benchmarks/salmon_quant/{sample}.tsv"
    priority: 50
    threads: 10
    conda: "../envs/salmon.yaml"
    shell:
        """
        salmon quant -i {input.index} -l A -p {threads} --gcBias --validateMappings -o {params.dir} -1 {input.fq1} -2 {input.fq2}
        """
