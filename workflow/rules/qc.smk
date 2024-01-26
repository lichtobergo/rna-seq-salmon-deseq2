rule fastqc_raw:
    input:
        unpack(fq_dict_from_sample)
    output:
        "results/qc/fastqc/{sample}_R1_fastqc.html",
        "results/qc/fastqc/{sample}_R1_fastqc.zip",
        "results/qc/fastqc/{sample}_R2_fastqc.html",
        "results/qc/fastqc/{sample}_R2_fastqc.zip",
    benchmark:
        "workflow/benchmarks/fastqc/{sample}.tsv"
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc {input.fq1}  --outdir results/qc/fastqc/ --quiet &&
        fastqc {input.fq2}  --outdir results/qc/fastqc/ --quiet
        """

rule fastqc_trimmed:
    input:
        unpack(get_fq)
    output:
        "results/qc/fastqc/trimmed/{sample}_R1_trimmed_fastqc.html",
        "results/qc/fastqc/trimmed/{sample}_R1_trimmed_fastqc.zip",
        "results/qc/fastqc/trimmed/{sample}_R2_trimmed_fastqc.html",
        "results/qc/fastqc/trimmed/{sample}_R2_trimmed_fastqc.zip",
    benchmark:
        "workflow/benchmarks/fastqc/{sample}.tsv"
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc {input.fq1}  --outdir results/qc/fastqc/trimmed/ --quiet &&
        fastqc {input.fq2}  --outdir results/qc/fastqc/trimmed/ --quiet
        """

rule multiqc:
    input:
        expand(["results/qc/fastqc/{sample}_{pair}_fastqc.html",
                "results/qc/fastqc/trimmed/{sample}_{pair}_trimmed_fastqc.html",
                "data/trimmed/{sample}_{pair}.trimming_report.txt",
                "results/pseudoaligment/quants/{sample}/quant.sf"],
                sample=samples["sample.name"], pair=["R1", "R2"]),
    output:
        "results/qc/multiqc/multiqc_report.html",
        directory("results/qc/multiqc_data"),
    params:
        extra="-ip -v",
    log:
        "workflow/logs/qc/multiqc.log"
    benchmark:
        "workflow/benchmarks/multiqc/multiqc.tsv"
    wrapper:
        "v3.3.3/bio/multiqc"