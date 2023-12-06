# An example of sample ung MultiQC in a snakemake workflow.
import pandas as pd

configfile: "config/config.yaml"
print("Config is: ", config)
samples = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample.name", drop = False)
    .sort_index()
)

samples = samples.iloc[:2]
READS = ["R1"",R2"]
print(samples)

# wildcard_constraints:
#     sample="|".join(samples["sample_name"])

def fq_dict_from_sample(wildcards):
    return {
        "fq1": samples.loc[wildcards.sample, "read1"],
        "fq2": samples.loc[wildcards.sample, "read2"]
    }


#print(sample)
rule all:
    input:
        expand("01_QC/Fastqc/{sample}_{reads}_fastqc.{ext}", 
                sample=samples["sample.name"], 
                reads=["R1", "R2"], 
                ext=["html", "zip"]),
         "01_QC/Multiqc/multiqc_report.html",
        
rule fastqc_file:
    input:
        unpack(fq_dict_from_sample)
    output:
        "01_QC/Fastqc/{sample}_R1_fastqc.html",
        "01_QC/Fastqc/{sample}_R1_fastqc.zip",
        "01_QC/Fastqc/{sample}_R2_fastqc.html",
        "01_QC/Fastqc/{sample}_R2_fastqc.zip"
    benchmark:
        "benchmarks/fastqc/{sample}.tsv"
    threads: 1
    conda: "fastqc"
    shell:
        """
        mkdir -p 01_QC/Fastqc
        fastqc {input.fq1} {input.fq2} -t {threads} --outdir 01_QC/Fastqc
        """

rule salmon_quant:
    input:
        unpack(fq_dict_from_sample),
        index = config["index"]["location"]
    output:
        "data/quants/{sample}/quant.sf"
    params:
        dir = "data/quants/{sample}"
    benchmark:
        "benchmarks/salmon_quant/{sample}.tsv"
    threads: 4
    conda:
        "salmon"
    shell:
        """
        salmon quant -i {input.index} -l A -p {threads} --gcBias --validateMappings -o {params.dir} -1 {input.fq1} -2 {input.fq2}
        """


rule multiqc:
    input:
        expand(["01_QC/Fastqc/{sample}_{reads}_fastqc.html",
                "data/quants/{sample}/quant.sf"],
                sample=samples["sample.name"], reads=["R1", "R2"]),
    output:
        "01_QC/Multiqc/multiqc_report.html",
        directory("01_QC/multiqc_data"),
    params:
        extra="-ip -v",
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    wrapper:
        "v3.0.2/bio/multiqc"
    


