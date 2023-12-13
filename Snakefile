# An example of sample ung MultiQC in a snakemake workflow.
import pandas as pd

configfile: "00_config/config.yaml"
print("Config is: ", config)
samples = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample.name", drop = False)
    .sort_index()
)

# samples = samples.iloc[:2]
READS = ["R1"",R2"]
print(samples)

# wildcard_constraints:
#     sample="|".join(samples["sample_name"])

def fq_dict_from_sample(wildcards):
    return {
        "fq1": samples.loc[wildcards.sample, "read1"],
        "fq2": samples.loc[wildcards.sample, "read2"]
    }

def get_fastqs(wildcards):
    return dict(
        zip(
            ["fq1", "fq2"],
            expand(
                "06_results/trimmed/{sample}_{pair}.fastq.gz",
                pair = ["R1", "R2"],
                **wildcards,
            ),
        )
    )
#print(sample)
rule all:
    input:
        expand("04_QC/Fastqc/{sample}_{reads}_fastqc.{ext}", 
                sample=samples["sample.name"], 
                reads=["R1", "R2"], 
                ext=["html", "zip"]),
         "04_QC/Multiqc/multiqc_report.html",
        
rule fastqc_file:
    input:
        unpack(fq_dict_from_sample)
    output:
        "04_QC/Fastqc/{sample}_R1_fastqc.html",
        "04_QC/Fastqc/{sample}_R1_fastqc.zip",
        "04_QC/Fastqc/{sample}_R2_fastqc.html",
        "04_QC/Fastqc/{sample}_R2_fastqc.zip"
    benchmark:
        "05_benchmarks/fastqc/{sample}.tsv"
    threads: 2
    conda: "fastqc"
    shell:
        """
        mkdir -p 04_QC/Fastqc
        fastqc {input.fq1} {input.fq2} -t {threads} --outdir 04_QC/Fastqc
        """

rule salmon_quant:
    input:
        unpack(get_fastqs),
        index = config["index"]["location"]
    output:
        "06_results/pseudoaligment/quants/{sample}/quant.sf"
    params:
        dir = "06_results/pseudoaligment/quants/{sample}"
    benchmark:
        "05_benchmarks/salmon_quant/{sample}.tsv"
    threads: 10
    conda:
        "salmon"
    shell:
        """
        salmon quant -i {input.index} -l A -p {threads} --gcBias --validateMappings -o {params.dir} -1 {input.fq1} -2 {input.fq2}
        """


rule multiqc:
    input:
        expand(["04_QC/Fastqc/{sample}_{reads}_fastqc.html",
                "06_results/pseudoaligment/quants/{sample}/quant.sf"],
                sample=samples["sample.name"], reads=["R1", "R2"]),
    output:
        "04_QC/Multiqc/multiqc_report.html",
        directory("04_QC/multiqc_data"),
    params:
        extra="-ip -v",
    benchmark:
        "05_benchmarks/multiqc/multiqc.tsv"
    wrapper:
        "v3.0.2/bio/multiqc"

rule trim_galore:
    input:
        unpack(fq_dict_from_sample)
    output:
        fasta_fwd = "06_results/trimmed/{sample}_R1.fastq.gz",
        report_fwd = "06_results/trimmed/{sample}_R1.trimming_report.txt",
        fasta_rev = "06_results/trimmed/{sample}_R2.fastq.gz",
        report_rev = "06_results/trimmed/{sample}_R2.trimming_report.txt",
    benchmark:
        "05_benchmarks/trim_galore/{sample}.tsv"
    threads: 4
    wrapper:
        "v3.0.3/bio/trim_galore/pe"
        
rule tximport:
    input:
        quant = expand("06_results/pseudoaligment/quants/{sample}/quant.sf", sample=samples["sample.name"]),
        tx_to_gene = "03_resources/tx2gene.tsv"
    output:
        txi = "06_results/txi.RDS"
    params:
        extra = "type='salmon', txOut=FALSE"
    #benchmark:
    #    "05_benchmarks/tximport/{sample}.tsv"
    threads: 1
    wrapper:
        "v3.1.0/bio/tximport"

rule test_DESeqDataSet_from_tximport:
    input:
        txi="06_results/txi.RDS",
        colData=config["samples"],
    output:
        "06_results/dds_txi.RDS",
    threads: 1
    log:
        "logs/DESeqDataSet/txi.log",
    params:
        formula="~group",  # Required R statistical formula
        # factor="condition", # Optionally used for relevel
        # reference_level="A", # Optionally used for relevel
        # tested_level="B", # Optionally used for relevel
        # min_counts=0, # Optionally used to filter low counts
        # extra="", # Optional parameters provided to import function
    wrapper:
        "v3.1.0/bio/deseq2/deseqdataset"

rule plot_PCA:
    input:
        "06_results/dds_txi.RDS"
    output:
        "06_results/plots/pca.png"
    params:
        pca_labels = config["pca"]["labels"]
    #conda: 
    #    "envs/pca.yaml"
    log:
        "logs/plots/pca/plot_pca.log"
    script:
        "02_scripts/plotPCA.R"

