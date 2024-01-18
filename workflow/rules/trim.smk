rule trim_galore:
    input:
        unpack(fq_dict_from_sample)
    output:
        fasta_fwd = "results/trimmed/{sample}_R1_trimmed.fastq.gz",
        report_fwd = "results/trimmed/{sample}_R1.trimming_report.txt",
        fasta_rev = "results/trimmed/{sample}_R2_trimmed.fastq.gz",
        report_rev = "results/trimmed/{sample}_R2.trimming_report.txt",
    benchmark:
        "workflow/benchmarks/trim_galore/{sample}.tsv"
    priority: 100
    threads: 4
    wrapper:
        "v3.3.3/bio/trim_galore/pe"
        