rule tximport:
    input:
        quant = expand("results/pseudoaligment/quants/{sample}/quant.sf", sample=samples["sample.name"]),
        tx_to_gene = "resources/tx2gene.tsv"
    output:
        txi = "results/deseq2/txi.RDS"
    params:
        extra = "type='salmon', txOut=FALSE"
    log: "workflow/logs/tximport/tximport.log"
    #benchmark:
    #    "workflow/benchmarks/tximport/{sample}.tsv"
    threads: 1
    wrapper:
        "v3.3.3/bio/tximport"

rule DESeqDataSet_from_tximport:
    input:
        txi="results/deseq2/txi.RDS",
        colData=config["samples"],
    output:
        temp("results/deseq2/dds_txi.RDS"),
    threads: 1
    log:
        "workflow/logs/DESeqDataSet/txi.log",
    params:
        formula=config["diffexp"]["model"],  # Required R statistical formula
        # factor="condition", # Optionally used for relevel
        # reference_level="A", # Optionally used for relevel
        # tested_level="B", # Optionally used for relevel
        # min_counts=0, # Optionally used to filter low counts
        # extra="", # Optional parameters provided to import function
    wrapper:
        "v3.3.3/bio/deseq2/deseqdataset"
        
rule gene_2_symbol:
    input:
        counts="{prefix}.tsv",
    output:
        symbol="{prefix}.symbol.tsv",
    params:
        species=get_bioc_species_name(),
    log:
        "workflow/logs/gene2symbol/{prefix}.log",
    conda:
        "../envs/biomart.yaml"
    script:
        "../scripts/gene2symbol.R"

rule deseq2_init:
    input:
        "results/deseq2/dds_txi.RDS"
    output:
        "results/deseq2/dds.RDS",
        "results/deseq2/normcounts.tsv",
    params:
        get_bioc_species_name(),
        get_count_threshold()
    threads: get_deseq2_threads()
    log:
        "workflow/logs/deseq2_init/init.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2-init.R"

rule plot_PCA:
    input:
        "results/deseq2/dds.RDS"
    output:
        # "results/plots/pca.png",
        report("results/plots/pca.{variable}.png", "../report/pca.rst"),
    params:
        pca_labels = config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "workflow/logs/plots/pca/pca.{variable}.log"
    script:
        "../scripts/plotPCA.R"

rule diffexp:
    input:
        "results/deseq2/dds.RDS",
    output:
        report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
    params:
        contrast = get_contrast
    log:
        "workflow/logs/deseq2/{contrast}.diffexp.log",
    conda:
        "../envs/deseq2.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2.R"

rule Enhanced_Volcano:
    input:
        "results/diffexp/{contrast}.diffexp.symbol.tsv"
    output:
        # "results/plots/pca.png",
        report("results/plots/volcano.{contrast}.png", "../report/volcano.rst"),
    params:
        contrast = get_contrast
    conda:
        "../envs/EnhancedVolcano.yaml"
    log:
        "workflow/logs/plots/volcano/volcano.{contrast}.log"
    script:
        "../scripts/EnhancedVolcano.R"

rule fgsea_GO:
    input:
        "results/diffexp/{contrast}.diffexp.symbol.tsv"
    output:
        # "results/plots/pca.png",
        report("results/plots/gseaGO.{contrast}.png", "../report/gsea.rst"),
    params:
        contrast = get_contrast
    conda:
        "../envs/fgsea.yaml"
    log:
        "workflow/logs/plots/fgsea/fgseaGO.{contrast}.log"
    script:
        "../scripts/fgseaGO.R"

rule fgsea_MSigDB:
    input:
        "results/diffexp/{contrast}.diffexp.symbol.tsv"
    output:
        # "results/plots/pca.png",
        report("results/plots/gseaMSigDB.{contrast}.png", "../report/gsea.rst"),
    params:
        contrast = get_contrast
    conda:
        "../envs/fgsea.yaml"
    log:
        "workflow/logs/plots/fgsea/fgseaMSigDB.{contrast}.log"
    script:
        "../scripts/fgseaMSigDB.R"

rule plotTopGenes:
    input:
        res = "results/diffexp/{contrast}.diffexp.symbol.tsv",
        dds = "results/deseq2/dds.RDS",
    output:
        upregulated = report("results/plots/top10Up.{contrast}.png", "../report/topGenes.rst"),
        downregulated = report("results/plots/top10Down.{contrast}.png", "../report/topGenes.rst"),
    params:
        contrast = get_contrast
    conda:
        "../envs/topGenes.yaml"
    log:
        "workflow/logs/plots/plotTopGenes/topGenes.{contrast}.log"
    script:
        "../scripts/plotTopGenes.R"