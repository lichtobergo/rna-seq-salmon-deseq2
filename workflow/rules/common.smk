import pandas as pd
import itertools

#### load config and sample sheets #####

samples = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample.name", drop = False)
    .sort_index()
)

units = (
    pd.read_csv(config["units"], sep="\t")
    .set_index("sample.name", drop = False)
    .sort_index()
)

simpleContrasts={}

def getSimpleContrasts(columns,rows):
    simpleContrasts={}
    i=0
    # simple Contrasts e.g. treated vs untreated
    for key in columns: # n*(n-1)/2 combinations
        for values in itertools.combinations(rows[i],2):
            if values[0]!=values[1]:
                name=values[0]+"-vs-"+values[1]
                simpleContrasts[name]=[key,values[0],values[1]]
        i+=1
    return simpleContrasts

def estimateContrasts():
    # columns=list(samples.keys())[1::] # get column names without overhead
    columns=list(["group"])
    rows=[]
    for key in columns:
        rows.append(list(samples.get(key).unique())) # get distinct rows values

    formula=config["diffexp"]["model"]
    if config["diffexp"]["generate-contrasts"]:
        global simpleContrasts
        simpleContrasts=getSimpleContrasts(columns,rows)
        print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::::::::: Trying to create all useful contrasts ::::::::::::::::::::::::::")
        print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
        config["diffexp"]["contrasts"]={**simpleContrasts}
        print(config["diffexp"]["contrasts"])
    else:
        print("::::::::::::::::::::::::::::::::::::::::::::::")
        print("::: Your custom formula design was choosen :::")
        print("::::::::::::::::::::::::::::::::::::::::::::::")

#samples = samples.iloc[:2]
#READS = ["R1"",R2"]
print(samples)

if config["diffexp"]["generate-contrasts"]:
    estimateContrasts()


def get_final_output():
    final_output = expand(
        "results/diffexp/{contrast}.diffexp.symbol.tsv",
        contrast=config["diffexp"]["contrasts"],
    )
    final_output.append("results/deseq2/normcounts.symbol.tsv")
    # final_output.append("results/counts/all.symbol.tsv")
    final_output.append("results/qc/multiqc/multiqc_report.html")

    if config["pca"]["activate"]:
        # final_output.append("results/plots/pca.png")
        # get all the variables to plot a PCA for
        pca_variables = list(config["diffexp"]["variables_of_interest"])
        if config["diffexp"]["batch_effects"]:
            pca_variables.extend(config["diffexp"]["batch_effects"])
        if config["pca"]["labels"]:
            pca_variables.extend(config["pca"]["labels"])
        final_output.extend(
            expand("results/plots/pca/pca.{variable}.png", variable=pca_variables)
        )
    final_output.append(
        expand("results/plots/volcano/volcano.{contrast}.png", contrast=config["diffexp"]["contrasts"])
    )
    # if config["enrichment"]["fgsea"]["activate"]:
    if config["enrichment"]["fgsea"]["GO"]["activate"]:
        final_output.append(
            expand("results/plots/gseaGO/gseaGO.{contrast}.png", contrast=config["diffexp"]["contrasts"])
        )
    if config["enrichment"]["fgsea"]["MSigDB"]["activate"]:
        final_output.append(
            expand("results/plots/gseaMSigDB/gseaMSigDB.{contrast}.png", contrast=config["diffexp"]["contrasts"])
        )
    final_output.append(
        expand("results/plots/topGenes/topGenes.{contrast}.png", contrast=config["diffexp"]["contrasts"])
    )
    final_output.append(
        expand(
            "results/plots/{type}_plots.pdf", type=["topGenes", "volcano", "gseaGO", "gseaMSigDB"]
        )
    )
    # final_output.append(
    #     "results/plots/volcano_plots.pdf"
    # )
    # final_output.append(
    #     "results/plots/topGenes_plots.pdf"
    # )
    return final_output      

def fq_dict_from_sample(wildcards):
    return {
        "fq1": units.loc[wildcards.sample, "read1"],
        "fq2": units.loc[wildcards.sample, "read2"]
    }

def get_fastqs(wildcards):
    return dict(
        zip(
            ["fq1", "fq2"],
            expand(
                "data/trimmed/{sample}_{pair}_trimmed.fastq.gz",
                pair = ["R1", "R2"],
                **wildcards,
            ),
        )
    )

def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return {
        "fq1": units.loc[wildcards.sample, "read1"],
        "fq2": units.loc[wildcards.sample, "read2"]
    }
    else:
        return dict(
            zip(
                ["fq1", "fq2"],
                expand(
                    "data/trimmed/{sample}_{pair}_trimmed.fastq.gz",
                    pair = ["R1", "R2"],
                    **wildcards
                )
            )
        )

def get_bioc_species_name():
    first_letter = config["ref"]["species"][0]
    subspecies = config["ref"]["species"].split("_")[1]
    return first_letter + subspecies

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]

def get_count_threshold():
    if config["diffexp"]["prefilter"]["custom"] == None:
        return config["diffexp"]["prefilter"]["threshold"]
    else:
        return config["diffexp"]["prefilter"]["custom"]
