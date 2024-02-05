# rule pca2pdf:
#     input:
#         expand()
#     output:
#         "results/plots/pca_plots.pdf"
#     conda:
#         "../envs/magick.yaml"
#     shell:
#         """
#         convert results/plots/pca/*.png {output}
#         """

rule topGenes2pdf:
    input:
        expand("results/plots/topGenes/topGenes.{contrast}.png", contrast=config["diffexp"]["contrasts"])
    output:
        "results/plots/topGenes_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/topGenes/*.png {output}
        """

rule volcano2pdf:
    input:
        expand("results/plots/volcano/volcano.{contrast}.png", contrast=config["diffexp"]["contrasts"])
    output:
        "results/plots/volcano_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/volcano/*.png {output}
        """

rule gseaGO2pdf:
    input:
        expand("results/plots/gseaGO/gseaGO.{contrast}.png", contrast=config["diffexp"]["contrasts"])
    output:
        "results/plots/gseaGO_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/gseaGO/*.png {output}
        """

rule gseaMSigDB2pdf:
    input:
        expand("results/plots/gseaMSigDB/gseaMSigDB.{contrast}.png", contrast=config["diffexp"]["contrasts"])
    output:
        "results/plots/gseaMSigDB_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/gseaMSigDB/*.png {output}
        """