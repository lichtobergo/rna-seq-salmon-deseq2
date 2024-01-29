rule pca2pdf:
    output:
        "results/plots/pca_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/pca/*.png {output}
        """
rule topGenes2pdf:
    output:
        "results/plots/topGenes_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/topGenes/*.png {output}
        """

rule volcano2pdf:
    output:
        "results/plots/volcano_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/volcano/*.png {output}
        """

rule gseaGO2pdf:
    output:
        "results/plots/gseaGO_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/gseaGO/*.png {output}
        """

rule gseaMSigDB2pdf:
    output:
        "results/plots/gseaMSigDB_plots.pdf"
    conda:
        "../envs/magick.yaml"
    shell:
        """
        convert results/plots/gseaMSigDB/*.png {output}
        """