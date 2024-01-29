log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

#_______________Loading packages______________________________#
suppressPackageStartupMessages(library(tidyverse))
library(readr)
library(DESeq2)
library(patchwork)
library(ggtext)



#_________________Set up parallel computing___________________#
parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  # setup parallelization
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}



#_______________Load input data_______________________________# 
res <- read_delim(
  snakemake@input[["res"]], 
  delim = "\t", 
  escape_double = FALSE,
  trim_ws = TRUE
)
snakemake@input[["res"]]

dds <- readRDS(snakemake@input[["dds"]])



##______________________VISUALIZATION______________________________#

# Tidy and annotate the data
splitTbl <- res %>% 
  dplyr::filter(
    !is.na(padj),
    pvalue != 1
  ) %>%
  group_by(sign(log2FoldChange)) %>% 
  group_split() %>% 
  purrr::map(
    \(x) x %>% 
      arrange(padj) %>% 
      slice_head(n = 10)
  ) %>% 
  set_names(c("downregulated", "upgregulated"))

topGenes <- splitTbl %>% 
  purrr::imap(
    \(x, idx) x %>% 
      pull(ensembl_gene_id, gene) %>% 
      purrr::imap(
        \(y, idy) plotCounts(
          dds = dds,
          gene = y,
          intgroup = c(names(snakemake@config[["diffexp"]][["variables_of_interest"]]), "group"),
          returnData = TRUE
        ) 
      ) %>% 
      list_rbind(names_to = "gene") %>% 
      dplyr::filter(
        group %in% snakemake@config[["diffexp"]][["contrasts"]][[
          snakemake@wildcards[["contrast"]]
          ]]
      ) %>% 
      ggplot(aes(x = sample.group, y = count, color = group)) +
      geom_boxplot(alpha = 0.1, show.legend = FALSE) +
      geom_jitter(size = 2, position = position_jitter(0.2)) +
      scale_y_continuous(labels = scales::label_comma()) +
      scale_color_brewer(palette = "Dark2") +
      facet_wrap(
        vars(gene), 
        nrow = 2,
        scales = "free",
        labeller = label_wrap_gen(width = 10, multi_line = TRUE)
        ) +
      labs(
        # title =paste(snakemake@params[[1]][2], " vs. ", snakemake@params[[1]][3]),
        title = paste0("Top 10 ", idx, " genes"),
        y = "Normalized counts"
        ) +
      theme_minimal(base_size = 8) +
      theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title.position = "plot",
        plot.background = element_blank(),
        panel.background = element_rect(color = "grey30"),
        strip.text = element_text(face = "bold"),
        legend.position = "top"
      )
  )

"grey20"


# ggsave(
#   filename = snakemake@output[["upregulated"]],
#   plot = topGenes$upgregulated,
#   dpi = 300
# )
# 
# ggsave(
#   filename = snakemake@output[["downregulated"]],
#   plot = topGenes$downregulated,
#   dpi = 300
# )

plot <- (topGenes$upgregulated + theme(plot.margin = unit(c(0,20,0,0), "pt"))) +
  (topGenes$downregulated + theme(plot.margin = unit(c(0,0,0,20), "pt"))) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste(snakemake@params[[1]][2], " vs. ", snakemake@params[[1]][3]),
    # subtitle = str_glue("The first 4 principal components of the data are plotted and samples are color coded for the covariate **{covariate}**"),
    # caption = 'Each dot represents a sample and is labelled with the sample name.',
    theme = theme(text = element_text(size = 10), plot.subtitle = element_textbox_simple())
  ) &
  theme(legend.position = "bottom")

ggsave(plot = plot, filename = snakemake@output[["topGenes"]], width=20, height = 11, units = "cm")





