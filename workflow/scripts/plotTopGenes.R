log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

#_______________Loading packages______________________________#
suppressPackageStartupMessages(library(tidyverse))
library(readr)
library(DESeq2)



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
      ggplot(aes(x = sample.group, y = count, color = sample.group)) +
      geom_boxplot(alpha = 0.1, show.legend = FALSE) +
      geom_jitter(size = 3, position = position_jitter(0.2)) +
      scale_y_continuous(labels = scales::label_comma()) +
      scale_color_brewer(palette = "Dark2") +
      facet_wrap(
        vars(gene), 
        nrow = 2,
        scales = "free",
        labeller = label_wrap_gen(width = 10, multi_line = TRUE)
        ) +
      labs(
        title = paste0("Top 10 ", idx, " genes"),
        y = "Normalized counts"
        ) +
      theme_minimal() +
      theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.background = element_rect(),
        panel.background = element_rect(color = "grey30"),
        strip.text = element_text(face = "bold"),
        legend.position = "top"
      )
  )

"grey20"

# save.image(file = "debugR.RData")

ggsave(
  filename = snakemake@output[["upregulated"]],
  plot = topGenes$upgregulated,
  dpi = 300
)

ggsave(
  filename = snakemake@output[["downregulated"]],
  plot = topGenes$downregulated,
  dpi = 300
)





