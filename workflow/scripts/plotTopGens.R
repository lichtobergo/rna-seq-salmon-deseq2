log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

#_______________Loading packages______________________________#
suppressPackageStartupMessages(library(tidyverse))
library(readr)
library(DESeq2)
#library(fgsea)



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
  snakemake@input[[1]], 
  delim = "\t", 
  escape_double = FALSE,
  trim_ws = TRUE
)



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
      dplyr::rename(
        ensembl_gene_id = gene,
        gene = external_gene_name
      ) %>% 
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
        group %in% snakemake@config[["diffexp"]][["contrasts"]][["M.blood-vs-R.blood"]]
      ) %>% 
      ggplot(aes(x = sample.group, y = count)) +
      geom_boxplot() +
      scale_y_continuous(labels = scales::label_comma()) +
      facet_wrap(vars(gene), scales = "free") +
      labs(title = paste0("Top 10 ", idx, " genes")) +
      theme_minimal()
  )








