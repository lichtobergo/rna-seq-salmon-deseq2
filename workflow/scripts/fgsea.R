log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

#_______________Loading packages______________________________#
suppressPackageStartupMessages(library(tidyverse))
library(readr)
library(msigdbr)
library(fgsea)



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



#_______________Assign variables_______________________________#
rankingMetric <- snakemake@config[["enrichment"]][["fgsea"]][["ranking_metric"]]
species <- str_to_sentence(
  str_replace(snakemake@config[["ref"]][["species"]])
)
GOterm <- snakemake@config[["enrichment"]][["fgsea"]][["GO"]][["GO_domain"]]
# subcategory <- snakemake@config[["enrichment"]][["fgsea"]][["MSigDB"]][["subcategory"]] 
sigLevel <- snakemake@config[["enrichment"]][["fgsea"]][["fdr_gene_set"]]



#_________________GSEA___________________#
# Steps toward doing gene set enrichment analysis (GSEA):

# 1- obtaining stats for ranking genes in your experiment,
# 2- creating a named vector out of the DESeq2 result
# 3- Obtaining a gene set from mysigbd
# 4- doing analysis
# run fgsea algorithm


# remove the NAs, averaging statistics for a multi-hit symbol
if (rankingMetric == "log2FoldChange") {
  res2 <- res %>% 
    dplyr::select(
      log2FoldChange,
      gene
    ) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(gene) %>% 
    summarize(stat = mean(log2FoldChange))
} else {
  res2 <- res %>% 
    mutate(
      metric = sign(log2FoldChange) * -log10(pvalue)
    ) %>% 
    dplyr::select(
      metric,
      gene
    ) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(gene) %>% 
    summarize(stat = mean(metric))
}

# creating  a named vector [ranked genes]
ranks <- deframe(res2)
head(ranks, 20)

# Load the pathway (gene set) into a named list
# download MsigDb gene set TODO: make it conditional 
GOset <- msigdbr(
  species = species, 
  category = "C5", 
  subcategory = "GOterm"
)
geneSet <- split(
  x = GOset$gene_symbol,
  f = GOset$gs_name
)

#Running fgsea algorithm:
fgseaRes <- fgsea(
  pathways = geneSet,
  stats = ranks
)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(1:10)

# To see what genes are in each of these pathways:
gene.in.pathway <- geneSet %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest(cols = c(SYMBOL)) %>% 
  inner_join(res, by=c("SYMBOL" = "gene"))



##______________________VISUALIZATION______________________________#

#__________bar plot _______________#
# Plot the normalized enrichment scores. 
# Color the bar indicating whether or not the pathway was significant:
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
plot <- fgseaResTidy %>% 
  mutate(
    pathway = str_remove(pathway, "GOBP_") %>% 
      str_to_sentence() %>% 
      str_replace_all("_", " ")
  ) %>% 
  ggplot(aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO:Biological processes pathways NES from GSEA",
       subtitle = "Top 10 up- and down regulated pathways") + 
  theme_minimal() +
  theme(
    plot.background = element_rect()
  )

ggsave(
  filename = snakemake@output[[1]],
  plot = plot,
  dpi = 300
)

