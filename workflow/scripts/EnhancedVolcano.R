log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(svglite)
# load results table 
res <- readr::read_delim(
  file = snakemake@input[[1]], 
  delim = "\t", escape_double = FALSE, 
  trim_ws = TRUE
)

# results table
# name vector for purrr mapping
# as.data.frame() %>% # convert to data frame
#   rownames_to_column(var = "ensembl_gene_id") %>% # make rownames a column
#   # join with tx2gene table for external gene names
#   left_join(distinct(
#     tx2gene[ , c(2,3)],
#     ensembl_gene_id,
#     .keep_all = TRUE
#   ))


# generate volcano plots
plot <- EnhancedVolcano::EnhancedVolcano(
  toptable = res,
  lab = res$gene,
  x = "log2FoldChange",
  y = "padj",
  ylab = bquote(~-Log[10] ~ "adjusted" ~ italic(P)),
  pCutoff = snakemake@config[["diffexp"]][["sig-level"]][["volcano-plot"]],
  pCutoffCol = "padj",
  title = paste(snakemake@params[[1]][2], " vs. ", snakemake@params[[1]][3])
)
# save plot as png file
ggplot2::ggsave(
  plot = plot,
  filename = snakemake@output[[1]],
  # width = 21.3,
  height = 18,
  units = "cm",
  dpi = 300
)

#save.image(file = "debugR.RData")