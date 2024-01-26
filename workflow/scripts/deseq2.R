log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("cli")
library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  # setup parallelization
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

contrast_config <- snakemake@config[["diffexp"]][["contrasts"]][[
  snakemake@wildcards[["contrast"]]
]]

# basic case of contrast specification, see:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
contrast <- c(
  contrast_config[1],
  contrast_config[2],
  contrast_config[3]
)

  # more complex contrast specification via list(c(), c()), see ?results docs of
  # the DESeq2 package and this tutorial (plus the linked seqanswers thread):
  # https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md

res <- results(
  dds,
  contrast = contrast,
  alpha = snakemake@config[["diffexp"]][["sig-level"]][["diffexp"]],
  parallel = parallel
)
# shrink fold changes for lowly expressed genes
# use ashr so we can use `contrast` as conversion to coef is not trivial, see:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
res <- lfcShrink(
  dds,
  contrast = contrast,
  res = res,
  type = "ashr"
)

# sort by p-value
res <- res[order(res$padj), ]
# TODO explore IHW usage


# store results
# svg(snakemake@output[["ma_plot"]])
# plotMA(res, ylim = c(-2, 2))
# dev.off()

write.table(
  data.frame(
    "gene" = rownames(res),
    res
  ),
  file = snakemake@output[[1]],
  row.names = FALSE,
  sep = "\t"
)
