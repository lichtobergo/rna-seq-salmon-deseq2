log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(biomaRt)
suppressPackageStartupMessages(library(tidyverse))
# useful error messages upon aborting
library("cli")

# this data set was aligned to the Ensembl mRatBN.7.2 assemly version 110
# ensembl <- useMart("ensembl")
# # listDatasets(ensembl)
# ensembl <- useDataset("rnorvegicus_gene_ensembl", mart = ensembl)
# attributes <- listAttributes(ensembl)
# attributes[1:10, ]
# tx2gene <- getBM(
#   attributes = c(
#     "ensembl_transcript_id_version",
#     "ensembl_gene_id",
#     "external_gene_name"
#   ),
#   mart = ensembl
# )

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because 
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while ( class(mart)[[1]] != "Mart" ) {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        # dataset = "rnorvegicus_gene_ensembl",
        dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
        version = snakemake@config[["ref"]][["release"]],
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        cli_abort(
          str_c(
            "Have tried all 4 available Ensembl biomaRt mirrors ",
            rounds,
            " times. You might have a connection problem, or no mirror is responsive.\n",
            "The last error message was:\n",
            message(e)
          )
        )
      }
      # hop to next mirror
      mart <- switch(mart,
                     useast = "uswest",
                     uswest = "asia",
                     asia = "www",
                     www = {
                       # wait before starting another round through the mirrors,
                       # hoping that intermittent problems disappear
                       Sys.sleep(30)
                       "useast"
                     }
      )
    }
  )
}


df <- read.table(snakemake@input[["counts"]], sep='\t', header=1)

g2g <- biomaRt::getBM(
  attributes = c( "ensembl_gene_id",
                  "external_gene_name"),
  filters = "ensembl_gene_id",
  values = df$gene,
  mart = mart
)
g2g$external_gene_name <- ifelse(g2g$external_gene_name == "", g2g$ensembl_gene_id, g2g$external_gene_name)

annotated <- g2g %>% 
  right_join(
    df,
    by = c("ensembl_gene_id" = "gene")
  ) %>%
  dplyr::rename(gene = external_gene_name) %>% 
  mutate(
    # gene = if_else(str_length(gene) == 0, ensembl_gene_id, gene)
    gene = if_else(str_length(gene) == 0, ensembl_gene_id, gene),
    gene = if_else(is.na(gene), ensembl_gene_id, gene)
  )

write.table(annotated, snakemake@output[["symbol"]], sep='\t', row.names=F)
