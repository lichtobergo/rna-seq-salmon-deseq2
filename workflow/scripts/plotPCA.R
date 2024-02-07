log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Define custom PCA plotting function
plotPCAcustom = function(object, intgroup="condition",
                         ntop=500, returnData=FALSE, pcsToUse=1:2)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  pcs <- paste0("PC", pcsToUse)
  d <- data.frame(V1=pca$x[,pcsToUse[1]],
                  V2=pca$x[,pcsToUse[2]],
                  group=group, intgroup.df, name=colnames(object))
  colnames(d)[1:2] <- pcs
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcsToUse]
    return(d)
  }
  
  ggplot(data=d, aes_string(x=pcs[1], y=pcs[2], color="group")) +
    geom_point(size=3) + 
    xlab(paste0(pcs[1],": ",round(percentVar[pcsToUse[1]] * 100),"% variance")) +
    ylab(paste0(pcs[2],": ",round(percentVar[pcsToUse[2]] * 100),"% variance")) +
    coord_fixed()
}


library(DESeq2)
suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
library(ggtext)

# load DESeq2 data
dds <- readRDS(snakemake@input[[1]])
# dds <- DESeq(dds)


rld <- rlog(
  object = dds,
  blind = FALSE
)
pca_12 <- plotPCAcustom(
  object = rld,
  
  intgroup = snakemake@wildcards[["variable"]],
  pcsToUse = 1:2,
  returnData = TRUE
)

pca_34 <- plotPCAcustom(
  object = rld,
  
  intgroup = snakemake@wildcards[["variable"]],
  pcsToUse = 3:4,
  returnData = TRUE
)

#save.image("pcaDebug.RData")

covariate <- snakemake@wildcards[["variable"]]
percent_var <- round(100 * attr(pca_12, "percentVar"))
if (length(levels(colData(rld)[[covariate]])) > 8) {
   color_pal <- RColorBrewer::brewer.pal(12, "Paired")
} else {
  color_pal <- RColorBrewer::brewer.pal(8, "Dark2")
}

plot_12<- ggplot(
  data = pca_12,
  aes_string(
    x = "PC1",
    y = "PC2",
    color = snakemake@wildcards[["variable"]]
    # label = "name"
  )
) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(
    aes(label = name) , 
    size = 3, #
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_color_manual(values = color_pal) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 8) +
  # ggtitle(
  #   label = str_glue("Principal components 1 and 2 color coded for covariate **{covariate}**")
  #   ) +
  theme(
    text = element_text(colour = "grey30"),
    axis.text = element_text(colour = "grey30"),
    line = element_line(colour = "grey30"),
    panel.grid = element_blank(),
    legend.justification = "top",
    plot.subtitle = element_textbox_simple(),
    plot.title = element_textbox_simple()
  ) 
  
  # ggtitle(
  #   label = NULL,
  #   subtitle = bquote(
  #     "Principal components 1 and 2 color coded for covariate" 
  #     ~ bold(.(snakemake@wildcards[["variable"]]))
  #   )
  # )
percent_var <- round(100 * attr(pca_34, "percentVar"))

plot_34 <- ggplot(
  data = pca_34,
  aes_string(
    x = "PC3",
    y = "PC4",
    color = snakemake@wildcards[["variable"]]
  )
) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(
    aes(label = name) , 
    size = 3, #
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_color_manual(values = color_pal) +
  xlab(paste0("PC3: ", percent_var[1], "% variance")) +
  ylab(paste0("PC4: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 8) + 
  # ggtitle(
  #   label = str_glue("Principal components 3 and 4 color coded for covariate **einseheshehehehehehehehehehehehehehehehehehehehr langer string**")
  # ) +
  theme(
    text = element_text(colour = "grey30"),
    axis.text = element_text(colour = "grey30"),
    line = element_line(colour = "grey30"),
    panel.grid = element_blank(),
    legend.justification = "top",
    plot.title = element_textbox_simple()
  )


plot <- plot_12 + plot_34 + plot_layout(guides = "collect") +
  plot_annotation(
    title = 'Principal Component Analysis',
    subtitle = str_glue("The first 4 principal components of the data are plotted and samples are color coded for the covariate **{covariate}**"),
    caption = 'Each dot represents a sample and is labelled with the sample name.',
    theme = theme(text = element_text(size = 8), plot.subtitle = element_textbox_simple())
  )
ggsave(plot = plot, filename = snakemake@output[[1]], width=20, height = 11, units = "cm")

# png(snakemake@output[[1]])
# print(plot)
# dev.off()

#save.image(file = "debugR.RData")
