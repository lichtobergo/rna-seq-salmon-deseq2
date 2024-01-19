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


percent_var <- round(100 * attr(pca_12, "percentVar"))

plot_12<- ggplot(
  data = pca_12,
  aes_string(
    x = "PC1",
    y = "PC2",
    color = snakemake@wildcards[["variable"]],
    label = "name"
  )
) +
  geom_point(size = 2.5) +
  ggrepel::geom_label_repel(size = 3, max.overlaps = 20) +
  #scale_color_manual(values = custom_colors$discrete) +
  # scale_color_manual(values = c("#66c2a5", "#fc8d62")) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    line = element_line(colour = "black"),
    panel.grid = element_blank(),
    legend.justification = "top",
    plot.title = element_text(hjust = 0.5)
  ) + 
  ggtitle("PCA-plot")

percent_var <- round(100 * attr(pca_34, "percentVar"))

plot_34 <- ggplot(
  data = pca_34,
  aes_string(
    x = "PC3",
    y = "PC4",
    color = snakemake@wildcards[["variable"]],
    label = "name"
  )
) +
  geom_point(size = 2.5) +
  ggrepel::geom_label_repel(size = 3, max.overlaps = 20) +
  #scale_color_manual(values = custom_colors$discrete) +
  # scale_color_manual(values = c("#66c2a5", "#fc8d62")) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    line = element_line(colour = "black"),
    panel.grid = element_blank(),
    legend.justification = "top",
    plot.title = element_text(hjust = 0.5)
  ) + 
  ggtitle("PCA-plot")

plot <- plot_12 + plot_34 + plot_layout(guides = "collect")
ggsave(plot = plot, filename = snakemake@output[[1]], width=12, height = 6)

# png(snakemake@output[[1]])
# print(plot)
# dev.off()

save.image(file = "debugR.RData")