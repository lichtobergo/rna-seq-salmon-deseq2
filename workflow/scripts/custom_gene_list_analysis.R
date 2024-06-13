library(DESeq2)
library(tidyverse)
# load data
setwd("~/mnt/projects_share/p1400")
dds <- readRDS("~/mnt/projects_share/p1400/results/deseq2/dds.RDS")

custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish)

custom_colors$lukas <- c('#ED7D31',
                         '#70AD47',
                         '#4472C4',
                         '#AEABAB')

# read xlsx file with sheets for each gene category
multipleSheets <- function(fname) {
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tbl <- purrr::map(sheets, \(x) readxl::read_xlsx(fname, sheet = x))
  data_frame <- purrr::map(tbl, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # print data frame
  print(data_frame)
}

path <- "resources/gene_lists_moritz.xlsx"
multipleSheets(path)

gene_lists <- multipleSheets(path)

# get the samples of interest from the data

samples <- c("08A",
             "09A",
             "11A",
             "08B",
             "09B",
             "11B",
             "01A",
             "03A",
             "07A")

dds <- dds[ , samples]


# get normalized count data for each gene and gene group

countData <- gene_lists %>% 
  purrr::imap(
    \(x, idx) x %>% 
      pull(ensembl, gene)
  ) %>% 
  purrr::imap(
    \(y, idy) y %>% 
      purrr::imap(
      \(z, idz)
      plotCounts(
        dds = dds,
        gene = z,
        intgroup = "group",
        returnData = TRUE
      )
    ) %>% 
      list_rbind(names_to = "gene") %>% 
      as_tibble()
  )

# make plots
plots <- countData %>% 
  imap(
    \(x, idx) x %>% 
      ggplot(
        aes(
          x = gene,
          y = count,
          fill = group
        )
      ) +
      # labs(title = str_replace_all(idx, "_", " ")) +
      geom_bar(position = "dodge", stat = "identity") +
      geom_point(size = 1.1, position = position_jitterdodge()) +
      scale_fill_manual(values = c('#AEABAB',
                                   '#ED7D31',
                                   '#70AD47')) +
      # scale_y_log10() +
      theme_classic() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  )

# save plots to disk
plots %>% 
  imap(
    \(x, idx)  
      ggsave(
        plot = x,
        device = ragg::agg_png(),
        path = "results/plots/custom_gene_counts", 
        filename = paste0(idx, ".png"),
        width = 12,
        height = 12*10/16,
        units = "cm",
        dpi = 300,
        scale = 1.5
        )
  )

# save raw count tables to disk
writexl::write_xlsx(countData,
                    path = "results/custom_gene_counts/custom_gene_tabel.xlsx")
