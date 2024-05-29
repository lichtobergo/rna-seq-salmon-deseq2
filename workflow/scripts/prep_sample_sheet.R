library(tidyverse)
library(readxl)

sample_sheet <- read_delim(
  "config/metadata.tsv", 
  delim = "\t", escape_double = FALSE, 
  col_types = cols(µl.of.RNA = col_skip(),
                   files = col_skip()), 
  trim_ws = TRUE
) %>% 
  rename(
    sample.name = names
  )

units <- tibble(
  sample.name = sample_sheet$sample.name,
  read1 = list.files(path = "data/fastq", pattern = "*_R1", full.names = TRUE),
  read2 = list.files(path = "data/fastq", pattern = "*_R2", full.names = TRUE)
)

write_delim(sample_sheet,
          file = "config/sample_sheet.tsv",
          delim = "\t")

write_delim(units,
            file = "config/units.tsv",
            delim = "\t")
