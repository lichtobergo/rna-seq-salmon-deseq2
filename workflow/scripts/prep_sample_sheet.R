library(tidyverse)
library(readxl)

sample_sheet <- read_delim(
  "config/sample_sheet.txt", 
  delim = "\t", escape_double = FALSE, 
  col_types = cols(Comments = col_skip()), 
  trim_ws = TRUE
) %>% 
  dplyr::filter(
    sample.name != "N_04"
  ) %>% 
  relocate(
    sample.name,
    .before = sample.number
  ) %>% 
  arrange(sample.name) %>% 
  mutate(
    sample.name = str_remove(sample.name, "_")
  )
# %>% 
#   dplyr::select(1:6) %>% 
#   rename(animal.number = "...6")
# colnames(sample_sheet_p1425) <- str_replace_all(colnames(sample_sheet_p1425), " ", ".")
# sample_sheet_p1425 <- sample_sheet_p1425 %>% 
#   mutate_all(
#     ~ str_replace(., " ", "_")
#     ) %>% 
#   mutate(
#     tissue = str_split_i(sample.group, "_", i = 2),
#     tissue = if_else(
#       tissue == "spinal cord",
#       "SC",
#       tissue
#     ),
#     sample.group = str_split_i(sample.group, "_", i = 1),
#     group = case_when(
#       str_detect(sample.group, "CY") ~ "TovaCAR.+AAV",
#       str_detect(sample.group, "TY") ~ "Tova.+AAV",
#       str_detect(sample.group, "CN") ~ "TovaCAR.noAAV",
#       str_detect(sample.group, "TN") ~ "Tova.noAAV"
#     )
#   )

units <- tibble(
  sample.name = sample_sheet$sample.name,
  read1 = list.files(path = "data/fastq", pattern = "*_R1", recursive = T, full.names = TRUE),
  read2 = list.files(path = "data/fastq", pattern = "*_R2", recursive = T, full.names = TRUE)
)

write_delim(sample_sheet,
          file = "config/sample_sheet.tsv",
          delim = "\t")

write_delim(units,
            file = "config/units.tsv",
            delim = "\t")
