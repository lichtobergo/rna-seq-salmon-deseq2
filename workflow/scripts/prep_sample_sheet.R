library(tidyverse)
library(readxl)

sample_sheet_p1425 <- read_excel(
  "~/mnt/projects_share/abSynT x cTOVA RNASeq total_LMK p1425.xlsx", 
  skip = 2
) %>% 
  dplyr::select(1:6) %>% 
  rename(animal.number = "...6")
colnames(sample_sheet_p1425) <- str_replace_all(colnames(sample_sheet_p1425), " ", ".")
sample_sheet_p1425 <- sample_sheet_p1425 %>% 
  mutate_all(
    ~ str_replace(., " ", "_")
    ) %>% 
  mutate(
    tissue = str_split_i(sample.group, "_", i = 2),
    tissue = if_else(
      tissue == "spinal cord",
      "SC",
      tissue
    ),
    sample.group = str_split_i(sample.group, "_", i = 1),
    group = case_when(
      str_detect(sample.group, "CY") ~ "TovaCAR.+AAV",
      str_detect(sample.group, "TY") ~ "Tova.+AAV",
      str_detect(sample.group, "CN") ~ "TovaCAR.noAAV",
      str_detect(sample.group, "TN") ~ "Tova.noAAV"
    ),
    group = paste(tissue, group, sep = ".")
  )

units <- tibble(
  sample.name = sample_sheet_p1425$sample.name,
  read1 = list.files(path = "data/fastq", pattern = "*_R1", recursive = T, full.names = TRUE),
  read2 = list.files(path = "data/fastq", pattern = "*_R2", recursive = T, full.names = TRUE)
)

write_delim(sample_sheet_p1425,
          file = "config/sample_sheet.tsv",
          delim = "\t")

write_delim(units,
            file = "config/units.tsv",
            delim = "\t")
