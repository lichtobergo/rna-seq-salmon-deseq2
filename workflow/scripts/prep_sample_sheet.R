library(tidyverse)
library(readxl)

sample_sheet_p1424 <- read_excel(
  "~/mnt/projects_share/abSynT x cTOVA RNASeq_LD p1424.xlsx",
  skip = 2
) %>%
  dplyr::select(1:8) %>%
  rename(
    sample.type = "...2",
    animal.number = "...7",
    sequencing.status = "...8"
  )

colnames(sample_sheet_p1424) <- str_replace_all(colnames(sample_sheet_p1424), " ", ".")

sample_sheet_p1424 <- sample_sheet_p1424 %>%
  mutate_all(
    ~ str_replace(., " ", "_")
  ) %>%
  mutate(
    tissue = case_when(
      str_detect(sample.name, pattern = "[1-3]") ~ "ctx",
      str_detect(sample.name, pattern = "[4-6]") ~ "blood"
    ),
    sample.group = str_split_i(sample.group, "_", i = 1),
    group = paste(sample.group, tissue, sep = ".")
  ) %>%
  dplyr::filter(
    is.na(sequencing.status)
  ) %>% 
  select(
    - sequencing.status
  ) %>% 
  arrange(sample.name)


write_delim(sample_sheet_p1424,
            file = "config/sample_sheet.tsv",
            delim = "\t")

write_delim(units,
            file = "config/units.tsv",
            delim = "\t")


units <- tibble(
  read1 = list.files(
    path = "data/fastq",
    pattern = "*_R1",
    recursive = TRUE,
    full.names = TRUE
  ),
  read2 = list.files(
    path = "data/fastq",
    pattern = "*_R2",
    recursive = TRUE,
    full.names = TRUE
  )
)

units <- units %>%
  mutate(
    sample.group = str_extract(
      string = read1,
      pattern = "[0-9](M|R)"
    ),
    sample.name = str_extract(
      string = read1,
      pattern = "[0-9](M|R)"
    )
  ) %>%
  select(
    sample.name,
    sample.group,
    read1,
    read2
  )

write_delim(sample_sheet_p1424,
            file = "config/sample_sheet.tsv",
            delim = "\t")

write_delim(units,
            file = "config/units.tsv",
            delim = "\t")
