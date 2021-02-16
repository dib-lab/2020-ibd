library(dplyr)
library(purrr)
library(readr)

files_default <- unlist(snakemake@input[['default']])

defaults <- files_default %>%
  set_names() %>%
  map_dfr(read_tsv, col_types = "cccddcccl", .id = "source") %>%
  mutate(source = gsub("outputs\\/abundtrim_singlem\\/", "", source)) %>%
  mutate(source = gsub("_otu_default\\.csv", "", source)) 

files16_R1 <- unlist(snakemake@input[['s16_R1']])
s16_R1 <- files16_R1 %>%
  set_names() %>%
  map_dfr(read_tsv, col_types = "cccddcccl", .id = "source") %>%
  mutate(source = gsub("outputs\\/abundtrim_singlem\\/", "", source)) %>%
  mutate(source = gsub("_otu_16s_R1\\.csv", "", source))

files16_R2 <- unlist(snakemake@input[['s16_R2']])
s16_R2 <- files16_R2 %>%
  set_names() %>%
  map_dfr(read_tsv, col_types = "cccddcccl", .id = "source") %>%
  mutate(source = gsub("outputs\\/abundtrim_singlem\\/", "", source)) %>%
  mutate(source = gsub("_otu_16s_R2\\.csv", "", source))

all <- rbind(defaults, s16_R1, s16_R2)
write_tsv(all, snakemake@output[['res']])
