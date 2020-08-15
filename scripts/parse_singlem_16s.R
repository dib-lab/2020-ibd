library(dplyr) 
library(purrr)
library(readr)
library(tidyr)


#files16 <- list.files(path = "outputs/sgc_genome_queries_singlem", 
#                    pattern = "_otu_16s.csv$", 
#                    recursive = T, full.names = T)
files16 <- unlist(snakemake@input[['s16']])
s16 <- files16 %>%
  set_names() %>% 
  map_dfr(read_tsv, col_types = "cccddcccl", .id = "source") %>%
  mutate(source = gsub("outputs\\/sgc_genome_queries_singlem\\/", "", source)) %>%
  mutate(source = gsub("_otu_16s\\.csv", "", source)) %>%
  separate(source, into = c("sample", "gather_genome"), remove = T, sep = "/")

write_tsv(s16, snakemake@output[['res']])
