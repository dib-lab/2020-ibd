library(dplyr)
library(readr)
library(purrr)

files <- unlist(snakemake@input[["fstat"]])
fstat <- files %>%
  set_names() %>%
  map_dfr(read_tsv, n_max = 1, col_names = c("k", "F", "num_kmers"), 
          .id = "library_name") %>%
  mutate(library_name = gsub("\\.fstat", "", library_name)) %>%
  mutate(library_name = gsub("outputs\\/ntcard\\/", "", library_name)) %>%
  select(library_name, num_kmers)

write_tsv(fstat, snakemake@output[["tsv"]])
