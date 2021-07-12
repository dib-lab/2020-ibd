library(dplyr)
library(readr)
library(purrr)

fstat <- Sys.glob("ntcard_kmer_count/*fstat") %>%
  set_names() %>%
  map_dfr(read_tsv, n_max = 1, col_names = c("k", "F", "num_kmers"), 
          .id = "library_name") %>%
  mutate(library_name = gsub(".abundtrim.fstat", "", library_name)) %>%
  mutate(library_name = gsub("ntcard_kmer_count\\/", "", library_name)) %>%
  select(library_name, num_kmers)

write_tsv(fstat, 'ntcard_kmer_count/all_kmer_count.tsv')
