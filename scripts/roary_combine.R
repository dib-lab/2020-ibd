library(readr)
library(dplyr)
library(purrr)

files <- unlist(snakemake@input[['filt']])
prefetch <- files %>%
  map_dfr(read_csv) %>%
  group_by(match_name) %>%
  slice(1) %>%
  ungroup()
write_tsv(prefetch, snakemake@output[["combined"]]

#accessions <- prefetch %>%
#  separate(match_name, into = "accession", sep = " ") %>%
#  select(accession) %>%
#  distinct()

#write_tsv(accessions, snakemake@output[['combined']], col_names = F)


