library(dplyr)
library(readr)
library(tidyr)

singlem <- read_tsv(snakemake@input[['res']], col_types = "ccccddcccl")

singlem <- singlem %>%
  select(source, sequence, num_hits) %>%
  group_by(source, sequence) %>%
  summarise(num_hits = sum(num_hits))

singlem <- pivot_wider(singlem, id_cols = source, names_from = sequence, values_from = num_hits)
singlem <- singlem %>%
  rename(sample = source)

write_tsv(singlem, snakemake@output[['counts']])
