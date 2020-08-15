library(dplyr)
library(readr)
library(tidyr)

singlem <- read_tsv(snakemake@input[['res']], col_types = "ccccddcccl")

singlem <- singlem %>%
  select(sample, sequence, num_hits) %>%
  group_by(sample, sequence) %>%
  summarise(num_hits = sum(num_hits))

singlem <- pivot_wider(singlem, id_cols = sample, names_from = sequence, values_from = num_hits)
write_tsv(singlem, snakemake@output[['counts']])
