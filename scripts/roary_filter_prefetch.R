library(readr)
library(dplyr)

prefetch <- read_csv(snakemake@input[['prefetch']]) %>%
  filter(jaccard >= .1)

write_csv(prefetch, snakemake@output[['filt']])
