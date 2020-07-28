library(dplyr)
library(readr)
library(purrr)

res <- read_csv(unlist(snakemake@input)) %>%
  distinct() %>%
  select(filename, record_name) %>%
  distinct() 

write_csv(snakemake@output[['names']], col_names = F)