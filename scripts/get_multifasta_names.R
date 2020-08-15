library(dplyr)
library(readr)
library(purrr)

res <- read_csv(unlist(snakemake@input)) %>%
  distinct() %>%
  select(filename, record_name) %>%
  distinct() %>%
  mutate(record_name = gsub(" .*", "", record_name)) 

write.table(res$record_name, snakemake@output[['names']], col.names = F, quote = F, row.names = F)
