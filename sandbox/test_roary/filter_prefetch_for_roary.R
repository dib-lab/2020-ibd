library(readr)
library(dplyr)

prefetch <- read_csv("GCF_008121495.1_genomic_prefetch.csv") %>%
  filter(jaccard >= .1) %>%
  rename(name = match_name)

#write_csv(prefetch, "GCF_008121495.1_genomic_prefetch_filtered.csv")
write_csv(prefetch, "outputs/genbank/roary_prefetch.x.genbank.gather.csv")
