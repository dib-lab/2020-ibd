library(dplyr)
library(readr)
library(tidyr)

# read in gather csv for shared assemblies and db lineages.
# join together, format accession to filename, and output. 
# required input for charcoal decontamination

gather_shared_assemblies <- read_csv(snakemake@input[["gather"]]) %>%
  select(name) %>%
  separate(name, into = "ident", sep = " ")
db_lineages <- read_csv(snakemake@input[["db_lineages"]])
gather_shared_assemblies_lineages <- left_join(gather_shared_assemblies, db_lineages, by = "ident") %>%
  mutate(ident = paste0(ident, "_genomic.fna.gz"))
write_csv(gather_shared_assemblies_lineages, snakemake@output[["lineages"]], 
          col_names = F)
