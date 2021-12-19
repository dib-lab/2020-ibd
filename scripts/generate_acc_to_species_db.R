library(dplyr)
library(readr)

lineages <- read_csv(snakemake@input[['lineages']], 
                     col_names = c("accession", "domain", "phylum", "class", 
                                   "order", "family", "genus", "species")) %>%
# lineages <- read_csv("outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.lineages.csv", 
#                      col_names = c("accession", "domain", "phylum", "class", 
#                                    "order", "family", "genus", "species")) %>%
  mutate(accession = gsub("_genomic.fna.gz", "", accession),
         species = gsub(" ", "_", species)) %>%
  select(accession, species)

write_csv(lineages, snakemake@output[['csv']])

# cat(paste0("/home/ntpierce/2021-rank-compare/output.rank-compare/sourmash-nodegraph/species/gtdb-rs202.", lineages$species, ".protein-k10.nodegraph"), sep = " ")
# 
# cp: cannot stat '/home/ntpierce/2021-rank-compare/output.rank-compare/sourmash-nodegraph/species/gtdb-rs202.s__CAG-170_sp900545925.protein-k10.nodegraph': No such file or directory
# cp: cannot stat '/home/ntpierce/2021-rank-compare/output.rank-compare/sourmash-nodegraph/species/gtdb-rs202.s__NK3B98_sp900758315.protein-k10.nodegraph': No such file or directory
# cp: cannot stat '/home/ntpierce/2021-rank-compare/output.rank-compare/sourmash-nodegraph/species/gtdb-rs202.s__Angelakisella_massiliensis.protein-k10.nodegraph': No such file or directory
