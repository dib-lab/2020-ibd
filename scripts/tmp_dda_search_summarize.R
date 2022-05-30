library(dplyr)
library(readr)
library(purrr)
library(tidyr)

# TO DO: ADD IN INFORMATION ABOUT THE NUMBER OF HASHES IN SIG, OR SOMETHING ABOUT THE SIZE

setwd("~/github/2020-ibd")

search_files <- Sys.glob("outputs/sgc_pangenome_catlases_corncob_sequences/*tsv")

# find files that only have one line (the header)
file_lengths <- character()
for(file in search_files){
  file_length <- length(readLines(file)) 
  file_lengths <- c(file_lengths, file_length)
}
# remove files that only have one line
file_lengths <- file_lengths > 1
search_files <- search_files[file_lengths]


search_results <- search_files %>%
  set_names() %>%
  map_dfr(read_csv, col_types = "dccccccc", .id = "accession")

search_results <- search_results %>%
  mutate(accession = gsub("_contigs_search_gtdb_genomic.tsv", "", accession)) %>%
  mutate(accession = gsub("outputs/sgc_pangenome_catlases_corncob_sequences/", "", accession)) %>%
  separate(accession, into = c("prefix", "accession_no_prefix", "disease", "abundance"), sep = "_", remove = F) %>%
  mutate(accession = gsub("_CD_increased", "", accession)) %>%
  filter(grepl("GCF", name))  %>%# only keep genomes in refseq
  filter(similarity >= 0.5) %>%
  filter(!grepl("ERS", name)) # filter out a suppressed refseq accession

tmp <- search_results %>%
  mutate(name_no_accession = gsub("GCF_[0-9\\.]{11} ", "", name)) %>%
  group_by(name_no_accession) %>%
  tally()

lineages <- read_csv("outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.lineages.csv",
                     col_names = c("accession", "domain", "phylum", "class",
                                   "order", "family", "genus", "species")) %>%
  mutate(accession = gsub("_genomic.fna.gz", "", accession))

search_results <- search_results %>%
  left_join(lineages) 

search_results <- search_results %>%
  select(accession, domain, phylum, class, order, family, genus, species, 
         disease, abundance, similarity, name)

tmp <- search_results %>%
  filter(species %in% c("s__Ruminococcus_B gnavus", "s__Enterocloster sp005845215",
                        "s__Enterocloster clostridioformis_A",
                        "s__Enterocloster clostridioformis",
                        "s__Enterocloster bolteae"))