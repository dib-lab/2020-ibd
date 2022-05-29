library(dplyr)
library(readr)
library(tidyr)
library(purrr)

setwd("~/github/2020-ibd/")

# read in summarized lineages from gather results
gather_lineage <- list.files("sandbox/test_gather_lineage_summarize", full.names = T, pattern = ".tsv$") %>%
  set_names() %>%
  map_dfr(read_delim, delim = ";", skip=3, col_names = c("tmp", "phylum", "class", "order", "family", "genus", "species"), .id = "source")  %>%
  separate(col = tmp, into = c("level", "fraction", "superkingdom"), sep = " ") %>%
  mutate(source = gsub("sandbox/test_gather_lineage_summarize/", "", source)) %>%
  mutate(source = gsub("_lineage_summarized.tsv", "", source)) %>%
  mutate(source = gsub("_vita_vars_genbank", "", source )) %>%
  separate(col = source, into = c("study", "seed"), sep = "_") %>%
  mutate(lineage = paste(superkingdom, phylum, class, order, family, genus, species, sep = ";"))

# filter to only contain species
species <- gather_lineage %>%
  filter(level == "species") 

# isolate species that were identified in all models (n = 155)
shared_gather_lineages <- species %>%
  group_by(superkingdom, phylum, class, order, family, genus, species) %>%
  tally() %>%
  filter(n == 36)
# replace the "unassigneds" with NAs to match formating with lineages DF
shared_gather_lineages[ shared_gather_lineages == "unassigned" ] <- NA

# read in lineage information and filter to the species that were identified
# in all models (n = 155) 
all_lineages <- read_csv("sandbox/test_gather_lineage_summarize/all_genbank_lineages.csv") %>%
  right_join(shared_gather_lineages)
# check for NAs in accession
table(is.na(all_lineages$accession))

# read in gather results and subset to accessions in the all_lineages df
gather_results <- Sys.glob("outputs/gather/*genbank*seed*csv") %>%
  set_names() %>%
  map_dfr(read_csv, .id = "source")  %>%
  separate(col = name, into = c("accession"), remove = F, sep = " ", extra = "drop") %>%
  mutate(accession = gsub("\\..*", "", accession)) %>%
  filter(accession %in% all_lineages$accession)

accessions <- gather_results %>%
  select(accession) %>%
  distinct()

accessions_to_lineage <- left_join(accessions, all_lineages)


tmp <- gather_results %>%
  group_by(accession) %>%
  tally() %>%
  filter(n == 36) %>%
  left_join(all_lineages)
tmp %>%
  select(species) %>%
  distinct()


# what fraction of hashes from each model is classifiable with ncbi --------
gather_frac <- gather_results %>%
  group_by(source) %>%
  summarize(total_original_query = sum(f_unique_to_query))