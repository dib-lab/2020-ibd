library(readr)
library(dplyr)
library(ggplot2)
library(vegan)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct()

generate_counts <- function(path){
  # read in an format counts
  counts <- read_tsv(path)
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$protein
  counts <- counts[ , -1]
}

generate_presence_absence <- function(path, threshold = 0){
  counts <- generate_counts(path)
  # replace any count > 0 with 1
  counts_pa <- ifelse(counts > threshold, 1, 0)
}

test_gene_abund <- function(counts_pa){
  n_proteins_observed_per_sample <- data.frame(library_name = colnames(counts_pa), 
                                               n_prots = colSums(counts_pa))
  
  n_proteins_observed_per_sample <- full_join(n_proteins_observed_per_sample, 
                                              info)
  tukey <- TukeyHSD(aov(n_prots~diagnosis, data = n_proteins_observed_per_sample))
  return(tukey)
}

mean_gene_abund <- function(counts_pa){
  n_proteins_observed_per_sample <- data.frame(library_name = colnames(counts_pa), 
                                               n_prots = colSums(counts_pa))
  
  n_proteins_observed_per_sample <- full_join(n_proteins_observed_per_sample, info)
  mean_proteins_observed_per_group <- n_proteins_observed_per_sample %>%
    group_by(diagnosis) %>%
    summarize(mean = mean(n_prots)) %>%
    ungroup() %>%
    column_to_rownames("diagnosis") %>%
    t() %>%
    as.data.frame()
  return(mean_proteins_observed_per_group)
}


tximport_files <- list.files("sandbox/tximport", ".tsv$", full.names = T)


counts_pa <- generate_presence_absence(path = tximport_files[1], threshold = 0)
tukey_counts <- test_gene_abund(counts_pa = counts_pa)
mean_counts <- mean_gene_abund(counts_pa = counts_pa)

all_mean_counts<- data.frame()
for(i in 1:length(tximport_files)){
  counts_pa <- generate_presence_absence(path = tximport_files[i], threshold = 0)
  mean_counts <- mean_gene_abund(counts_pa = counts_pa)
  mean_counts$pangenome <- tximport_files[i]
  all_mean_counts <- rbind(all_mean_counts, mean_counts)
}

###  format table
all_mean_counts <- all_mean_counts %>%
  mutate(pangenome = gsub("_counts_raw\\.tsv", "", pangenome)) %>%
  mutate(pangenome = gsub("sandbox\\/tximport\\/", "", pangenome)) %>%
  distinct() %>%
  mutate(pangenome = gsub("\\.fna", "", pangenome)) %>%
  mutate(pangenome = gsub("\\.fa", "", pangenome))


# join with taxonomy info 
gtdb <- read_tsv("sandbox/gtdbtk_41_genomes/gtdbtk.bac120.summary.tsv") %>%
  mutate(accession = user_genome) %>%
  select(accession, classification) %>%
  mutate(classification = gsub("[dkpcofgs]__", "", classification)) %>%
  separate(data = ., col = classification, sep = ";",
           into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
  mutate(strain = NA) %>%
  mutate(source = "gtdbtk")%>%
  #transmute(labels = ifelse(species == "", genus, species), across(everything()))
  mutate(labels = ifelse(species == "", genus, species)) %>%
  select(accession, labels)

tmp <-left_join(all_mean_counts, gtdb, by = c("pangenome" = "accession")) %>%
  select(labels, CD, UC, nonIBD) %>%
  mutate(CD = round(CD)) %>%
  mutate(UC = round(UC)) %>%
  mutate(nonIBD = round(nonIBD)) %>%
  arrange(desc(nonIBD))

write_tsv(tmp, "tmp_mean_gene_counts.tsv")
  