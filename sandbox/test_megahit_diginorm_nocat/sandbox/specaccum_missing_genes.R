library(dplyr)
library(readr)
library(purrr)
library(clusterProfiler)

# Pathways etc -----

# pathways <- "http://rest.kegg.jp/link/pathway/ko"
# download.file(url = pathways, destfile = "inputs/kegg_pathways.tsv")
pathways <- read_tsv("inputs/kegg_pathways.tsv", col_names = c("KO", "path"))
pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

# pathway_names <- "http://rest.kegg.jp/list/pathway"
# download.file(url = pathway_names, destfile = "inputs/kegg_pathway_names.tsv")
pathway_names <- read_tsv("inputs/kegg_pathway_names.tsv", col_names = c("path", "name"))

# Metadata -----
info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct()

info_uc <- info %>%
  filter(diagnosis == "UC")

eggnog_files <- list.files("sandbox/test_megahit_diginorm_nocat/eggnog", 
                           ".annotations$", full.names = T)

eggnog <- eggnog_files %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, .id = "genome", skip = 3) %>%
  mutate(genome = gsub(".emapper.annotations", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/eggnog/", "", genome)) %>%
  rename(query_name = "#query_name")

# UC c bolt --------------

cbolt <- read_tsv("sandbox/tximport/GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna_counts_raw.tsv")
cbolt_uc <- cbolt %>%
  select(protein, info_uc$library_name)
cbolt_uc_sums <- rowSums(cbolt_uc[ , -1])
table(cbolt_uc_sums == 0)
cbolt_uc_missing <- cbolt_uc[cbolt_uc_sums == 0 , ]
View(cbolt_uc_missing)
cbolt_uc_missing_eggnog <- eggnog %>%
  filter(query_name %in% cbolt_uc_missing$protein) %>%
  filter(genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna")
table(is.na(cbolt_uc_missing_eggnog$KEGG_ko))
en <- enricher(gene = gsub(",.*", "", cbolt_uc_missing_eggnog$KEGG_ko), 
               TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(en)


# no accum cd -------------------------------------------------------------

no_accum_cd <- c('ERS235530_10.fna', 'ERS396297_11.fna', 'ERS608576_22.fna', 
                 'LiJ_2014__O2.UC28-1__bin.61.fa', 'LiSS_2016__FAT_DON_8-22-0-0__bin.28.fa', 
                 'NielsenHB_2014__MH0094__bin.44.fa', 'SRR6028281_bin.3.fa', 
                 'SRS103987_37.fna', 'SRS104400_110.fna',  'SRS476209_42.fna')
genome <- read_tsv(paste0("sandbox/tximport/", 'ERS235530_10.fna', "_counts_raw.tsv"))
info_subset <- info %>%
  filter(diagnosis == "CD")
genome_subset <- genome %>%
  select(protein, info_subset$library_name)
genome_subset_sums <- rowSums(genome_subset[ , -1])
min(genome_subset_sums)
print(table(genome_subset_sums == 0))
genome_subset_missing <- genome_subset[genome_subset_sums == 0 , ]

## note these functions are written for CD only
count_no_accum <- function(info, genome){
  genome_counts <- read_tsv(paste0("sandbox/tximport/", genome, "_counts_raw.tsv"))
  info_subset <- info %>%
    filter(diagnosis == "CD")
  genome_subset <- genome_counts %>%
    select(protein, info_subset$library_name)
  genome_subset_sums <- rowSums(genome_subset[ , -1])
  print(table(genome_subset_sums == 0))
  genome_subset_missing <- genome_subset[genome_subset_sums == 0 , ]
  return(genome_subset_missing)
}


return_eggnog_no_accum <- function(info, genome, eggnog){
  genome_subset_missing <- count_no_accum(info, genome) 
  genome_subset_missing_eggnog <- eggnog %>%
    filter(query_name %in% genome_subset_missing$protein) %>%
    filter(genome == genome)
  return(genome_subset_missing_eggnog)
}


# test on one
cag1024 <- return_eggnog_no_accum(info = info, genome= no_accum_cd[1], eggnog)
#View(cag1024)

# apply to all
# tmp <- lapply(X = no_accum_cd, FUN = function(x) count_no_accum(info = info, genome = x))
tmp <- lapply(X = no_accum_cd, FUN = function(x) return_eggnog_no_accum(info = info, genome = x, eggnog = eggnog))


# concatenate kegg
tmp <- do.call(rbind, tmp)
en <- enricher(gene = gsub(",.*", "", tmp$KEGG_ko), 
               TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(en)





# make a table with num genes in pangenome, and num missing in CD  -------------

tximport_files <- list.files("sandbox/tximport", ".tsv$", full.names = T)
info_subset_cd <- info %>%
  filter(diagnosis == "CD")
info_subset_uc <- info %>%
  filter(diagnosis == "UC")
info_subset_nonibd <- info %>%
  filter(diagnosis == "nonIBD")
all_missing <- data.frame()

for(i in 1:length(tximport_files)){
  # read in pangenome counts
  tximport <- tximport_files[i]
  tximport <- read_tsv(tximport)
  num_genes_in_pangenome <- nrow(tximport)
  # subset to disease subtype
  cd_subset <- tximport %>%
    select(protein, info_subset_cd$library_name)
  uc_subset <- tximport %>%
    select(protein, info_subset_uc$library_name)
  nonibd_subset <- tximport %>%
    select(protein, info_subset_nonibd$library_name)
  # sum across rows 
  cd_subset_sums <- rowSums(cd_subset[ , -1])
  uc_subset_sums <- rowSums(uc_subset[ , -1])
  nonibd_subset_sums <- rowSums(nonibd_subset[ , -1])
  # subset to zero
  cd_subset_missing <- cd_subset[cd_subset_sums == 0 , ]
  uc_subset_missing <- uc_subset[uc_subset_sums == 0 , ]
  nonibd_subset_missing <- nonibd_subset[nonibd_subset_sums == 0 , ]
  # count number of genes that are unobserved in each subtype
  cd_n_missing <- nrow(cd_subset_missing)
  uc_n_missing <- nrow(uc_subset_missing)
  nonibd_n_missing <- nrow(nonibd_subset_missing)
  # make a dataframe recording this information
  df <- data.frame("pangenome" = tximport_files[i],
                   "pangenome_genes" = num_genes_in_pangenome,
                   "CD" = cd_n_missing,
                   "UC" = uc_n_missing,
                   "nonIBD" = nonibd_n_missing)
  all_missing <- rbind(all_missing, df)
}

all_missing <- all_missing %>%
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

tmp <-left_join(all_missing, gtdb, by = c("pangenome" = "accession")) %>%
  select(labels, pangenome_genes, CD, UC, nonIBD)

mean(tmp$UC)
mean(tmp$CD)
mean(tmp$nonIBD)

mean(tmp$UC/tmp$pangenome_genes) * 100
mean(tmp$CD/tmp$pangenome_genes)* 100
mean(tmp$nonIBD/tmp$pangenome_genes)* 100
