library(dplyr)
library(readr)
library(purrr)
library(clusterProfiler)
# Metadata
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
dotplot(en)


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


count_no_accum <- function(info, genome, eggnog){
  genome_counts <- read_tsv(paste0("sandbox/tximport/", genome, "_counts_raw.tsv"))
  info_subset <- info %>%
    filter(diagnosis == "CD")
  genome_subset <- genome_counts %>%
    select(protein, info_subset$library_name)
  genome_subset_sums <- rowSums(genome_subset[ , -1])
  print(table(genome_subset_sums == 0))
  genome_subset_missing <- genome_subset[genome_subset_sums == 0 , ]
  genome_subset_missing_eggnog <- eggnog %>%
    filter(query_name %in% genome_subset_missing$protein) %>%
    filter(genome == genome)
  return(genome_subset_missing_eggnog)
}

cag1024 <- count_no_accum(info = info, genome= no_accum_cd[1], eggnog)
View(cag1024)
tmp <- lapply(X = no_accum_cd, FUN = function(x) count_no_accum(info = info, genome = x, eggnog = eggnog))
en <- enricher(gene = gsub(",.*", "", tmp$KEGG_ko), 
               TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(en)

mean(c(2089, 330, 314, 122, 113, 95, 144, 76,279, 97))
