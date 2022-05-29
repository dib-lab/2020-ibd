library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(rjson)
library(clusterProfiler)
# read in library info
info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, diagnosis) %>%
  mutate(library_name = gsub("-", "\\.", library_name)) %>%
  distinct()


# HASH RESULTS ------------------------------------------------------------

# read in corncob hash differential abundance results
corncob_hash <- read_tsv("sandbox/hash_corncob/sig_ccs_hashes.tsv") %>%
  mutate(hashval = as.character(aa_seq)) %>%
  select(-aa_seq)

# read in sgc multifasta query annotations
sgc_annot <- read_csv("outputs/gather_matches_loso_multifasta/all-multifasta-query-results.csv") %>%
  filter(hashval != "hashval") %>%
  separate(record_name, sep = " ", into = c("sequence_name", "prokka_annotation"), extra = "merge") %>%
  mutate(hashval = as.character(hashval)) %>%
  select(-catlas_base) %>%
  distinct()

hash_results <- left_join(corncob_hash, sgc_annot, by = "hashval")

# read in eggnog gather genome prokka annotations
eggnog_hash <- read_tsv("outputs/gather_matches_loso_multifasta/all-multifasta-query-results.emapper.annotations", 
                   comment = "#", 
                   col_names = c('query_name', 'seed_eggNOG_ortholog',	
                                 'seed_ortholog_evalue', 'seed_ortholog_score',
                                 'best_tax_level', 'Preferred_name', 'GOs', 'EC',
                                 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',	'KEGG_TC',
                                 'CAZy', 'BiGG_Reaction', 'taxonomic_scope', 
                                 'eggNOG_OGs', 'best_eggNOG_OG', 'COG_functional_cat',	
                                 'eggNOG'))
eggnog_hash$KEGG_ko <- gsub(",.*", "", eggnog_hash$KEGG_ko) # select first KO annotation


hash_results <- left_join(hash_results, eggnog_hash, by = c("sequence_name" = "query_name"))

# Read in variable importance
read_varimp <- function(path_optimal_rf){
  study <- gsub("outputs\\/optimal_rf\\/", "", path_optimal_rf)
  study <- gsub("_optimal_rf\\.RDS", "", study)
  optimal_rf <- readRDS(path_optimal_rf)
  varimp <- data.frame(hash = names(optimal_rf$variable.importance), 
                       importance = optimal_rf$variable.importance,
                       study = study)
  rownames(varimp) <- seq(1:nrow(varimp))
  # add a column where varimp is normalized by the total var imp
  # e.g., divide by the sum of all variable importances
  # this will make all variable importance measures sum to 1
  varimp$norm <- varimp$importance / sum(varimp$importance)
  return(varimp)
}

ihmp_varimp <- read_varimp("outputs/optimal_rf/iHMP_optimal_rf.RDS")
prjeb2054_varimp <- read_varimp("outputs/optimal_rf/PRJEB2054_optimal_rf.RDS")
prjna237362_varimp <- read_varimp("outputs/optimal_rf/PRJNA237362_optimal_rf.RDS")
prjna385949_varimp <- read_varimp("outputs/optimal_rf/PRJNA385949_optimal_rf.RDS")
prjna400072_varimp <- read_varimp("outputs/optimal_rf/PRJNA400072_optimal_rf.RDS")
srp057027_varimp <- read_varimp("outputs/optimal_rf/SRP057027_optimal_rf.RDS")

varimp <- rbind(ihmp_varimp, prjeb2054_varimp, prjna237362_varimp, 
                prjna385949_varimp, prjna400072_varimp, srp057027_varimp)

varimp_cum <- varimp %>%
  group_by(hash) %>%
  summarise(total_imp = sum(norm)) %>%
  left_join(varimp) %>%
  arrange(desc(total_imp)) %>%
  select(hash, total_imp) %>%
  distinct() %>%
  mutate(total_imp = total_imp/6)
# varimp_cum$hash <- as.numeric(as.character(varimp_cum$hash))

hash_results <- left_join(hash_results, varimp_cum, by = c("hashval" = "hash"))

unassembled_hashes <- read_csv("sandbox/test_megahit_diginorm_nocat/megahit_gather/at_least_5_vita_vars_vs_megahit_assemblies_tbp0.un.csv") %>%
  mutate(hashval = as.character(minhash)) %>%
  select(hashval) %>%
  mutate(assembled = "unassembled") 

hash_results <- left_join(hash_results, unassembled_hashes, by = "hashval")
hash_results$assembled <- ifelse(is.na(hash_results$assembled), "assembled", hash_results$assembled)

# read in assembled/unassembled information
hash_results <- hash_results %>%
  mutate(gather_genome = gsub(".fna.*", "", sequence_name)) %>%
  mutate(gather_genome = gsub(".fa.*", "", gather_genome))%>%
  mutate(gather_genome = gsub("_232_.*", "", gather_genome)) %>%
  mutate(direction = ifelse(estimate > 0, "up", "down")) %>%
  mutate(set = gsub("mu\\.diagnosis", "", paste0(mu, "_", direction))) 


# PANGENOME RESULTS -------------------------------------------------------

eggnog_files <- list.files("sandbox/test_megahit_diginorm_nocat/eggnog", 
                           ".annotations$", full.names = T)



eggnog_pan <- eggnog_files %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, .id = "genome", comment = "#", 
          col_names = c('query_name', 'seed_eggNOG_ortholog',	
                        'seed_ortholog_evalue', 'seed_ortholog_score',
                        'best_tax_level', 'Preferred_name', 'GOs', 'EC',
                        'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                        'KEGG_Reaction', 'KEGG_rclass', 'BRITE',	'KEGG_TC',
                        'CAZy', 'BiGG_Reaction', 'taxonomic_scope', 
                        'eggNOG_OGs', 'best_eggNOG_OG', 'COG_functional_cat',	
                        'eggNOG')) %>%
  mutate(genome = gsub(".emapper.annotations", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/eggnog/", "", genome))

eggnog_pan$KEGG_ko <- gsub(",.*", "", eggnog_pan$KEGG_ko) # select first KO annotation


corncob_files <- list.files("sandbox/test_megahit_diginorm_nocat/corncob",
                            ".sig_ccs.tsv$", full.names = T) 
corncob_pan <- corncob_files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genome") %>%
  mutate(genome = gsub("_sig_ccs.tsv", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/corncob/", "", genome)) 

pan_results <- left_join(corncob_pan, eggnog_pan, by = c("genome", "aa_seq" = "query_name")) %>%
  mutate(genome = gsub(".fa", "", genome)) %>%
  mutate(genome = gsub(".fna", "", genome))

species <- read_tsv("sandbox/test_megahit_diginorm_nocat/gtdbtk/41genomes_out/gtdbtk.bac120.summary.tsv") %>%
  select(user_genome, classification) %>%
  mutate(classification = gsub("[dkpcofgs]__", "", classification)) %>%
  separate(data = ., col = classification, sep = ";",
           into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))

pan_results <- left_join(pan_results, species, by = c("genome" = "user_genome"))


# HASH:PANGENOME GENE MAP -------------------------------------------------

sigs <- list.files("sandbox/test_megahit_diginorm_nocat/sandbox/comp_hashes_to_assemblies", 
                   pattern = "matches.sig*", full.names = T)
hash_to_pangenome <- data.frame()
for(j in 1:length(sigs)){
  sig_json <- fromJSON(file = sigs[j])
  for(i in 1:length(sig_json)){
    sig <- sig_json[i]
    df <- as.data.frame(sig) %>%
      select(filename, name, signatures.mins) %>%
      mutate(hash = as.character(signatures.mins)) %>%
      select(-signatures.mins) %>%
      separate(name, into = c("aa_seq", "prokka"), sep = " ", extra = "merge") %>%
      mutate(source = gsub("sandbox\\/test_megahit_diginorm_nocat\\/sandbox/comp_hashes_to_assemblies\\/", "", sigs[j])) %>%
      mutate(source = gsub("_matches\\.sig", "", source)) %>%
      mutate(genome = gsub(".fna.cdhit.ffn", "", filename))
    hash_to_pangenome <- rbind(hash_to_pangenome, df)
  }
}


# only hashes that were from hash-level analysis
hash_to_pangenome <- hash_to_pangenome %>%
  filter(hash %in% hash_results$hashval)

# join to pangenome results for annotations/corncob estimates
pan_results_only_hashes <- left_join(hash_to_pangenome, pan_results, by = c("genome", "aa_seq")) %>%
  arrange(source, hash)

# combine with hash-level results to compare

pan_results_only_hashes %>%
  select(hash, aa_seq, prokka, estimate, seed_eggNOG_ortholog, KEGG_ko)

# hash_results_only_hashes <- hash_results %>%
#   filter(gather_genome == "GCF_900036035.1_RGNV35913_genomic") %>%
#   select(hashval, sequence_name, prokka_annotation, estimate, seed_eggNOG_ortholog, KEGG_ko, mu, bonferroni) %>%
#   filter(hashval %in% pan_results_only_hashes$hash) %>%
#   arrange(hashval) %>%
#   select(hashval, estimate) %>%
#   distinct()


tmp <- pan_results_only_hashes %>%
  mutate(sig_in_pangenome = ifelse(is.na(estimate), "no", "yes")) %>%
  select(hash, source, sig_in_pangenome) %>%
  distinct() %>%
  group_by(source, sig_in_pangenome) %>%
  tally()

pan_results_only_hashes %>%
  select(hash, source) %>%
  distinct() %>%
  group_by(source) %>%
  tally()

# what's going on with the one's that are getting dropped as not significant?
corncob_files_all <- list.files("sandbox/test_megahit_diginorm_nocat/corncob",
                                ".all_ccs.tsv$", full.names = T) 
corncob_pan_all <- corncob_files_all %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genome") %>%
  mutate(genome = gsub("_all_ccs.tsv", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/corncob/", "", genome))

corncob_pan_all_cd <- corncob_pan_all %>%
  filter(mu %in% c("mu.diagnosisCD")) %>%
  mutate(bonferroni = p.adjust(p_value, method = "bonferroni"))
corncob_pan_all_uc <- corncob_pan_all %>%
  filter(mu %in% c("mu.diagnosisUC")) %>%
  mutate(bonferroni = p.adjust(p_value, method = "bonferroni"))
corncob_pan_all <- rbind(corncob_pan_all_cd, corncob_pan_all_uc)
corncob_pan_all$genome <- gsub(".fa", "", corncob_pan_all$genome)
corncob_pan_all$genome <- gsub(".fna", "", corncob_pan_all$genome)


tmp2 <- pan_results_only_hashes %>%
  select(filename, aa_seq, prokka, hash, source, genome) %>%
  left_join(corncob_pan_all, by = c('aa_seq', 'genome'))
View(tmp2)
tmp2 %>% 
  mutate(pangenome_direction = ifelse(estimate < 0, "down", "up")) %>%
  select(source, hash, pangenome_direction) %>%
  distinct() %>%
  group_by(source, pangenome_direction) %>%
  tally()
  
View(tmp2 %>% 
  mutate(pangenome_direction = ifelse(estimate < 0, "down", "up")) %>%
  select(source, hash, pangenome_direction) %>%
  distinct())

