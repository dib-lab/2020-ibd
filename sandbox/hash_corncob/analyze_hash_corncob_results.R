library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(clusterProfiler)
# read in library info
info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, diagnosis) %>%
  mutate(library_name = gsub("-", "\\.", library_name)) %>%
  distinct()

# read in corncob hash differential abundance results
corncob <- read_tsv("sandbox/hash_corncob/sig_ccs_hashes.tsv") %>%
  mutate(hashval = as.character(aa_seq)) %>%
  select(-aa_seq)

# read in sgc multifasta query annotations
sgc_annot <- read_csv("outputs/gather_matches_loso_multifasta/all-multifasta-query-results.csv") %>%
  filter(hashval != "hashval") %>%
  separate(record_name, sep = " ", into = c("sequence_name", "prokka_annotation"), extra = "merge") %>%
  mutate(hashval = as.character(hashval)) %>%
  select(-catlas_base) %>%
  distinct()

hash_results <- left_join(corncob, sgc_annot, by = "hashval")

# read in eggnog gather genome prokka annotations
eggnog <- read_tsv("outputs/gather_matches_loso_multifasta/all-multifasta-query-results.emapper.annotations", 
                   comment = "#", 
                   col_names = c('query_name', 'seed_eggNOG_ortholog',	
                                 'seed_ortholog_evalue', 'seed_ortholog_score',
                                 'best_tax_level', 'Preferred_name', 'GOs', 'EC',
                                 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',	'KEGG_TC',
                                 'CAZy', 'BiGG_Reaction', 'taxonomic_scope', 
                                 'eggNOG_OGs', 'best_eggNOG_OG', 'COG_functional_cat',	
                                 'eggNOG'))
eggnog$KEGG_ko <- gsub(",.*", "", eggnog$KEGG_ko) # select first KO annotation


hash_results <- left_join(hash_results, eggnog, by = c("sequence_name" = "query_name"))


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

hash_results_small <- hash_results %>%
  select(hashval, mu, estimate, gather_genome, best_tax_level, prokka_annotation, Preferred_name, 
         GOs, EC, KEGG_ko, CAZy, sequence_name, assembled) %>%
  distinct()

# separate CD/UC up/down --------------------------------------------------

cd_up <- hash_results_small %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate > 0)

cd_down <- hash_results_small %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate < 0)

uc_up <- hash_results_small %>%
  filter(mu == "mu.diagnosisUC") %>%
  filter(estimate > 0)

uc_down <- hash_results_small %>%
  filter(mu == "mu.diagnosisUC") %>%
  filter(estimate < 0)

# plot fractions assembled for up/down
hashval_by_assembled <- hash_results %>%
  select(hashval, assembled, set) %>%
  distinct()
ggplot(hashval_by_assembled, aes(x = set, fill = assembled)) +
  geom_bar()+
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(y = "number of shared hashes") +
  scale_fill_manual(values = c(assembled = "#dddddd", unassembled = "#858585")) +
  scale_x_discrete(labels = c("CD_down" = "CD down", "CD_up" = "CD up",
                              "UC_down" = "UC down", "UC_up" = "UC up"))


# plot per disease num hashes annotated per org in up and down ------------

hash_by_assembled_genome <- hash_results %>%
  select(mu, set, direction, hashval, gather_genome, assembled) %>%
  distinct()

lca <- read_csv("inputs/at_least_5_studies_vita_vars_gather_all_lca.csv") %>%
  select(name_no_gz, GTDB) %>%
  mutate(name = gsub(".fa", "", name_no_gz)) %>%
  mutate(name = gsub(".fna", "", name)) %>%
  select(-name_no_gz)

# cat(paste0("'", lca$name, "' = '", lca$GTDB, "', "))

plt_hash_count <- ggplot(hash_by_assembled_genome, aes(x = forcats::fct_infreq(gather_genome), fill = direction)) +
  geom_bar() +
  facet_wrap(~mu, labeller = as_labeller(c("mu.diagnosisCD" = "CD", "mu.diagnosisUC" = "UC"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(size = 14)) +
  labs(x = "genome", y = "number of hashes") +
  scale_fill_manual(values = c("down" = "#FF7F00", "up" = "#6A3D9A")) +
  scale_x_discrete(labels = c('ERS537353_12' = 'Flavonifractor',  
                              'SRS1719577_6' = 'Faecalibacterium prausnitzii_D',  
                              'SRS294916_20' = 'Romboutsia timonensis', 
                              'SRS1719498_9' = 'Acetatifactor sp900066565',  
                              'ERS396297_11' = 'Lachnospira eligens_B',  
                              'ERS473255_26' = 'Faecalibacterium prausnitzii_G',  
                              'ERS608576_22' = 'Ruminococcus_E bromii_B',  
                              'ERS537328_30' = 'Faecalibacterium prausnitzii_K',  
                              'SRS143598_15' = 'Lachnospira sp000437735',  
                              'SRS476209_42' = 'Ruminococcus_D bicirculans',  
                              'ERS235603_16' = 'Agathobacter rectale',
                              'SRS1719112_8' = 'Oscillibacter sp900066435',
                              'ERS235531_43' = 'Faecalibacterium prausnitzii_F',
                              'SRS1735506_4' = 'Bacteroides ovatus',  
                              'SRS104400_110' = 'Lachnospira sp900316325',
                              'SRR4305229_bin.5' = 'Roseburia inulinivorans',  
                              'SRS103987_37' = 'ER4 sp000765235',  
                              'XieH_2016__YSZC12003_37172__bin.63' = 'Acutalibacter sp000435395',  
                              'LiJ_2014__O2.UC28-1__bin.61' = 'Ruminiclostridium_E siraeum',  
                              'LeChatelierE_2013__MH0074__bin.19' = 'CAG-45 sp900066395', 
                              'GCF_900036035.1_RGNV35913_genomic' = 'Faecalicatena gnavus', 
                              'GCF_000371685.1_Clos_bolt_90B3_V1_genomic' = 'Clostridium_M bolteae', 
                              'ERS608524_37' = 'Gemmiger formicilis',  
                              'LoombaR_2017__SID1050_bax__bin.11' = 'Anaeromassilibacillus sp002159845',
                              'SRS1735645_19' = 'Acetatifactor sp900066365', 
                              'ERS396519_11' = 'Lawsonibacter asaccharolyticus',
                              'SRR6028281_bin.3' = 'Lachnospira eligens_B', 
                              'VatanenT_2016__G80445__bin.9' = 'Faecalibacterium prausnitzii_D', 
                              'ZeeviD_2015__PNP_Main' = 'Blautia_A sp900066165', 
                              'GCF_000508885.1_ASM50888v1_genomic' = 'Flavonifractor sp000508885',
                              'SRS075078_49' = 'TF01-11 sp003529475',  
                              'LiSS_2016__FAT_DON_8-22-0-0__bin.28' = 'CAG-170 sp000432135', 
                              'NielsenHB_2014__MH0094__bin.44' = 'Prevotella copri',  
                              'SRR5127401_bin.3' = 'UBA11774 sp003507655', 
                              'ERS537235_19' = 'Bacteroides_B massiliensis',  
                              'SRR5558047_bin.10' = 'Alistipes putredinis', 
                              'VogtmannE_2016__MMRS43563715ST-27-0-0__bin.70' = 'CAG-81 sp900066785',  
                              'ERS537218_9' = 'Gemmiger sp003476825', 
                              'QinJ_2012__CON-091__bin.20' = 'Faecalibacterium prausnitzii_G', 
                              'GCF_001405615.1_13414_6_47_genomic' = 'Agathobacter faecis',  
                              'ERS235530_10' = 'CAG-1024 sp000432015')) +
  coord_flip()

hash_by_assembled_genome_imp <- hash_results %>%
  select(mu, set, direction, hashval, gather_genome, assembled, total_imp) %>%
  distinct()

plt_var_imp <- ggplot(hash_by_assembled_genome_imp, aes(x = forcats::fct_infreq(gather_genome), y = total_imp, fill = direction)) +
  geom_col() +
  facet_wrap(~mu, 
             labeller = as_labeller(c("mu.diagnosisCD" = "CD", "mu.diagnosisUC" = "UC"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 7),
        #strip.text.x = element_text(size = 14)) +
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x = "genome", y = "cumulative importance", fill = "") +
  scale_fill_manual(values = c("down" = "#FF7F00", "up" = "#6A3D9A")) +
  scale_x_discrete(labels = c('ERS537353_12' = 'Flavonifractor',  
                              'SRS1719577_6' = 'Faecalibacterium prausnitzii_D',  
                              'SRS294916_20' = 'Romboutsia timonensis', 
                              'SRS1719498_9' = 'Acetatifactor sp900066565',  
                              'ERS396297_11' = 'Lachnospira eligens_B',  
                              'ERS473255_26' = 'Faecalibacterium prausnitzii_G',  
                              'ERS608576_22' = 'Ruminococcus_E bromii_B',  
                              'ERS537328_30' = 'Faecalibacterium prausnitzii_K',  
                              'SRS143598_15' = 'Lachnospira sp000437735',  
                              'SRS476209_42' = 'Ruminococcus_D bicirculans',  
                              'ERS235603_16' = 'Agathobacter rectale',
                              'SRS1719112_8' = 'Oscillibacter sp900066435',
                              'ERS235531_43' = 'Faecalibacterium prausnitzii_F',
                              'SRS1735506_4' = 'Bacteroides ovatus',  
                              'SRS104400_110' = 'Lachnospira sp900316325',
                              'SRR4305229_bin.5' = 'Roseburia inulinivorans',  
                              'SRS103987_37' = 'ER4 sp000765235',  
                              'XieH_2016__YSZC12003_37172__bin.63' = 'Acutalibacter sp000435395',  
                              'LiJ_2014__O2.UC28-1__bin.61' = 'Ruminiclostridium_E siraeum',  
                              'LeChatelierE_2013__MH0074__bin.19' = 'CAG-45 sp900066395', 
                              'GCF_900036035.1_RGNV35913_genomic' = 'Faecalicatena gnavus', 
                              'GCF_000371685.1_Clos_bolt_90B3_V1_genomic' = 'Clostridium_M bolteae', 
                              'ERS608524_37' = 'Gemmiger formicilis',  
                              'LoombaR_2017__SID1050_bax__bin.11' = 'Anaeromassilibacillus sp002159845',
                              'SRS1735645_19' = 'Acetatifactor sp900066365', 
                              'ERS396519_11' = 'Lawsonibacter asaccharolyticus',
                              'SRR6028281_bin.3' = 'Lachnospira eligens_B', 
                              'VatanenT_2016__G80445__bin.9' = 'Faecalibacterium prausnitzii_D', 
                              'ZeeviD_2015__PNP_Main' = 'Blautia_A sp900066165', 
                              'GCF_000508885.1_ASM50888v1_genomic' = 'Flavonifractor sp000508885',
                              'SRS075078_49' = 'TF01-11 sp003529475',  
                              'LiSS_2016__FAT_DON_8-22-0-0__bin.28' = 'CAG-170 sp000432135', 
                              'NielsenHB_2014__MH0094__bin.44' = 'Prevotella copri',  
                              'SRR5127401_bin.3' = 'UBA11774 sp003507655', 
                              'ERS537235_19' = 'Bacteroides_B massiliensis',  
                              'SRR5558047_bin.10' = 'Alistipes putredinis', 
                              'VogtmannE_2016__MMRS43563715ST-27-0-0__bin.70' = 'CAG-81 sp900066785',  
                              'ERS537218_9' = 'Gemmiger sp003476825', 
                              'QinJ_2012__CON-091__bin.20' = 'Faecalibacterium prausnitzii_G', 
                              'GCF_001405615.1_13414_6_47_genomic' = 'Agathobacter faecis',  
                              'ERS235530_10' = 'CAG-1024 sp000432015')) +
  coord_flip()

pdf("panel_hum_hashes_var_imp_by_up_down.pdf", height = 10.5, width = 6)
ggarrange(plt_hash_count, plt_var_imp, nrow = 2,
          common.legend = T, legend = "bottom", heights = c(1, .9))
dev.off()




plt_hash_count2 <- ggplot(hash_by_assembled_genome, aes(x = forcats::fct_infreq(gather_genome), 
                                     fill = paste0(direction, "_", assembled))) +
  geom_bar() +
  facet_wrap(~mu, labeller = as_labeller(c("mu.diagnosisCD" = "CD", "mu.diagnosisUC" = "UC"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(size = 14)) +
  labs(x = "genome", y = "number of shared hashes", fill = "") + 
  scale_fill_manual(values = c("down_assembled" = "#FDBF6F", 
                               "down_unassembled" = "#FF7F00",
                               "up_assembled" = "#CAB2D6", 
                               "up_unassembled" = "#6A3D9A"),
                    labels = c("down_assembled" = "down, assembled", 
                               "down_unassembled" = "down, unassembled",
                               "up_assembled" = "up, assembled", 
                               "up_unassembled" = "up, unassembled")) +
  scale_x_discrete(labels = c('ERS537353_12' = 'Flavonifractor',  
                              'SRS1719577_6' = 'Faecalibacterium prausnitzii_D',  
                              'SRS294916_20' = 'Romboutsia timonensis', 
                              'SRS1719498_9' = 'Acetatifactor sp900066565',  
                              'ERS396297_11' = 'Lachnospira eligens_B',  
                              'ERS473255_26' = 'Faecalibacterium prausnitzii_G',  
                              'ERS608576_22' = 'Ruminococcus_E bromii_B',  
                              'ERS537328_30' = 'Faecalibacterium prausnitzii_K',  
                              'SRS143598_15' = 'Lachnospira sp000437735',  
                              'SRS476209_42' = 'Ruminococcus_D bicirculans',  
                              'ERS235603_16' = 'Agathobacter rectale',
                              'SRS1719112_8' = 'Oscillibacter sp900066435',
                              'ERS235531_43' = 'Faecalibacterium prausnitzii_F',
                              'SRS1735506_4' = 'Bacteroides ovatus',  
                              'SRS104400_110' = 'Lachnospira sp900316325',
                              'SRR4305229_bin.5' = 'Roseburia inulinivorans',  
                              'SRS103987_37' = 'ER4 sp000765235',  
                              'XieH_2016__YSZC12003_37172__bin.63' = 'Acutalibacter sp000435395',  
                              'LiJ_2014__O2.UC28-1__bin.61' = 'Ruminiclostridium_E siraeum',  
                              'LeChatelierE_2013__MH0074__bin.19' = 'CAG-45 sp900066395', 
                              'GCF_900036035.1_RGNV35913_genomic' = 'Faecalicatena gnavus', 
                              'GCF_000371685.1_Clos_bolt_90B3_V1_genomic' = 'Clostridium_M bolteae', 
                              'ERS608524_37' = 'Gemmiger formicilis',  
                              'LoombaR_2017__SID1050_bax__bin.11' = 'Anaeromassilibacillus sp002159845',
                              'SRS1735645_19' = 'Acetatifactor sp900066365', 
                              'ERS396519_11' = 'Lawsonibacter asaccharolyticus',
                              'SRR6028281_bin.3' = 'Lachnospira eligens_B', 
                              'VatanenT_2016__G80445__bin.9' = 'Faecalibacterium prausnitzii_D', 
                              'ZeeviD_2015__PNP_Main' = 'Blautia_A sp900066165', 
                              'GCF_000508885.1_ASM50888v1_genomic' = 'Flavonifractor sp000508885',
                              'SRS075078_49' = 'TF01-11 sp003529475',  
                              'LiSS_2016__FAT_DON_8-22-0-0__bin.28' = 'CAG-170 sp000432135', 
                              'NielsenHB_2014__MH0094__bin.44' = 'Prevotella copri',  
                              'SRR5127401_bin.3' = 'UBA11774 sp003507655', 
                              'ERS537235_19' = 'Bacteroides_B massiliensis',  
                              'SRR5558047_bin.10' = 'Alistipes putredinis', 
                              'VogtmannE_2016__MMRS43563715ST-27-0-0__bin.70' = 'CAG-81 sp900066785',  
                              'ERS537218_9' = 'Gemmiger sp003476825', 
                              'QinJ_2012__CON-091__bin.20' = 'Faecalibacterium prausnitzii_G', 
                              'GCF_001405615.1_13414_6_47_genomic' = 'Agathobacter faecis',  
                              'ERS235530_10' = 'CAG-1024 sp000432015')) +
  coord_flip()

plt_var_imp2 <- ggplot(hash_by_assembled_genome_imp, aes(x = forcats::fct_infreq(gather_genome), y = total_imp, 
                                     fill = paste0(direction, "_", assembled))) +
  geom_col() +
  facet_wrap(~mu, labeller = as_labeller(c("mu.diagnosisCD" = "CD", "mu.diagnosisUC" = "UC"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        #strip.text.x = element_text(size = 14)) +
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x = "genome", y = "cumulative importance", fill = "") + 
  scale_fill_manual(values = c("down_assembled" = "#FDBF6F", 
                               "down_unassembled" = "#FF7F00",
                               "up_assembled" = "#CAB2D6", 
                               "up_unassembled" = "#6A3D9A"),
                    labels = c("down_assembled" = "down, assembled", 
                               "down_unassembled" = "down, unassembled",
                               "up_assembled" = "up, assembled", 
                               "up_unassembled" = "up, unassembled")) +
  scale_x_discrete(labels = c('ERS537353_12' = 'Flavonifractor',  
                              'SRS1719577_6' = 'Faecalibacterium prausnitzii_D',  
                              'SRS294916_20' = 'Romboutsia timonensis', 
                              'SRS1719498_9' = 'Acetatifactor sp900066565',  
                              'ERS396297_11' = 'Lachnospira eligens_B',  
                              'ERS473255_26' = 'Faecalibacterium prausnitzii_G',  
                              'ERS608576_22' = 'Ruminococcus_E bromii_B',  
                              'ERS537328_30' = 'Faecalibacterium prausnitzii_K',  
                              'SRS143598_15' = 'Lachnospira sp000437735',  
                              'SRS476209_42' = 'Ruminococcus_D bicirculans',  
                              'ERS235603_16' = 'Agathobacter rectale',
                              'SRS1719112_8' = 'Oscillibacter sp900066435',
                              'ERS235531_43' = 'Faecalibacterium prausnitzii_F',
                              'SRS1735506_4' = 'Bacteroides ovatus',  
                              'SRS104400_110' = 'Lachnospira sp900316325',
                              'SRR4305229_bin.5' = 'Roseburia inulinivorans',  
                              'SRS103987_37' = 'ER4 sp000765235',  
                              'XieH_2016__YSZC12003_37172__bin.63' = 'Acutalibacter sp000435395',  
                              'LiJ_2014__O2.UC28-1__bin.61' = 'Ruminiclostridium_E siraeum',  
                              'LeChatelierE_2013__MH0074__bin.19' = 'CAG-45 sp900066395', 
                              'GCF_900036035.1_RGNV35913_genomic' = 'Faecalicatena gnavus', 
                              'GCF_000371685.1_Clos_bolt_90B3_V1_genomic' = 'Clostridium_M bolteae', 
                              'ERS608524_37' = 'Gemmiger formicilis',  
                              'LoombaR_2017__SID1050_bax__bin.11' = 'Anaeromassilibacillus sp002159845',
                              'SRS1735645_19' = 'Acetatifactor sp900066365', 
                              'ERS396519_11' = 'Lawsonibacter asaccharolyticus',
                              'SRR6028281_bin.3' = 'Lachnospira eligens_B', 
                              'VatanenT_2016__G80445__bin.9' = 'Faecalibacterium prausnitzii_D', 
                              'ZeeviD_2015__PNP_Main' = 'Blautia_A sp900066165', 
                              'GCF_000508885.1_ASM50888v1_genomic' = 'Flavonifractor sp000508885',
                              'SRS075078_49' = 'TF01-11 sp003529475',  
                              'LiSS_2016__FAT_DON_8-22-0-0__bin.28' = 'CAG-170 sp000432135', 
                              'NielsenHB_2014__MH0094__bin.44' = 'Prevotella copri',  
                              'SRR5127401_bin.3' = 'UBA11774 sp003507655', 
                              'ERS537235_19' = 'Bacteroides_B massiliensis',  
                              'SRR5558047_bin.10' = 'Alistipes putredinis', 
                              'VogtmannE_2016__MMRS43563715ST-27-0-0__bin.70' = 'CAG-81 sp900066785',  
                              'ERS537218_9' = 'Gemmiger sp003476825', 
                              'QinJ_2012__CON-091__bin.20' = 'Faecalibacterium prausnitzii_G', 
                              'GCF_001405615.1_13414_6_47_genomic' = 'Agathobacter faecis',  
                              'ERS235530_10' = 'CAG-1024 sp000432015')) +
  coord_flip()


pdf("panel_hum_hashes_var_imp_by_up_down_assembled.pdf", height = 10.5, width = 6)
ggarrange(plt_hash_count2, plt_var_imp2, nrow = 2,
          common.legend = T, legend = "bottom", heights = c(1, .9))
dev.off()


# count num abund hashes --------------------------------------------------
cd_up %>%
  select(hashval) %>%
  distinct() %>% 
  nrow()

cd_down  %>%
  select(hashval) %>%
  distinct() %>% 
  nrow()

uc_up  %>%
  select(hashval) %>%
  distinct() %>% 
  nrow()

uc_down  %>%
  select(hashval) %>%
  distinct() %>% 
  nrow()


# enrichment --------------------------------------------------------------

# download kegg dbs
pathways <- read_tsv("http://rest.kegg.jp/link/pathway/ko", col_names = c("KO", "path"))
pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

pathway_names <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))



cd_up_enriched<- enricher(gene = cd_up$KEGG_ko, 
                          TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(cd_up_enriched)

cd_down_enriched<- enricher(gene = cd_down$KEGG_ko, 
                          TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(cd_down_enriched)


uc_up_enriched<- enricher(gene = uc_up$KEGG_ko, 
                          TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(uc_up_enriched)

uc_down_enriched<- enricher(gene = unique(uc_down$KEGG_ko), 
                            TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(uc_down_enriched)

# per species enrichment ------------------------------------------------

# RUMINOCCUS GNAVUS
rgnv_cd_up <- cd_up %>%
  filter(gather_genome == "GCF_900036035.1_RGNV35913_genomic")
write.table(rgnv_cd_up, "sandbox/hash_corncob/GCF_900036035.1_RGNV35913_genomic_cd_up_hashes.txt", 
           col.names = F, row.names = F, quote = F)
rgnv_cd_up %>%
  select(hashval) %>%
  distinct()
rgnv_cd_up_enriched<- enricher(gene = rgnv_cd_up$KEGG_ko, 
                          TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(rgnv_cd_up_enriched)

rgnv_cd_down <- cd_down %>%
  filter(gather_genome == "GCF_900036035.1_RGNV35913_genomic")
rgnv_cd_down %>%
  select(hashval) %>%
  distinct() 
write.table(rgnv_cd_down, "sandbox/hash_corncob/GCF_900036035.1_RGNV35913_genomic_cd_down_hashes.txt", 
            col.names = F, row.names = F, quote = F)

rgnv_cd_down_enriched<- enricher(gene = rgnv_cd_down$KEGG_ko, 
                               TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(rgnv_cd_down_enriched)

## CLOS BOLT

cbolt_cd_up <- cd_up %>%
  filter(gather_genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic")

cbolt_cd_up_enriched<- enricher(gene = cbolt_cd_up$KEGG_ko, 
                               TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(cbolt_cd_up_enriched)

cbolt_cd_down <- cd_down %>%
  filter(gather_genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic")


cbolt_cd_down_enriched<- enricher(gene = cbolt_cd_down$KEGG_ko, 
                                 TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(cbolt_cd_down_enriched)

## Clostridiales bacterium VE202-03 (firmicutes) 	GCF_000508885.1_ASM50888v1_genomic
ve202_cd_up <- cd_up %>%
  filter(gather_genome == "GCF_000508885.1_ASM50888v1_genomic")


ve202_cd_up_enriched<- enricher(gene = ve202_cd_up$KEGG_ko, 
                                TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(ve202_cd_up_enriched)

ve202_cd_down <- cd_down %>%
  filter(gather_genome == "GCF_000508885.1_ASM50888v1_genomic")

ve202_cd_down_enriched<- enricher(gene = ve202_cd_down$KEGG_ko, 
                                  TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(ve202_cd_down_enriched)

ve202_uc_up <- uc_up %>%
  filter(gather_genome == "GCF_000508885.1_ASM50888v1_genomic")
ve202_uc_down <- uc_down %>%
  filter(gather_genome == "GCF_000508885.1_ASM50888v1_genomic")
# ERS537353_12

bin12_cd_up <- cd_up %>%
  filter(gather_genome == "ERS537353_12")

bin12_cd_up_enriched<- enricher(gene = bin12_cd_up$KEGG_ko, 
                                TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(bin12_cd_up_enriched)

bin12_cd_down <- cd_down %>%
  filter(gather_genome == "ERS537353_12")

bin12_cd_down_enriched<- enricher(gene = bin12_cd_down$KEGG_ko, 
                                  TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(bin12_cd_down_enriched)



bin12_uc_up <- uc_up %>%
  filter(gather_genome == "ERS537353_12")

bin12_uc_up_enriched<- enricher(gene = bin12_uc_up$KEGG_ko, 
                                TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(bin12_uc_up_enriched)

bin12_uc_down <- uc_down %>%
  filter(gather_genome == "ERS537353_12")

bin12_uc_down_enriched<- enricher(gene = bin12_uc_down$KEGG_ko, 
                                  TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(bin12_uc_down_enriched)



# write out hash sets -----------------------------------------------------
# hash sets just for CD
write.table(cbolt_cd_up$hashval, "sandbox/hash_corncob/GCF_000371685.1_Clos_bolt_90B3_V1_genomic_cd_up_hashes.txt", 
            col.names = F, row.names = F, quote = F)
write.table(cbolt_cd_down$hashval, "sandbox/hash_corncob/GCF_000371685.1_Clos_bolt_90B3_V1_genomic_cd_down_hashes.txt", 
            col.names = F, row.names = F, quote = F)
write.table(rgnv_cd_up$hashval, "sandbox/hash_corncob/GCF_900036035.1_RGNV35913_genomic_cd_up_hashes.txt", 
            col.names = F, row.names = F, quote = F)
write.table(rgnv_cd_down$hashval, "sandbox/hash_corncob/GCF_900036035.1_RGNV35913_genomic_cd_down_hashes.txt", 
            col.names = F, row.names = F, quote = F)

# hash sets for CD and UC
write.table(ve202_cd_up$hashval, "sandbox/hash_corncob/GCF_000508885.1_ASM50888v1_genomic_cd_up_hashes.txt", 
            col.names = F, row.names = F, quote = F)
write.table(ve202_cd_down$hashval, "sandbox/hash_corncob/GCF_000508885.1_ASM50888v1_genomic_cd_down_hashes.txt", 
            col.names = F, row.names = F, quote = F)
write.table(ve202_uc_up$hashval, "sandbox/hash_corncob/GCF_000508885.1_ASM50888v1_genomic_uc_up_hashes.txt", 
            col.names = F, row.names = F, quote = F)
write.table(ve202_uc_down$hashval, "sandbox/hash_corncob/GCF_000508885.1_ASM50888v1_genomic_uc_down_hashes.txt", 
            col.names = F, row.names = F, quote = F)


write.table(bin12_cd_down$hashval, "sandbox/hash_corncob/ERS537353_12_cd_down_hashes.txt", 
            col.names = F, row.names = F, quote = F)
write.table(bin12_cd_up$hashval, "sandbox/hash_corncob/ERS537353_12_cd_up_hashes.txt", 
            col.names = F, row.names = F, quote = F)

write.table(bin12_uc_down$hashval, "sandbox/hash_corncob/ERS537353_12_uc_down_hashes.txt", 
            col.names = F, row.names = F, quote = F)
write.table(bin12_uc_up$hashval, "sandbox/hash_corncob/ERS537353_12_uc_up_hashes.txt", 
            col.names = F, row.names = F, quote = F)

# plot per species enrichment ------------------------

# crohn's disease

cd_up_enriched_orgs <- bin12_cd_up_enriched@result %>%
  filter(qvalue < 0.05) %>%
  select(Description, GeneRatio, Count) %>%
  mutate(genome = "Flavonifractor")

cd_up_enriched_orgs <- cbolt_cd_up_enriched@result %>%
  filter(qvalue < 0.05) %>%
  select(Description, GeneRatio, Count) %>%
  mutate(genome = "Clostridium_M bolteae") %>%
  rbind(cd_up_enriched_orgs)

cd_up_enriched_orgs <- ve202_cd_up_enriched@result %>%
  filter(qvalue < 0.05) %>%
  select(Description, GeneRatio, Count) %>%
  mutate(genome = "Flavonifractor sp000508885")%>%
  rbind(cd_up_enriched_orgs)


cd_up_enriched_orgs <- rgnv_cd_up_enriched@result %>%
  filter(qvalue < 0.05) %>%
  select(Description, GeneRatio, Count) %>%
  mutate(genome = "Faecalicatena gnavus")%>%
  rbind(cd_up_enriched_orgs)

common_descriptions <- cd_up_enriched_orgs %>% 
  group_by(Description) %>% 
  tally() %>% 
  filter(n >= 2) %>% 
  ungroup()

cd_up_enriched_orgs <- cd_up_enriched_orgs %>% 
  filter(Description %in% common_descriptions$Description) %>%
  separate(GeneRatio, into = c("Count", "Total"), sep = "\\/") %>%
  mutate(GeneRatio = as.numeric(Count)/as.numeric(Total))

ggplot(cd_up_enriched_orgs, aes(x = Description, y = GeneRatio, color = genome, size = Count)) +
  geom_point(alpha = .5) +
  theme_minimal() +
  coord_flip()+
  theme(legend.position = "right",
        axis.title.y = element_blank()) +
  labs(x = "gene ratio") + 
  guides(colour = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values = c("#BEBADA", "#FCCDE5", "#D9D9D9", "#BC80BD"),
                     labels = c("Clostrium bolteae", "Faecalicatena gnavus", 
                                "Flavonifractor plautii", "Flavonifractor sp000508885"))


bin12_cd_down_enriched
ve202_cd_down_enriched
cbolt_cd_down_enriched
rgnv_cd_down_enriched



bin12_uc_up_enriched
bin12_uc_down_enriched

