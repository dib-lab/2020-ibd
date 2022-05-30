setwd("~/github/2020-ibd")

library(dplyr)
library(readr)
library(tidyr)
library(ggthemes)
library(ggtree)
# library(tidytree) 
# library(treeio)
library(ape)
# library(phytools)
#library(ggnewscale)
#library(ggstance)
#library(cowplot)

# read the phylogenetic tree ---------------------------------------------
tree <- read.tree("sandbox/gtdbtk_41_genomes/gtdbtk.bac120.classify.tree")

# drop tips from GTDB -- otherwise is a very large tree
drop_tips <- tree$tip.label[grepl(pattern = "^GB_", tree$tip.label)]
tree <- drop.tip(tree, drop_tips)
drop_tips <- tree$tip.label[grepl(pattern = "^RS_", tree$tip.label)]
tree <- drop.tip(tree, drop_tips)

# read in and format taxonomy ---------------------------------------------

# Read in GTDB taxonomy and reorder to match tree

gtdb <- read_tsv("sandbox/gtdbtk_41_genomes/gtdbtk.bac120.summary.tsv") %>%
  mutate(accession = user_genome) %>%
  select(accession, classification) %>%
  mutate(classification = gsub("[dkpcofgs]__", "", classification)) %>%
  separate(data = ., col = classification, sep = ";",
           into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
  mutate(strain = NA) %>%
  mutate(source = "gtdbtk")%>%
  #transmute(labels = ifelse(species == "", genus, species), across(everything()))
  mutate(labels = ifelse(species == "", genus, species))

# node names are stored in T$tip.label.
gtdb <- gtdb[order(match(gtdb$accession, tree$tip.label)), ] # sort Tax according to the tree
stopifnot(all.equal(gtdb$accession, tree$tip.label)) # make sure the sorting worked
# tree$tip.label <- gtdb$labels

# prepare to label clades
ruminococcaceae <- findMRCA(tree, c("Ruminococcus_E bromii_B", "Gemmiger sp003476825"))
oscillospiraceae <- findMRCA(tree, c("CAG-170", "Flavonifractor sp000508885"))
lachnospiraceae <- findMRCA(tree, c("Agathobacter rectale", "Clostridium_M bolteae"))
bacteroidaceae <- findMRCA(tree, c("Prevotella copri", "Bacteroides_B massiliensis"))

pdf("new_figure_scripts/gtdbtk_41_genome_tree.pdf", height = 6.5, width = 3)
ggtree(tree)+
  geom_tiplab(size = 2) +
  geom_cladelabel(node = ruminococcaceae, barsize = 1, label = "Ruminococcaceae", 
                  offset = 1.35, angle = 90, hjust = .5, vjust = 1)+
  geom_cladelabel(node = oscillospiraceae, barsize = 1, label = "Oscillospiraceae", 
                  offset = 1.3, angle = 90, hjust = .5, vjust = 1) +
  geom_cladelabel(node = lachnospiraceae, barsize = 1, label = "Lachnospiraceae", 
                  offset = 1.3, angle = 90, hjust = .5, vjust = 1) +
  geom_cladelabel(node = bacteroidaceae, barsize = 1, label = "Bacteroidaceae", 
                  offset = 1.1, angle = 90, hjust = .3, vjust = 1) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()


# read in aux data --------------------------------------------------------

library(tidyr)

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

# # read in assembled/unassembled information
# unassembled_hashes <- read_csv("sandbox/test_megahit_diginorm_nocat/megahit_gather/at_least_5_vita_vars_vs_megahit_assemblies_tbp0.un.csv") %>%
#   mutate(hashval = as.character(minhash)) %>%
#   select(hashval) %>%
#   mutate(assembled = "unassembled") 
# 
# hash_results <- left_join(hash_results, unassembled_hashes, by = "hashval")
# hash_results$assembled <- ifelse(is.na(hash_results$assembled), "assembled", hash_results$assembled)
# 

# remove protein numbers to generate a gather genome column
hash_results <- hash_results %>%
  mutate(gather_genome = gsub(".fna.*", "", sequence_name)) %>%
  mutate(gather_genome = gsub(".fa.*", "", gather_genome))%>%
  mutate(gather_genome = gsub("_232_.*", "", gather_genome)) %>%
  mutate(gather_genome = gsub("ZeeviD_2015__PNP_Main", "ZeeviD_2015__PNP_Main_232__bin.27", gather_genome)) %>%
  mutate(direction = ifelse(estimate > 0, "up", "down")) %>%
  mutate(set = gsub("mu\\.diagnosis", "", paste0(mu, "_", direction)))

hash_results_small <- hash_results %>%
  select(hashval, mu, estimate, gather_genome,  direction) %>%
  distinct()


# add data -----------------------------------------------------------------
hash_results_small$gather_genome <- as.factor(hash_results_small$gather_genome)
kmers_CD <- hash_results_small %>%
  filter(mu == "mu.diagnosisCD") %>%
  group_by(gather_genome, direction) %>% 
  tally() %>%
  filter(!is.na(gather_genome)) %>%
  left_join(gtdb, by = c("gather_genome" = "accession"))
kmers_CD <- kmers_CD[order(match(kmers_CD$gather_genome, gtdb$accession)), ]

plot_kmers_cd <- ggplot(kmers_CD, aes(x = gather_genome, y = n, fill = direction)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(hjust=0.5)) + 
  coord_flip() +
  labs(y = "number of k-mers", title = "CD") +
  scale_fill_manual(values = c("down" = "#DCDCDC", "up" = "#808080"), drop = F) +
  scale_x_discrete(drop = F,
                   labels = c('ERS537353_12' = 'Flavonifractor plautti',  
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
                              'ZeeviD_2015__PNP_Main_232__bin.27' = 'Blautia_A sp900066165', 
                              'GCF_000508885.1_ASM50888v1_genomic' = 'Flavonifractor sp000508885',
                              'SRS075078_49' = 'TF01-11 sp003529475',  
                              'LiSS_2016__FAT_DON_8-22-0-0__bin.28' = 'CAG-170', 
                              'NielsenHB_2014__MH0094__bin.44' = 'Prevotella copri',  
                              'SRR5127401_bin.3' = 'UBA11774 sp003507655', 
                              'ERS537235_19' = 'Bacteroides_B massiliensis',  
                              'SRR5558047_bin.10' = 'Alistipes putredinis', 
                              'VogtmannE_2016__MMRS43563715ST-27-0-0__bin.70' = 'CAG-81 sp900066785',  
                              'ERS537218_9' = 'Gemmiger sp003476825', 
                              'QinJ_2012__CON-091__bin.20' = 'Faecalibacterium prausnitzii_G', 
                              'GCF_001405615.1_13414_6_47_genomic' = 'Agathobacter faecis',  
                              'ERS235530_10' = 'CAG-1024 sp000432015'))
plot_kmers_cd
kmers_UC <- hash_results_small %>%
  filter(mu == "mu.diagnosisUC") %>%
  group_by(gather_genome, direction) %>% 
  tally() %>%
  filter(!is.na(gather_genome)) %>%
  left_join(gtdb, by = c("gather_genome" = "accession"))
# back add dropped level
tmp <- kmers_CD[71, ] %>%
  mutate(n = sum(n) - sum(n))
kmers_UC <- rbind(kmers_UC, tmp)

plot_kmers_uc <- ggplot(kmers_UC, aes(x = gather_genome, y = n, fill = direction)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(hjust=0.5),
        legend.position = "bottom") + 
  coord_flip() +
  labs(y = "number of k-mers", title = "UC") +
  ylim(0, 350) +
  scale_fill_manual(values = c("down" = "#DCDCDC", "up" = "#808080"), drop = F) +
  scale_x_discrete(drop = F, labels = c('ERS537353_12' = 'Flavonifractor plautti',  
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
                              'ZeeviD_2015__PNP_Main_232__bin.27' = 'Blautia_A sp900066165', 
                              'GCF_000508885.1_ASM50888v1_genomic' = 'Flavonifractor sp000508885',
                              'SRS075078_49' = 'TF01-11 sp003529475',  
                              'LiSS_2016__FAT_DON_8-22-0-0__bin.28' = 'CAG-170', 
                              'NielsenHB_2014__MH0094__bin.44' = 'Prevotella copri',  
                              'SRR5127401_bin.3' = 'UBA11774 sp003507655', 
                              'ERS537235_19' = 'Bacteroides_B massiliensis',  
                              'SRR5558047_bin.10' = 'Alistipes putredinis', 
                              'VogtmannE_2016__MMRS43563715ST-27-0-0__bin.70' = 'CAG-81 sp900066785',  
                              'ERS537218_9' = 'Gemmiger sp003476825', 
                              'QinJ_2012__CON-091__bin.20' = 'Faecalibacterium prausnitzii_G', 
                              'GCF_001405615.1_13414_6_47_genomic' = 'Agathobacter faecis',  
                              'ERS235530_10' = 'CAG-1024 sp000432015'))


plot_kmers_uc

library(ggplot2)
library(aplot)
# color tree by family
ruminococcaceae <- phytools::findMRCA(tree, c("SRS476209_42", "ERS537218_9"))
oscillospiraceae <- phytools::findMRCA(tree, c("LiSS_2016__FAT_DON_8-22-0-0__bin.28", "GCF_000508885.1_ASM50888v1_genomic"))
lachnospiraceae <- phytools::findMRCA(tree, c("ERS235603_16", "GCF_000371685.1_Clos_bolt_90B3_V1_genomic"))
bacteroidaceae <- phytools::findMRCA(tree, c("NielsenHB_2014__MH0094__bin.44", "ERS537235_19"))
acutalibacteraceae <- phytools::findMRCA(tree, c("ERS608576_22", "XieH_2016__YSZC12003_37172__bin.63"))


plot_tree <- ggtree(tree) %<+% gtdb + 
  geom_tippoint(aes(color = family), size = 0) +
  geom_hilight(ruminococcaceae, aes(fill =  "darkorchid4")) +
  geom_hilight(oscillospiraceae, fill = "chartreuse4") +
  geom_hilight(lachnospiraceae, fill = "steelblue") +
  geom_hilight(acutalibacteraceae, fill = "darkgoldenrod4") +  
  geom_hilight(bacteroidaceae, fill = "pink") +
  guides(color = guide_legend(override.aes = list(size=8, alpha = .5))) +
  scale_color_manual(name="family",
                    values=c(Ruminococcaceae="darkorchid4",
                             Oscillospiraceae="chartreuse4",
                             Lachnospiraceae="steelblue",
                             Acutalibacteraceae = "darkgoldenrod4",
                             Bacteroidaceae = "pink",
                             `CAG-138`= "white",
                             Peptostreptococcaceae = "white",
                             Rikenellaceae ="white"), 
                    breaks = c("Ruminococcaceae", "Oscillospiraceae", "Lachnospiraceae",
                               "Acutalibacteraceae", "Bacteroidaceae")) +
  theme(legend.text = element_text(face = "italic"),
        legend.position = "bottom")

plot_tree
plot_kmers_cd %>% insert_left(plot_tree)
plot_kmers_uc %>% insert_left(plot_tree)


