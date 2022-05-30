library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

sgc <- read_csv("outputs/gather_matches_loso_multifasta/all-multifasta-query-results.csv") %>%
  filter(hashval != "hashval") %>%
  separate(record_name, sep = " ", into = c("sequence_name", "prokka_annotation"), extra = "merge") %>%
  mutate(hashval = as.character(hashval))
length(unique(sgc$hashval))
length(unique(sgc$hashval))/3859



# read in var imp ---------------------------------------------------------

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
head(varimp_cum)

varimp_cum_sgc <- varimp_cum %>%
  filter(hash %in% sgc$hashval)
sum(varimp_cum_sgc$total_imp)
sum(varimp_cum$total_imp)
dim(varimp_cum_sgc)

hash_info <- left_join(sgc, varimp_cum, by = c("hashval" = "hash"))
length(unique(hash_info$hashval))
length(hash_info$hashval)
length(hash_info$hashval)
sum(hash_info$total_imp) # variable imp held by hashes after multifasta query

# read in eggnog ----------------------------------------------------------
# read in eggnog results
eggnog <- read_tsv("outputs/gather_matches_loso_multifasta/all-multifasta-query-results.emapper.annotations", 
                   comment = "#", 
                   col_names = c('query_name', 'seed_eggNOG_ortholog',	
                                 'seed_ortholog_evalue', 'seed_ortholog_score',
                                 'best_tax_level', 'Preferred_name', 'GOs', 'EC',
                                 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',	'KEGG_TC',
                                 'CAZy', 'BiGG_Reaction', 'taxonomic_scope', 
                                 'eggNOG_OGs', 'best_eggNOG_OG', 'COG_functional_cat',	
                                 'eggNOG', 'free_text_desc'))
eggnog$KEGG_ko <- gsub(",.*", "", eggnog$KEGG_ko) # select first KO annotation

hash_info <- left_join(hash_info, eggnog, by = c("sequence_name" = "query_name"))


# read in corncob ---------------------------------------------------------

corncob <- read_tsv("sandbox/hash_corncob/sig_ccs_hashes.tsv") %>%
  mutate(hashval = as.character(aa_seq)) %>%
  select(-aa_seq)

hash_info <- left_join(hash_info, corncob, by = "hashval")

# designate marker genes -------------------------------------------------

# sources: singlem, checkm, 
# https://link.springer.com/content/pdf/10.1007/s12275-018-8014-6.pdf
marker_genes <- c("rpsK", "rpsQ", "rpsB", "rpsH", "rpsP", "rpsD", "rpsO", "rpsT", "rpsS", 
                  "rpsE", "rpsC", "rpsI", "rpsJ", "rpsU", "rpsR", "rpsF", "rpsN", "rpsM", 
                  "rpsL", "rpsG", "rplU", "rplQ", "rplP", "rplS", "rplL", "rplD", "rplB", 
                  "rplC", "rplW", "rplE", "rplV", "rplI", "rplM", "rplA", "rplX", "rplJ", "rplR",
                  "rplN", "rplO", "rplT", "rpoB", "rpoA", "rpoC", "rpoN", "rpoD",  "rpmC",  
                  "rpmJ",  "rpmB",  "rpmA",  "rpmE2", "rpmD", "gyrA", "gyrB", "ispH", "rsgA",
                  "lepA", "yqgF", "nusA", "nusG", "uvrC", "secA", "secY", "secG", "recR", "recA",
                  "rumA", "rumA1", "spoU", "rrmJ", "rrmA", "infB", "alaS", "argS", 
                  "ileS", "ctgA", "coaE", "cysS", "dnaA", "dnaG", "dnaX", "engA", "ffh", 
                  "frr", "ftsY", "gmk", "hisS", "infC", "leuS", "ligA", 
                  "pgk", "pheS", "pheT", "prfA", "pyrG", "rbfA", "rnc", "serS", "tsaD",
                  "uvrB", "ybeY", "ychF", "aspS", "ksgA", "fmt", "tilS", "tsf", "trmD",
                  "tig", "smpB", "truB")


ribo <- c("16S ribosomal 23S ribosomal 50S ribosomal")

hash_info <- hash_info %>%
  mutate(marker =  ifelse(str_detect(string = prokka_annotation, pattern = "ribosomal"), "marker", 
                          ifelse(Preferred_name %in% marker_genes, "marker", "other"))) 

hash_info <- hash_info %>%
  mutate(direction =  ifelse(estimate < 0, "down", "up")) 

hash_info_simple <- hash_info %>%
  select(hashval, total_imp, marker, direction) %>%
  distinct()

ggplot(hash_info_simple, aes(x = marker, y = total_imp)) +
  geom_col() +
  theme_minimal() +
  labs(x = "genetic element", y = "cumulative importance") + coord_flip()


table(hash_info$marker)
hash_info_simple %>%
  group_by(marker) %>%
  summarise(sum(total_imp))
nrow(hash_info_simple)
length(unique(hash_info_simple$hashval))/3859

View(hash_info_simple)
table(hash_info_simple$marker, hash_info_simple$direction)
# add assembly vs. not ----------------------------------------------------

unassembled_hashes <- read_csv("sandbox/test_megahit_diginorm_nocat/megahit_gather/at_least_5_vita_vars_vs_megahit_assemblies_tbp0.un.csv") %>%
  mutate(hashval = as.character(minhash)) %>%
  select(hashval) %>%
  mutate(assembled = "unassembled") 

hash_info <- left_join(hash_info, unassembled_hashes, by = "hashval")
hash_info$assembled <- ifelse(is.na(hash_info$assembled), "assembled", hash_info$assembled)



hash_info_simple <- hash_info %>%
  select(hashval, total_imp, marker, assembled) %>%
  distinct()

ggplot(hash_info_simple, aes(x = marker, y = total_imp, fill = assembled)) +
  geom_col() +
  theme_minimal() +
  labs(x = "genetic element", y = "cumulative importance") +
  scale_fill_manual(values = c(assembled = "#dddddd", unassembled = "#858585"))
  
