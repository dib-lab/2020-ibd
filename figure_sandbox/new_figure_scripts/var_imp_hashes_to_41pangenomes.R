setwd("~/github/2020-ibd")

library(dplyr)
library(readr)
library(ggplot2)
library(ranger)
library(yarrr)

# the purpose of this script is to plot the proportion of each of predictive 
# genome from vita variable selection based on when they appear
# in importance rankings. i.e. only 1 genome contains the most
# predictive hash

# read in ranger and get variable importance

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
  arrange(desc(total_imp))

# scale_fill_manual(values = c("#D57A6DFF", "#E8B762FF", "#9CCDDFFF", "#525052FF", "#E6CEAFFF", "#BA9570FF"))


# plot genome importance for 41 pangenomes -------------------------------

hash_to_pangenome <- read_csv("outputs/sgc_pangenome_gather/hash_to_genome_map_at_least_5_studies_pangenome.csv") %>%
  select(-X1) %>%
  mutate(sig = gsub("_renamed\\.sig", "", sig)) %>%
  mutate(sig = gsub("outputs\\/sgc_pangenome_sigs\\/", "", sig))
varimp_cum$hash <- as.numeric(as.character(varimp_cum$hash))

all <- left_join(hash_to_pangenome, varimp_cum, by = 'hash')

# summarize to total importance by pangenome
by_pangenome <- all %>%
  select(sig, hash, total_imp) %>%
  distinct() %>%
  group_by(sig) %>%
  summarize(pangenome_imp = sum(total_imp, na.rm = T)) %>%
  mutate(pangenome_imp_div_6 = pangenome_imp / 6)

# add in annotations
lca <- read_csv("inputs/at_least_5_studies_vita_vars_gather_all_lca.csv")

by_pangenome_labelled <- left_join(by_pangenome, lca, by = c("sig" = "name_no_gz"))
colnames(by_pangenome_labelled) <- c("genome", "pangenome_imp", "pangenome_importance", "genome2", "GTDB", 
                                  "NCBI", "contaminated_with", "contaminated")                
by_pangenome_labelled <- by_pangenome_labelled %>%
  select(genome, pangenome_importance, GTDB, NCBI, contaminated, contaminated_with) %>%
  arrange(desc(pangenome_importance)) %>%
  mutate(rank = 1:nrow(.))
write_tsv(by_pangenome_labelled, "~/Desktop/41pangenomes_gather_matches_ranked.tsv")

ggplot(by_pangenome, aes(x = reorder(sig, pangenome_imp_div_6), y = pangenome_imp_div_6)) +
  geom_col() +
  theme_minimal() +
  coord_flip()

# summarize by study
by_pangenome_study <- all %>%
  filter(!is.na(study)) %>%
  select(sig, hash, norm, study) %>%
  distinct() %>%
  group_by(sig, study) %>%
  summarize(pangenome_imp_study = sum(norm, na.rm = T)) %>%
  mutate(pangenome_imp_study_div_6 = pangenome_imp_study / 6)

ggplot(by_pangenome_study, aes(x = reorder(sig, by_pangenome_study$pangenome_imp_study_div_6), 
                            y = pangenome_imp_study_div_6, fill = study)) +
  geom_col() +
  theme_minimal() +
  coord_flip() +
  labs(y = "fraction of total importance", x ="pangenome") +
  scale_fill_manual(values = c("#D57A6DFF", "#E8B762FF", "#9CCDDFFF", "#525052FF", "#E6CEAFFF", "#BA9570FF"))



by_pangenome_study_labeled <- left_join(by_pangenome_study, lca, by = c("sig" = "name_no_gz")) %>%
  arrange(desc(pangenome_imp_study_div_6))

a <- ifelse(by_pangenome_study_labeled$`Contaminated? (GTDB)` == "yes", "red", "black")
ggplot(by_pangenome_study_labeled, aes(x = reorder(GTDB, by_pangenome_study_labeled$pangenome_imp_study_div_6), 
                                    y = pangenome_imp_study_div_6, fill = study)) +
  geom_col() +
  theme_minimal() +
  coord_flip() +
  labs(y = "fraction of total importance", x = "genome") +
  theme(axis.text.y = element_text(colour = a)) +
  scale_fill_manual(values = c("#D57A6DFF", "#E8B762FF", "#9CCDDFFF", "#525052FF", "#E6CEAFFF", "#BA9570FF"))
