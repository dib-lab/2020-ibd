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
  
varimp_cum_top_hashes <- head(unique(varimp_cum$hash), n = 20)

# almost 8000 hashes have a variable importance of 0 across all studies, 
# meaning they are not important. 
pdf(file = "new_figure_scripts/var_imp_hashes.pdf", width = 6, height = 4)
ggplot(varimp_cum %>%
         filter(hash %in% varimp_cum_top_hashes),
       aes(x = reorder(hash, total_imp), y = norm/6, fill = study)) +
  geom_col() +
  theme_minimal() +
  labs(x = "hash", y = "normalized importance (sums to 1)", title = "Twenty hashes with greatest cumulative normalized variable importance") +
  coord_flip() +
  scale_fill_manual(values = c("#D57A6DFF", "#E8B762FF", "#9CCDDFFF", "#525052FF", "#E6CEAFFF", "#BA9570FF")) +
  theme(plot.title.position = "plot")
dev.off()
# piratepal(palette = "decision")

# test how much var imp is contained in the shared hashes -----------------

tmp <- varimp_cum %>%
  filter(hash %in% at_least_5_hashes) # this vector is created in count_intersect_var_imp_hashes.R
sum(tmp$norm)
(sum(tmp$norm))/6

tmp <- varimp_cum %>%
  filter(hash %in% intersect_hashes_all) # this vector is created in count_intersect_var_imp_hashes.R
sum(tmp$norm)
sum(tmp$norm) / 6
# hash to genome map ------------------------------------------------------

varimp_cum$hash <- as.numeric(as.character(varimp_cum$hash))

hash_to_genome_mapping <- read_csv("outputs/gather_matches/vita_hash_to_genome_mapping.csv")
hash_to_genome_mapping$genome <- gsub("outputs\\/gather_matches\\/", "", hash_to_genome_mapping$genome)
hash_to_genome_mapping$genome <- gsub("\\.fna\\.sig", "", hash_to_genome_mapping$genome)
hash_to_genome_mapping$genome <- gsub(".sig", "", hash_to_genome_mapping$genome)
hash_to_genome_mapping$genome <- gsub(", whole genome shotgun sequence", "", hash_to_genome_mapping$genome)
hash_to_genome_mapping$genome <- gsub("\\.fa\\.sig", "", hash_to_genome_mapping$genome)
hash_to_genome_mapping$genome <- gsub("\\.fa\\.gz\\.sig", "", hash_to_genome_mapping$genome)
hash_to_genome_mapping <- filter(hash_to_genome_mapping, !is.na(genome))

all <- left_join(varimp_cum, hash_to_genome_mapping, by = c("hash" = "index"))

genome_study_imp <- all %>%
  group_by(genome, study) %>%
  summarize(study_genome_varimp_cum = sum(norm)) %>%
  arrange(desc(study_genome_varimp_cum)) 

genome_imp <- all %>%
  group_by(genome) %>%
  summarize(genome_varimp_cum = sum(norm)) %>%
  arrange(desc(genome_varimp_cum)) 

imp <- left_join(genome_study_imp, genome_imp)

# get the XX genomes with the highest cumulative variable importance
top_genomes <- imp %>% 
  select(genome, genome_varimp_cum) %>%
  distinct() %>%
  arrange(desc(genome_varimp_cum)) %>%
  head(n = 50)

pdf(file = "new_figure_scripts/var_imp_genomes.pdf", width = 11, height = 6)
ggplot(imp %>%
         filter(genome %in% top_genomes$genome), 
       aes(x = reorder(genome, genome_varimp_cum), y = study_genome_varimp_cum, fill = study)) +
  geom_col() +
  theme_minimal() +
  coord_flip() +
  scale_fill_manual(values = c("#D57A6DFF", "#E8B762FF", "#9CCDDFFF", "#525052FF", "#E6CEAFFF", "#BA9570FF")) +
  labs(y = "cumulative normalized genome importance (sums to 6)", x = "top 50 genomes")
dev.off()


# test how much var imp is contained in the shared genomes -----------------
source("new_figure_scripts/plot_gather_output_var_imp.R")
at_least_5_genomes <- gsub(".*\\/","", at_least_5_genomes) # this is from plot_gather_output_var_imp.R script
at_least_5_genomes <- gsub(".fna","", at_least_5_genomes)
at_least_5_genomes <- gsub(", whole genome shotgun sequence","", at_least_5_genomes)


table(at_least_5_genomes %in% imp$genome)
at_least_5_genomes[!at_least_5_genomes %in% imp$genome]
at_least_5_genomes_df <- imp %>%
  filter(genome %in% at_least_5_genomes) %>% # this vector is created in count_intersect_var_imp_hashes.R
  select(genome, genome_varimp_cum) %>%
  distinct()
sum(at_least_5_genomes_df$genome_varimp_cum) / 6
# 4 out of 6 var imp contained in the shared 41 genomes