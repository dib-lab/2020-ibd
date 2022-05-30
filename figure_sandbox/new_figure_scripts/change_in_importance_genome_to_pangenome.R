setwd("~/github/2020-ibd")

library(dplyr)
library(readr)
library(ggplot2)
library(ranger)
library(yarrr)
library(ggrepel)

# the purpose of this script is to plot the change in variable importance
# held by a genome after it becomes a pangenome via spacegraphcats queries

# read in ranger and get variable importance --------------------------------
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

# generate data from of importance 41 pangenomes -------------------------------

hash_to_pangenome <- read_csv("outputs/sgc_pangenome_gather/hash_to_genome_map_at_least_5_studies_pangenome.csv") %>%
  select(-X1) %>%
  mutate(sig = gsub("_renamed\\.sig", "", sig)) %>%
  mutate(sig = gsub("outputs\\/sgc_pangenome_sigs\\/", "", sig))

varimp_cum$hash <- as.numeric(as.character(varimp_cum$hash))

pangenome_importance <- left_join(hash_to_pangenome, varimp_cum, by = 'hash')

# summarize to total importance by pangenome
by_pangenome <- pangenome_importance %>%
  select(sig, hash, total_imp) %>%
  distinct() %>%
  group_by(sig) %>%
  summarize(pangenome_imp = sum(total_imp, na.rm = T)) %>%
  mutate(pangenome_imp_div_6 = pangenome_imp / 6)

# generate data from of importance 41 genomes -------------------------------

hash_to_genome <- read_csv("sandbox/41genome_sigs/41-hash-to-gather-match.csv") %>%
  select(-X1) %>%
  mutate(sig = gsub("\\.sig", "", sig))

varimp_cum$hash <- as.numeric(as.character(varimp_cum$hash))

genome_importance <- left_join(hash_to_genome, varimp_cum, by = 'hash')

# summarize to total importance by genome
by_genome <- genome_importance %>%
  select(sig, hash, total_imp) %>%
  distinct() %>%
  group_by(sig) %>%
  summarize(genome_imp = sum(total_imp, na.rm = T)) %>%
  mutate(genome_imp_div_6 = genome_imp / 6)



# join and label and plot  ----------------------------------------------------
lca <- read_csv("~/Desktop/at_least_5_studies_vita_vars_gather_all_lca.csv")
by_pangenome_labelled <- left_join(by_pangenome, lca, by = c("sig" = "name_no_gz"))
by_genome_labelled <- left_join(by_genome, lca, by = c("sig" = "name"))

all_labelled <- left_join(by_genome_labelled, by_pangenome_labelled, by = c("sig" = "name"))

# all log-log plot
ggplot() +
  geom_point(data = all_labelled, aes(x = genome_imp_div_6, y = pangenome_imp_div_6), size = 3, alpha = .5) +
  theme_minimal() +
  #xlim(c(0, .13)) + ylim(c(0, .13))+
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  labs(x = "log cumulative importance of genome", y = "log cumulative importance of pangenome",
       title = "Log cumulative importance of genomes before and after spacegraphcats queries") +
  theme(plot.title.position = "plot",
        panel.grid.major = element_line(color = "grey25"),
        panel.grid.minor = element_line(color = "grey25")) +
  geom_label_repel(data = subset(all_labelled, pangenome_imp_div_6 > .05), 
                   aes(x = genome_imp_div_6, y = pangenome_imp_div_6, label = GTDB.x), force =5) +
  geom_label_repel(data = subset(all_labelled, genome_imp_div_6 > .045), 
                   aes(x = genome_imp_div_6, y = pangenome_imp_div_6, label = GTDB.x), force = 5)

# no log log plot
ggplot() +
  geom_point(data = all_labelled, aes(x = genome_imp_div_6, y = pangenome_imp_div_6), size = 3, alpha = .5) +
  theme_minimal() +
  #xlim(c(0, .13)) + ylim(c(0, .13))+
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey25") +
  labs(x = "cumulative importance of genome", y = "cumulative importance of pangenome",
       title = "Cumulative importance of genomes before and after spacegraphcats queries") +
  theme(plot.title.position = "plot",
        panel.grid.major = element_line(color = "grey25"),
        panel.grid.minor = element_line(color = "grey25")) +
  geom_label_repel(data = subset(all_labelled, pangenome_imp_div_6 > .05), 
                   aes(x = genome_imp_div_6, y = pangenome_imp_div_6, label = GTDB.x), force =5) +
  geom_label_repel(data = subset(all_labelled, genome_imp_div_6 > .045), 
                   aes(x = genome_imp_div_6, y = pangenome_imp_div_6, label = GTDB.x), force = 5)


