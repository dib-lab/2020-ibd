setwd("~/github/2020-ibd")
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(ComplexUpset)

files_all <- list.files("outputs/gather", "vita_vars_all.csv$", full.names = T)
files_all
all <- files_all %>%
  set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/gather\\/", "", source)) %>%
  mutate(source = gsub("_vita_vars_all\\.csv", "", source))

iHMP <- all %>% filter(source == "iHMP")
iHMP <- as.character(iHMP$name)

PRJEB2054 <- all %>% filter(source == "PRJEB2054")
PRJEB2054 <- as.character(PRJEB2054$name)

PRJNA237362 <- all %>% filter(source == "PRJNA237362")
PRJNA237362 <- as.character(PRJNA237362$name)

PRJNA237362 <- all %>% filter(source == "PRJNA237362")
PRJNA237362 <- as.character(PRJNA237362$name)

PRJNA385949 <- all %>% filter(source == "PRJNA385949")
PRJNA385949 <- as.character(PRJNA385949$name)

PRJNA400072  <- all %>% filter(source == "PRJNA400072")
PRJNA400072  <- as.character(PRJNA400072$name)

SRP057027  <- all %>% filter(source == "SRP057027")
SRP057027  <- as.character(SRP057027$name)

upset_df <- UpSetR::fromList(list(iHMP = iHMP,
                          PRJEB2054 = PRJEB2054,
                          PRJNA237362 = PRJNA237362,
                          PRJNA385949 = PRJNA385949,
                          PRJNA400072 = PRJNA400072,
                          SRP057027 = SRP057027))

upset_df$genome <- unique(c(iHMP, PRJEB2054, PRJNA237362, PRJNA385949, 
                          PRJNA400072,  SRP057027))
# remove folder directory structure from genome name
upset_df$genome <- gsub(".*\\/", "", upset_df$genome)

# read in variable importance, summarize to variable importance per genome
# (cummulative over study), and add to upset dataframe


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
varimp_cum$hash <- as.numeric(as.character(varimp_cum$hash))

# hash to genome map ------------------------------------------------------
hash_to_genome_mapping <- read_csv("outputs/gather_matches/vita_hash_to_genome_mapping.csv")
hash_to_genome_mapping$genome <-gsub(".*\\/", "", hash_to_genome_mapping$genome)
hash_to_genome_mapping$genome <-gsub("\\.sig", "", hash_to_genome_mapping$genome)

hash_to_genome_mapping <- filter(hash_to_genome_mapping, !is.na(genome))

all <- left_join(varimp_cum, hash_to_genome_mapping, by = c("hash" = "index"))

genome_imp <- all %>%
  group_by(genome) %>%
  summarize(genome_varimp_cum = sum(norm)) %>%
  arrange(desc(genome_varimp_cum)) 

# check 1:1 match
table(upset_df$genome %in% genome_imp$genome)

upset_df <- left_join(upset_df, genome_imp)
write_tsv(upset_df, "tmp_upset_df.tsv")

# try upset ---------------------------------------------------------------
library(readr)
upset_df<- read_tsv("tmp_upset_df.tsv")
upset_df <- upset_df %>%
  select(genome, iHMP, PRJEB2054, PRJNA237362, PRJNA385949, 
         PRJNA400072, SRP057027, genome_varimp_cum)
studies = c("iHMP", "PRJEB2054", "PRJNA237362", "PRJNA385949", "PRJNA400072", "SRP057027")
upset(upset_df, studies, sort_intersections_by='degree',
      dot_size = 2, set_sizes=FALSE, name=NULL, width_ratio = 0.1,
      base_annotations=list('k-mers'=intersection_size(counts=FALSE)),            
      annotations = list('cumulative importance'=list(
        aes=aes(x=intersection, genome_varimp_cum/6),
        geom=list(geom_bar(stat='identity')))))


# upset small -------------------------------------------------------------
upset_df_small <- upset_df[rowSums(upset_df[ , 1:6]) >= 5, ]

upset(upset_df_small, studies, sort_intersections_by='degree',
      dot_size = 2, set_sizes=FALSE, name=NULL, width_ratio = 0.1,
      base_annotations=list('genomes'=intersection_size(counts=FALSE)),            
      annotations = list('cumulative importance'=list(
        aes=aes(x=intersection, genome_varimp_cum/6),
        geom=list(geom_bar(stat='identity')))))

sum(upset_df_small$genome_varimp_cum) / 6
