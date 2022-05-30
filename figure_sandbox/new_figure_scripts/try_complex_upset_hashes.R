setwd("~/github/2020-ibd")
library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(purrr)
library(readr)

prjeb2054 <- read_tsv("outputs/vita_rf/PRJEB2054_vita_vars.txt", col_names = "hash")
prjna237362 <- read_tsv("outputs/vita_rf/PRJNA237362_vita_vars.txt", col_names = "hash")
ihmp <- read_tsv("outputs/vita_rf/iHMP_vita_vars.txt", col_names = "hash")
srp057027 <- read_tsv("outputs/vita_rf/SRP057027_vita_vars.txt", col_names = "hash")
prjna385949 <- read_tsv("outputs/vita_rf/PRJNA385949_vita_vars.txt", col_names = "hash")
prjna400072 <- read_tsv("outputs/vita_rf/PRJNA400072_vita_vars.txt", col_names = "hash")

upset_df_hash <- UpSetR::fromList(list(PRJEB2054 = prjeb2054$hash, 
                 PRJNA237362 = prjna237362$hash,
                 iHMP = ihmp$hash,
                 SRP057027 = srp057027$hash,
                 PRJNA400072 = prjna400072$hash, 
                 PRJNA385949 = prjna385949$hash))

upset_df_hash$hash <- unique(c(prjeb2054$hash, prjna237362$hash, ihmp$hash,
                               srp057027$hash, prjna400072$hash, prjna385949$hash))

# read in variable importance and add to upset dataframe

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
  arrange(desc(total_imp))
varimp_cum$hash <- as.numeric(as.character(varimp_cum$hash))

# check 1:1 match
table(upset_df_hash$hash %in% varimp_cum$hash)

upset_df_hash <- left_join(upset_df_hash, varimp_cum)

# try upset ---------------------------------------------------------------
studies = c("iHMP", "PRJEB2054", "PRJNA237362", "PRJNA385949", "PRJNA400072", "SRP057027")
ComplexUpset::upset(upset_df_hash, studies, sort_intersections_by='degree',
                    dot_size = 2, set_sizes=FALSE, name=NULL, width_ratio = 0.1,
                    base_annotations=list('k-mers'=intersection_size(counts=FALSE)),            
                    annotations = list('cumulative importance'=list(
                      aes=aes(x=intersection, total_imp/6),
                      geom=list(geom_bar(stat='identity')))))


# plot small --------------------------------------------------------------
# filter df to hashes that are in at least 5 of 6 studies
upset_df_hash_small <- upset_df_hash[rowSums(upset_df_hash[ , 1:6]) >= 5, ]


ComplexUpset::upset(upset_df_hash_small, studies, sort_intersections_by='degree',
                    dot_size = 2, set_sizes=FALSE, name=NULL, width_ratio = 0.1,
                    base_annotations=list('k-mers'=intersection_size(counts=FALSE)),            
                    annotations = list('cumulative importance'=list(
                      aes=aes(x=intersection, total_imp/6),
                      geom=list(geom_bar(stat='identity')))))
sum(upset_df_hash_small$total_imp) / 6
