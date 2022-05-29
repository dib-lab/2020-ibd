# figure out which hashes don't make it into the megahit assemblies and 
# what the identity of those hashes are

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

#  scp -P 2022 -i ~/.ssh/key_name.pem tereiter@farm.cse.ucdavis.edu:/home/tereiter/github/2020-ibd/all-query-results.csv .
query_res <- read_csv("~/Downloads/all-query-results.csv") %>%
  mutate(filename = gsub("/home/tereiter/github/2020-ibd/outputs/gather_matches_loso_prokka/", "", filename)) %>%
  mutate(filename = gsub(".ffn", "", filename)) %>%
  mutate(sample = gsub("_k31_r1", "", catlas_base)) %>%
  separate(col = record_name, into = c("orf_name", "annotation"), sep = " ", remove = F, extra = "merge")

# join with GTDBTK annotation information ------------------
lca <- read_csv("inputs/at_least_5_studies_vita_vars_gather_all_lca.csv")
query_res <- left_join(query_res, lca, by = c("filename" = "name_no_gz"))

# join with hash importance information --------------------
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

query_res <- left_join(query_res, varimp_cum, by = c("hashval" = "hash"))

# read in hashes that didn't make it into the megahit assemblies
# scp -P 2022 -i ~/.ssh/key_name.pem tereiter@farm.cse.ucdavis.edu:~/github/2020-ibd/sandbox/test_megahit_diginorm_nocat/megahit_sigs/at_least_5_studies_vita_vars_minus_all_megahit_sigs.csv .
unassembled_hashes <- read_csv("~/Downloads/at_least_5_studies_vita_vars_minus_all_megahit_sigs.csv")


# compare -----------------------------------------------------------------

unassembled_hashes <- query_res %>%
  filter(hashval %in% unassembled_hashes$minhash) %>%
  select(hashval, annotation) %>%
  distinct()

length(unique(unassembled_hashes$hashval))
View(unassembled_hashes)

unassembled_hashes2 <- query_res %>%
  filter(hashval %in% unassembled_hashes$hashval)
View(unassembled_hashes2)
assembled_hashes <- query_res %>%
  filter(!hashval %in% unassembled_hashes$hashval) %>%
  select(hashval, annotation) %>%
  distinct()

View(assembled_hashes)
