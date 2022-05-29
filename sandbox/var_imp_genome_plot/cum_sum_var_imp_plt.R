# the purpose of this script is to plot the proportion of each of the 129
# predictive genomes from vita variable selection based on when they appear
# in importance rankings. i.e. only 1 genome contains the most
# predictive hash

setwd("~/github/2020-ibd/")
library(ranger)
library(dplyr)
library(ggplot2)
library(readr)

# read in ranger and get variable importance
optimal_ranger <- readRDS("outputs/optimal_rf/optimal_ranger.RDS")
imp <- optimal_ranger$variable.importance
imp <- t(data.frame(as.list(imp)))
imp <- as.data.frame(imp)
imp$hash <- rownames(imp)
imp$hash <- gsub("^X", "", imp$hash)
imp <- imp %>%
  filter(hash != "diagnosis")
imp$hash <- as.numeric(imp$hash)
imp <- imp[order(imp$importance, decreasing = T), ]
imp$rank <- 1:nrow(imp)
colnames(imp) <- c("importance", "hash", "rank")
  
hash_to_genome_mapping <- read_csv("sandbox/var_imp_genome_plot/vita_hash_to_129_genome_mapping.csv")
hash_to_genome_mapping$genome <- gsub("gather_genome_129_sigs\\/", "", hash_to_genome_mapping$genome)
hash_to_genome_mapping$genome <- gsub("\\.fna\\.sig", "", hash_to_genome_mapping$genome)

all <- full_join(imp, hash_to_genome_mapping, by = c("hash" = "index"))

# keep only the first observation of each genome. 
# note the df needs to be ordered by rank for this to behave as expected.
all_rank <- all[!duplicated(all$genome),]

all_rank_summarized <- all_rank %>%
  group_by(rank) %>%
  tally()

# calculate the cumulative sum at each rank
cumulative_sum_all <- vector()
for(row in 1:nrow(all_rank_summarized)){
  if(row == 1){
    cumulative_sum <- all_rank_summarized$n[row]
  } else {
    cumulative_sum <- cumulative_sum + all_rank_summarized$n[row]
  }
  cumulative_sum_all[row] <- cumulative_sum
}

all_rank_summarized$cumulative_sum <- cumulative_sum_all
  
ggplot(all_rank_summarized, aes(x = rank, y = cumulative_sum)) +
  geom_point() +
  theme_minimal()
