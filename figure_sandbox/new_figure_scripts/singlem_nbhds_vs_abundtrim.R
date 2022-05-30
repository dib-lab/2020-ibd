setwd("~/github/2020-ibd")
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

singlem_sgc <- read_tsv("outputs/sgc_genome_queries_singlem/combined.tsv")
singlem_abundtrim <- read_tsv("outputs/abundtrim_singlem/combined.tsv") %>%
  mutate(sample = source)

# sample, sequence
singlem_sgc <- singlem_sgc %>%
  select(sample, sequence) %>%
  distinct() %>%
  mutate(source = "sgc")

singlem_abundtrim <- singlem_abundtrim %>%
  select(sample, sequence) %>%
  distinct() %>%
  mutate(source = "abundtrim")

singlem <- rbind(singlem_sgc, singlem_abundtrim)

singlem_tally <- singlem %>%
  group_by(sample, source) %>%
  tally() %>%
  pivot_wider(id_cols = sample, names_from = source, values_from = n) %>%
  mutate(fraction = sgc/abundtrim)

ggplot(singlem_tally, aes(x = abundtrim, y = sgc)) + 
  geom_point(alpha = .1) +
  ylim(0, 4500) +
  theme_minimal() + 
  labs(x = "error-trimmed metagenome", y = "microbial core")

ggplot(singlem_tally, aes(x = fraction)) +
  geom_density() +
  theme_minimal() + 
  xlim(0, 1) +
  ggtitle("fraction of marker genes in microbial core compared to whole metagenome")

