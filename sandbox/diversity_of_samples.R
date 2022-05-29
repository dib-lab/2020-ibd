# tmp test diversity in non-IBD, UC, and CD. 
library(readr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
abundtrim <- read_tsv("~/Downloads/diversity/mqc_fastp_filtered_reads_plot_1.txt",
                      col_names = c("sample", "passed", "low_qual", "too_many_n", "too_short"), 
                      skip = 1) %>%
  mutate(sample = gsub("\\.abundtrim", "", sample)) %>%
  mutate(total_reads = rowSums(.[2:5])) %>%
  select(sample, total_reads)
filt_sigs <- read_csv("~/Downloads/diversity/sig_describe_filt_named_sig.csv") %>%
  mutate(signature_file = gsub("_filt_named\\.sig", "", signature_file))

all <- left_join(filt_sigs, abundtrim, by = c('signature_file' = 'sample'))

all$diversity <- all$n_hashes / all$total_reads
all <- all %>%
  select(signature_file, n_hashes, total_reads, diversity)

info <- read_tsv("~/github/2020-ibd/inputs/working_metadata.tsv") %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct

all <- left_join(all, info, by = c("signature_file" = "library_name"))

ggplot(all, aes(x = diagnosis, y = diversity)) +
  geom_boxplot() +
  theme_minimal() 

fm1  = aov(diversity ~ diagnosis, data = all)
anova(fm1)

TukeyHSD(fm1)



ggplot(all, aes(x = diagnosis, y = diversity)) +
  geom_beeswarm(aes(color = total_reads)) +
  theme_minimal() 

ggplot(all, aes(x = diagnosis, y = n_hashes)) +
  geom_beeswarm(aes(color = total_reads)) +
  theme_minimal() 


# trying camille's method -------------------------------------------------

# 1 - (number of hashes over estimated unique k-mers)
# 7376151 greater_than_one_count_hashes.txt

all$diversity <- (all$n_hashes / 7376151)

ggplot(all, aes(x = diagnosis, y = diversity)) +
  geom_beeswarm(aes(color = total_reads)) +
  theme_minimal() 


fm1  = aov(diversity ~ diagnosis, data = all)
anova(fm1)

TukeyHSD(fm1)
