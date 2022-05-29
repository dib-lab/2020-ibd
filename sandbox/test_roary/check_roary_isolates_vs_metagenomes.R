library(dplyr)
library(readr)
library(ggplot2)
library(tibble)
setwd("~/github/2020-ibd")

pa <- read_csv("sandbox/test_roary/outputs/roary_with_megahit_and_isolates/gene_presence_absence.csv")

ggplot(pa, aes(x = `No. isolates`)) +
  geom_histogram() +
  scale_y_sqrt()

# filter out information columns
# Genes that were not observed in a sample are marked with NA;
# replace NAs with 0, where 0 becomes unobserved.
# Convert character columns into numeric, preserving the 0s but coercing 
# character strings into NAs. Replace NAs with 1s, where 1 becomes observed
pa_no_info <- pa %>%
  column_to_rownames("Gene") %>%
  select(-'Non-unique Gene name', -'Annotation', -'No. isolates', -'No. sequences', 
         -'Avg sequences per isolate', -'Genome Fragment',-'Order within Fragment',
         -'Accessory Fragment', -'Accessory Order with Fragment', -'QC',
         -'Min group size nuc', -'Max group size nuc', -'Avg group size nuc') %>%
  replace(is.na(.), 0) %>%
  mutate_all(., as.numeric) %>% # introduce NAs by coercion 
  replace(is.na(.), 1)
# check that that worked; all cases should be true
table(tmp == pa$`No. isolates`)
# mostly off by one, good enough for now
tmp[tmp != pa$`No. isolates`]
pa$`No. isolates`[pa$`No. isolates` != tmp]

pa_assemblies <- pa_no_info %>%
  select(-starts_with("GCF"))
pa_assemblies_0 <- rowSums(pa_assemblies)
# the isolates add 14,458 gene sequences not found in the metagenomes
table(pa_assemblies_0 == 0)

# 
pa_isolates <- pa_no_info %>%
  select(starts_with("GCF"))
pa_isolates_0 <- rowSums(pa_isolates)
# the assemblies add 16,700 gene sequences not found in the isolates
table(pa_isolates_0 == 0)
