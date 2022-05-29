setwd("~/github/ibd")
library(vegan)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ggrepel)
library(stringr)
library(plotly)
library(viridis)


# read in comp and format  ----------------------------------------------
comp <- read_csv("sandbox/02_filt_comp_all/comp_all_jaccard.csv") # read in mat
colnames(comp) <- gsub("_filt", "", colnames(comp))
rownames(comp) <-colnames(comp) 


# read in and formate metadata -------------------------------------------

info <- read_tsv("inputs/working_metadata.tsv") %>%
  filter(library_name %in% colnames(comp)) %>%
  group_by(library_name, study_accession, diagnosis) %>%
  summarise(read_count = sum(read_count)) %>%
  as.data.frame()

info_hmp <- read_tsv("inputs/hmp2_mgx_metadata.tsv") %>%
  mutate(study_accession = "hmp") %>%
  select(External.ID, study_accession, diagnosis, reads_raw) %>%
  distinct() %>%
  filter(External.ID %in% colnames(comp))
colnames(info_hmp) <- c("library_name", "study_accession", "diagnosis", "read_count")

info <- rbind(info, info_hmp)


# dist and permanova ------------------------------------------------------

dist <- dist(comp)                                     # compute distances
# sort info by colnames
info <- info[match(colnames(comp), info$library_name), ]
info$diagnosis <- as.factor(info$diagnosis)
info$study_accession <- as.factor(info$study_accession)

perm <- adonis(dist ~ diagnosis*study_accession*read_count, 
               data = info, 
               permutations = 1000)

perm
# Call:
#   adonis(formula = dist ~ diagnosis * study_accession * read_count,      data = info, permutations = 1000) 
# 
# Permutation: free
# Number of permutations: 1000
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model
# diagnosis                               2     319.7 159.832  48.260
# study_accession                         5    1124.9 224.972  67.929
# read_count                              1      62.0  62.029  18.729
# diagnosis:study_accession               7      93.1  13.301   4.016
# diagnosis:read_count                    2      20.1  10.070   3.041
# study_accession:read_count              5      91.7  18.330   5.535
# diagnosis:study_accession:read_count    7      24.5   3.496   1.055
# Residuals                            2262    7491.4   3.312        
# Total                                2291    9227.3                
# R2   Pr(>F)    
# diagnosis                            0.03464 0.000999 ***
#   study_accession                      0.12191 0.000999 ***
#   read_count                           0.00672 0.000999 ***
#   diagnosis:study_accession            0.01009 0.000999 ***
#   diagnosis:read_count                 0.00218 0.007992 ** 
#   study_accession:read_count           0.00993 0.000999 ***
#   diagnosis:study_accession:read_count 0.00265 0.329670    
# Residuals                            0.81187             
# Total                                1.00000             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




