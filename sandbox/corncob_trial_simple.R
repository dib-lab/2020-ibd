setwd("~/github/ibd")
library(readr)
library(phyloseq)
library(dplyr)
library(corncob)
library(ggplot2)
library(feather)


# read in metadata --------------------------------------------------------

info_libs <- read_tsv("inputs/working_metadata.tsv") %>%
  select(study_accession, library_name, subject, diagnosis) %>%
  distinct()

info_hmp <- read_tsv("inputs/hmp2_mgx_metadata.tsv") %>%
  mutate(study_accession = "hmp") %>%
  rename(library_name = External.ID, subject = Participant.ID) %>%
  select(study_accession, library_name, subject, diagnosis)

all_info <- rbind(info_libs, info_hmp)

# remove library_name from the info data frame
info <- sample_data(select(all_info, -library_name))
# Pull out library_name to use as sample names
sample_names(info) <- pull(all_info, library_name)

# read in count data ------------------------------------------------------

count_info <- read_csv("outputs/hashtables/all_unnormalized_abund_hashes_filtbyexpr.csv")
counts <- otu_table(select(count_info, -X1), taxa_are_rows = TRUE)
# Pull out the gene ID to use as "taxa" names
taxa_names(counts) <- pull(count_info, X1)

# create phyloseq object --------------------------------------------------

ibd <- phyloseq(counts, info)


# run corncob -------------------------------------------------------------

corncob <- system.time(differentialTest(formula = ~ study_accession + subject + diagnosis,
                            formula_null = ~ study_accession + subject,
                            phi.formula = ~ 1,
                            phi.formula_null = ~ 1,
                            data = ibd, 
                            test = "Wald", boot = FALSE))
# Error in differentialTest(formula = ~study_accession + subject + diagnosis,  : 
#                             All models failed to converge! 
#                             
# If you are seeing this, it is likely that your model is overspecified. This occurs when your sample size is not large enough to estimate all the parameters of your model. This is most commonly due to categorical variables that include many categories. To confirm this, try running a model for a single taxon with bbdml.
# Timing stopped at: 1.667e+05 1452 1.983e+05

corncob <- differentialTest(formula = ~ study_accession + diagnosis,
                            formula_null = ~ study_accession,
                            phi.formula = ~ 1,
                            phi.formula_null = ~ 1,
                            data = ibd, 
                            test = "Wald", boot = FALSE)
#     user   system  elapsed 
# 2466.805  192.198 2679.383 