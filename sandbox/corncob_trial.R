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

#count_info<- read_feather("all_unnormalized XXX .feather")
count_info <- read_csv("outputs/hashtables/all_unnormalized_abund_hashes_filtbyexpr.csv")
#counts <- otu_table(select(count_info, -hash), taxa_are_rows = TRUE)
counts <- otu_table(select(count_info, -X1), taxa_are_rows = TRUE)
# Pull out the gene ID to use as "taxa" names
#taxa_names(counts) <- pull(count_info, hash)
taxa_names(counts) <- pull(count_info, X1)

# create phyloseq object --------------------------------------------------

ibd <- phyloseq(counts, info)


# run corncob -------------------------------------------------------------

# For purposes of this tutorial, we will test for differential abundance across 
# BMI, controlling for the effect of Sex on the mean relative abundance, and 
# controlling for the effect of both Sex and BMI on the dispersion of the 
# relative abundances. For a more detailed explanation of how we take models 
# such as this and translate them into corncob, see corncob_tutorial.
corncob <- differentialTest(formula = ~ study_accession + diagnosis,
                            formula_null = ~ study_accession,
                            phi.formula = ~ study_accession + diagnosis,
                            phi.formula_null = ~ study_accession,
                            data = ibd, 
                            test = "Wald", boot = FALSE)
saveRDS(corncob, "~/github/ibd/sandbox/corncob_trial.RDS")

