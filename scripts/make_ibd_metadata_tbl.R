setwd("~/github/2020-ibd")

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)

# PRJEB2054 has unique library names, even when two libraries correspond to one 
# sample. Change the "library_name" to the "sample_alias", which will designate
# which libraries need to be combined to be the same sample. 
# do this first so that changing library name doesn't mess with the other samples
prj2054 <- read_tsv("inputs/ENA/PRJEB2054.txt") %>%
  filter(library_strategy != "AMPLICON") %>%
  filter(library_strategy != "RNA-Seq") %>%
  filter(read_count >= 1000000) %>%
  filter(library_layout  == "PAIRED")
prj2054$library_name <- prj2054$sample_alias


files <- list.files("inputs/ENA", full.names = T)
files <- files[!files %in%  "inputs/ENA/PRJEB2054.txt"]
prj <- files %>%
  map(read_tsv) %>%
  reduce(rbind) %>%
  filter(library_strategy != "AMPLICON") %>%
  filter(library_strategy != "RNA-Seq") %>%
  filter(read_count >= 1000000) %>%
  filter(library_layout  == "PAIRED")
prj <- rbind(prj, prj2054)

# sample_alias has disease status for PRJEB2054
# PRJNA385949 time series (use as validation)
# pediatric crohn's data set: https://raw.githubusercontent.com/louiejtaylor/sbx_lewis2015/master/metadata/metadata_with_SRR.csv
# PRJNA400072 matches library name

# fix sample alias for PRJNA237362
prj$sample_alias <- gsub("PRJNA237362\\.", "", prj$sample_alias)
# join with study-specific metadata
# PRJEB2054 -- sample_alias
# PRJNA400072 -- library_name
# PRJNA237362 -- sample_alias
# SRP057027 -- run_accession, library_name
# PRJNA385949 -- library_name

prjeb2054 <- read_csv("inputs/metadata/PRJEB2054_metadata.csv")
prj <- left_join(prj, prjeb2054, by = "sample_alias")
prjna237362 <- read_csv("inputs/metadata/PRJNA237362_metadata.csv")
prj <- left_join(prj, prjna237362, by = c("library_name", "sample_alias"))
# prj %>% filter(study_accession == "PRJNA237362") %>% select(library_name, sample_alias, diagnosis, individual)
prjna400072 <- read_csv("inputs/metadata/PRJNA400072_metadata.csv")
prj <- left_join(prj, prjna400072, by = "library_name") 
srp057027 <- read_csv("inputs/metadata/SRP057027_metadata.csv")
prj <- left_join(prj, srp057027, by = c("run_accession", "library_name"))
prjna385949 <- read_csv("inputs/metadata/PRJNA385949_metadata.csv")
prj <- left_join(prj, prjna385949, by = "library_name")

# coalesce colnames
prj <- prj %>%
  mutate(diagnosis = coalesce(diagnosis.x.x, diagnosis.y.y, diagnosis.x, diagnosis.y, diagnosis)) %>%
  select(-diagnosis.x.x, -diagnosis.y.y, -diagnosis.x, -diagnosis.y) %>%
  mutate(age = coalesce(age.x, age.y)) %>%
  select(-age.x, -age.y) %>%
  mutate(sample = coalesce(sample.x, sample.y)) %>%
  select(-sample.x, -sample.y) %>%
  mutate(steroids = coalesce(steroids.x, steroids.y)) %>%
  select(-steroids.x, -steroids.y) %>%
  mutate(PCDAI = coalesce(as.character(PCDAI.x), as.character(PCDAI.y))) %>%
  select(-PCDAI.x, -PCDAI.y) %>%
  mutate(antibiotic = coalesce(antibiotic.x, antibiotic.y)) %>%
  select(-antibiotic.x, -antibiotic.y) %>% 
  mutate(subject = coalesce(as.character(subject), as.character(patient_id), individual)) %>%
  select(-patient_id, -individual) %>%
  mutate(fecal_calprotectin = coalesce(fecal_calprotectin.x, as.character(fecal_calprotectin.y))) %>%
  select(-fecal_calprotectin.x, -fecal_calprotectin.y) %>%
  filter(diagnosis != "IC") # also removes NAs

wrk <- select(prj, study_accession, run_accession, library_name, read_count, 
              fastq_ftp, sample_alias, diagnosis, subject)
wrk <- separate(data = wrk, 
                col = fastq_ftp,
                into = c("fastq_ftp_1", "fastq_ftp_2", "fastq_ftp_3"), 
                sep = ";", remove = F, extra = "merge", fill = "right")
# some fastq files have 3 download links, where the first one contains what is
# probably a merged or interleaved file. Remove these from fastq download links. 
# fastq 1 download:

fastq_1 <- vector()
for(i in 1:length(wrk$fastq_ftp)){
  if(grepl(pattern = "_1", x = wrk$fastq_ftp_1[i])){
    fastq_1[i] <- wrk$fastq_ftp_1[i] 
  } else {
    fastq_1[i] <- wrk$fastq_ftp_2[i]
  }
}

fastq_2 <- vector()
for(i in 1:length(wrk$fastq_ftp)){
  if(grepl(pattern = "_1", x = wrk$fastq_ftp_1[i])){
    fastq_2[i] <- wrk$fastq_ftp_2[i] 
  } else {
    fastq_2[i] <- wrk$fastq_ftp_3[i]
  }
}

wrk$fastq_ftp_1 <- fastq_1
wrk$fastq_ftp_2 <- fastq_2

wrk <- select(wrk, study_accession, run_accession, library_name, read_count, 
              fastq_ftp_1, fastq_ftp_2, sample_alias, diagnosis, subject)
wrk$subject <- ifelse(is.na(wrk$subject), wrk$library_name, wrk$subject)

# make sure that there are 500 subjects
all.equal(wrk %>% select(subject) %>% distinct() %>% nrow(), 500)


# subsample time series ---------------------------------------------------

# subsample so that there is only one observation for each individual
# e.g., eliminate time series
# sort by library name, as lib name for time series samples is "alphabetical", 
# with the first samples alphabetically are those that were taken first. 
# Do this because some patients were treated, and we want to exclude treatment
# samples. 
wrk <- wrk[order(wrk$library_name, decreasing = F), ]
wrk <- wrk[!duplicated(wrk$library_name), ]
# check that we still have 500 subjects
all.equal(wrk %>% select(subject) %>% distinct() %>% nrow(), 500)

# combine with ihmp data --------------------------------------------------

hmp <- read_tsv("inputs/hmp2_mgx_metadata.tsv") %>%
  mutate(study_accession = "iHMP") %>%
  mutate(library_name = External.ID) %>%
  mutate(run_accession = NA) %>%
  mutate(sample_alias = NA) %>%
  mutate(fastq_ftp_1 = NA) %>%
  mutate(fastq_ftp_2 = NA) %>%
  mutate(read_count = reads_raw) %>%
  mutate(subject = Participant.ID) %>%
  arrange(subject, week_num) %>%
  select(study_accession, run_accession, library_name, read_count, fastq_ftp_1, fastq_ftp_2, sample_alias, diagnosis, subject)

hmp <- hmp[!duplicated(hmp$library_name), ]
all.equal(hmp %>% select(subject) %>% distinct() %>% nrow(), 106)


all.equal(colnames(hmp), colnames(wrk))
wrk <- rbind(wrk, hmp)

write_tsv(x = wrk, path = "inputs/working_metadata.tsv")
