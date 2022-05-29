library(dplyr)
library(readr)
library(tidyr)

ihmp_samples <- read_tsv("~/Downloads/iHMP_106_samples.txt") %>%
  select(library_name, subject, original_read_count = read_count, diagnosis)
bbduk <- read_delim("~/Downloads/file_sizes_bbduk.txt", delim = " ", 
                    col_names = c("permissions", "rm", "owner", "owner2", 
                                  "size", "month", "day", "year", "file")) %>%
  mutate(library_name = gsub(".fq.gz", "", file)) %>%
  separate(library_name, into = c("library_name", "pair"), sep = "_") %>%
  separate(pair, into = c("pair", "type"), sep = "\\.") %>%
  select(library_name, pair, size, type)

trim <- read_delim("~/Downloads/file_sizes_trim.txt", delim = " ", 
                   col_names = c("permissions", "rm", "owner", "owner2", 
                                 "size", "month", "day", "year", "file")) %>%
  mutate(library_name = gsub(".fq.gz", "", file)) %>%
  separate(library_name, into = c("library_name", "pair"), sep = "_") %>%
  separate(pair, into = c("pair", "type"), sep = "\\.") %>%
  select(library_name, pair, size, type)

raw <- read_delim("~/Downloads/file_sizes_raw.txt", delim = " ", 
                   col_names = c("permissions", "rm", "owner", "owner2", 
                                 "size", "month", "day", "year", "file")) %>%
  mutate(library_name = gsub(".fq.gz", "", file)) %>%
  separate(library_name, into = c("library_name", "pair"), sep = "_") %>%
  separate(pair, into = c("pair", "type"), sep = "\\.") %>%
  select(library_name, pair, size, type) %>%
  mutate(pair = paste0("R", pair)) %>%
  mutate(type = gsub("fastq", "raw", type))

sizes <- rbind(bbduk, trim, raw)
sizes <- ihmp_samples %>%
  left_join(sizes, by = "library_name") %>%
  filter(pair %in% c("R1", "R2"))

tmp <- pivot_wider(sizes, id_cols = library_name:pair, names_from = type, values_from = size) %>%
  select(library_name, subject, original_read_count, diagnosis, pair, raw, trim, nohost, human) 
View(tmp)

write_tsv(tmp, "~/Downloads/ihmp_contamination.tsv")
