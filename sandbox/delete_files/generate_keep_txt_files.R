library(readr)
library(dplyr)

setwd("github/2020-ibd")

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name) %>%
  distinct()

# trim
trim_R1 <- info %>%
  mutate(library_name = paste0(library_name, "_R1.trim.fq.gz"))
trim_R2 <- info %>%
  mutate(library_name = paste0(library_name, "_R2.trim.fq.gz"))
trim_o1 <- info %>%
  mutate(library_name = paste0(library_name, "_o1.trim.fq.gz"))
trim_o2 <- info %>%
  mutate(library_name = paste0(library_name, "_o2.trim.fq.gz"))
write_tsv(rbind(trim_R1, trim_R2, trim_o1, trim_o2), 
          "sandbox/delete_files/keep_trim.txt", col_names = F)

# cat
cat_R1 <- info %>%
  mutate(library_name = paste0(library_name, "_1.fastq.gz"))
cat_R2 <- info %>%
  mutate(library_name = paste0(library_name, "_2.fastq.gz"))
write_tsv(rbind(cat_R1, cat_R2), "sandbox/delete_files/keep_cat.txt", col_names = F)  

# bbduk
bbduk_R1_nohost <- info %>%
  mutate(library_name = paste0(library_name, "_R1.nohost.fq.gz"))
bbduk_R2_nohost <- info %>%
  mutate(library_name = paste0(library_name, "_R2.nohost.fq.gz"))
bbduk_R1_human <- info %>%
  mutate(library_name = paste0(library_name, "_R1.human.fq.gz"))
bbduk_R2_human <- info %>%
  mutate(library_name = paste0(library_name, "_R2.human.fq.gz"))
write_tsv(rbind(bbduk_R1_nohost, bbduk_R2_nohost, bbduk_R1_human, bbduk_R2_human), 
          "sandbox/delete_files/keep_bbduk.txt", col_names = F)

# abundtrim
info %>%
  mutate(library_name = paste0(library_name, ".abundtrim.fq.gz")) %>%
  write_tsv("sandbox/delete_files/keep_abundtrim.txt", col_names = F)

# cut
cut_R1 <- info %>%
  mutate(library_name = paste0(library_name, "_R1.cut.fq.gz"))
cut_R2 <- info %>%
  mutate(library_name = paste0(library_name, "_R2.cut.fq.gz"))
write_tsv(rbind(cut_R1, cut_R2), "sandbox/delete_files/keep_cut.txt", col_names = F)  
