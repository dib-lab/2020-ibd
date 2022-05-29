library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct

generate_counts <- function(path){
  # read in an format counts
  counts <- read_tsv(path)
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$protein
  counts <- counts[ , -1]
}

generate_presence_absence <- function(path, threshold = 0){
  counts <- generate_counts(path)
  # replace any count > 0 with 1
  counts_pa <- ifelse(counts > threshold, 1, 0)
}

counts <- generate_presence_absence(path = "sandbox/tximport/SRS294916_20.fna_counts_raw.tsv", threshold = 0)

annot_df <- data.frame(library_name = colnames(counts))
annot_df <- left_join(annot_df, info, by = "library_name") %>%
  select(library_name, diagnosis)
annot_df <- as.data.frame(annot_df)
rownames(annot_df) <- annot_df$library_name
annot_df <- annot_df["diagnosis"]
pheatmap(as.matrix(counts), cluster_rows = T, cluster_cols = T, show_rownames = F,
         show_colnames = F, annotation_col = annot_df)



counts <- generate_presence_absence(path = "sandbox/tximport/SRS1735645_19.fna_counts_raw.tsv", threshold = 0)

annot_df <- data.frame(library_name = colnames(counts))
annot_df <- left_join(annot_df, info, by = "library_name") %>%
  select(library_name, diagnosis)
annot_df <- as.data.frame(annot_df)
rownames(annot_df) <- annot_df$library_name
annot_df <- annot_df["diagnosis"]
pheatmap(as.matrix(counts), cluster_rows = T, cluster_cols = T, show_rownames = F,
         show_colnames = F, annotation_col = annot_df)
