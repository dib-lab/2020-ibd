library(readr)
library(dplyr)
library(ggplot2)
library(vegan)
library(pheatmap)
library(microbiome)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct()

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

annot_col_df <- data.frame(library_name = colnames(rt)) %>%
  left_join(info) %>%
  select(library_name, diagnosis) %>%
  column_to_rownames("library_name")

rt <- generate_counts("sandbox/test_megahit_diginorm_nocat/tximport/SRS294916_20.fna_counts_raw.tsv")
rt <- transform(rt, transform = "clr")
pheatmap(rt, cluster_rows = T, cluster_cols = T, annotation_col = annot_col_df, 
         show_rownames = F, show_colnames = F)

cbolt <- generate_counts("sandbox/test_megahit_diginorm_nocat/tximport/GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna_counts_raw.tsv")
cbolt <- transform(cbolt, transform = "clr")
min(cbolt)
annot_col_df <- data.frame(library_name = colnames(cbolt)) %>%
  left_join(info) %>%
  select(library_name, diagnosis) %>%
  column_to_rownames("library_name")
pheatmap(cbolt, cluster_rows = T, cluster_cols = T, annotation_col = annot_col_df, 
         show_rownames = F, show_colnames = F)

generate_presence_absence <- function(path, threshold = 0){
  counts <- generate_counts(path)
  # replace any count > 0 with 1
  counts_pa <- ifelse(counts > threshold, 1, 0)
}

cbolt_pa <- generate_presence_absence("sandbox/test_megahit_diginorm_nocat/tximport/GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna_counts_raw.tsv",
                                      threshold = 0)
annot_col_df <- data.frame(library_name = colnames(cbolt_pa)) %>%
  left_join(info) %>%
  select(library_name, diagnosis) %>%
  column_to_rownames("library_name")
pheatmap(cbolt_pa, cluster_rows = T, cluster_cols = T, annotation_col = annot_col_df, 
         show_rownames = F, show_colnames = F)

generate_presence_absence <- function(path, threshold = 0){
  counts <- generate_counts(path)
  # replace any count > 0 with 1
  counts_pa <- ifelse(counts > threshold, 1, 0)
}

cbolt_pa <- generate_presence_absence("sandbox/test_megahit_diginorm_nocat/tximport/GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna_counts_raw.tsv",
                                      threshold = 3)
annot_col_df <- data.frame(library_name = colnames(cbolt_pa)) %>%
  left_join(info) %>%
  select(library_name, diagnosis) %>%
  column_to_rownames("library_name")
pheatmap(cbolt_pa, cluster_rows = T, cluster_cols = T, annotation_col = annot_col_df, 
         show_rownames = F, show_colnames = F)



rgnv_pa <- generate_presence_absence("sandbox/test_megahit_diginorm_nocat/tximport/GCF_900036035.1_RGNV35913_genomic.fna_counts_raw.tsv",
                                      threshold = 0)
annot_col_df <- data.frame(library_name = colnames(rgnv_pa)) %>%
  left_join(info) %>%
  select(library_name, diagnosis) %>%
  column_to_rownames("library_name")
pheatmap(rgnv_pa, cluster_rows = T, cluster_cols = T, annotation_col = annot_col_df, 
         show_rownames = F, show_colnames = F)
