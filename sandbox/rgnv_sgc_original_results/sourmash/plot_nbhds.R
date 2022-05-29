library(readr)
library(dplyr)
library(pheatmap)
library(tibble)
library(RColorBrewer)
setwd("~/github/2020-ibd")

comp <- read_csv("sandbox/rgnv_sgc_original_results/sourmash/rgnv_k31_scaled2000.comp.csv")
colnames(comp) <- gsub("_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz",
                       "", colnames(comp))
rownames(comp) <- colnames(comp)  

info <- read_tsv("inputs/working_metadata.tsv")
annot_col <- info %>%
  select(library_name, diagnosis) %>%
  distinct() %>%
  as.data.frame() %>%
  column_to_rownames('library_name')

annot_col$diagnosis <- factor(annot_col$diagnosis, levels = c("nonIBD", "CD", "UC"))

mycolors <- list(diagnosis = c("#29908b", "#430753", "#fce640"))
names(mycolors$diagnosis) <- levels(annot_col$diagnosis)

pheatmap(as.matrix(comp), annotation_col = annot_col, annotation_row = annot_col,
         show_rownames = F, show_colnames = F, annotation_colors = mycolors,
         breaks = seq(0, 0.4, length.out = 100))
