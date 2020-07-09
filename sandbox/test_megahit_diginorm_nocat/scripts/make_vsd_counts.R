library(purrr)
library(readr)
library(dplyr)
library(DESeq2)

files <- unlist(snakemake@input[['tximport']])

counts <- files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genome") %>%
  mutate(genome = gsub("_counts_raw.tsv", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/tximport", "", genome)) %>%
  mutate(protein = paste0(genome, ":", protein)) %>%
  select(-genome)

counts <- as.data.frame(counts)
rownames(counts) <- counts$protein
counts <- counts[ , -1]
counts <- apply(counts, 1:2, round)
vsd <- vst(as.matrix(counts))

vsd$protein <- rownames(vsd)
write_tsv(as.data.frame(vsd), snakemake@output[['vsd']])
