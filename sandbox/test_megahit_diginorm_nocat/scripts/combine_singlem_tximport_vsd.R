library(purrr)
library(readr)
library(dplyr)
library(DESeq2)

# import singlem
singlem_counts <- read_tsv("../../outputs/sgc_genome_queries_singlem/singlem_counts.tsv")
singlem_counts[is.na(singlem_counts)] <- 0                      # replace 0 with NA
singlem_sample_names <- singlem_counts$sample               # set sample as rownames
singlem_counts <- singlem_counts[ , -1]                         # remove the samples column
singlem_counts <- as.data.frame(t(singlem_counts))
colnames(singlem_counts) <- singlem_sample_names
print("singlem imported")
# import tximport
files <- unlist(snakemake@input[['tximport']])
counts <- files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genome") %>%
  mutate(genome = gsub("_counts_raw.tsv", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/tximport/", "", genome)) %>%
  mutate(protein = paste0(genome, ":", protein)) %>%
  select(-genome) %>%
  select(protein, colnames(singlem_counts)) # some of the samples are dropped by singlem, so filt to those
print("counts imported")
# make a dataframe, put protein name as gene name, and then remove that column
counts <- as.data.frame(counts)
rownames(counts) <- counts$protein
counts <- counts[ , -1]

# join by columnname
stopifnot(all.equal(colnames(counts), colnames(singlem_counts)))
counts <- rbind(counts, singlem_counts)
print("count combined")
# round counts to nearest integer
counts <- apply(counts, 1:2, round)
print("counts rounded")
# variance stabilize transform
vsd <- vst(as.matrix(counts))
print("counts vsd")

vsd <- as.data.frame(vsd)
vsd$protein <- rownames(vsd)
print("writing counts")
write_tsv(vsd, snakemake@output[['vsd']])
