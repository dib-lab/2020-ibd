library(readr)
library(dplyr)
library(tibble)

# use cbind approach ------------------------------------------------------

# It's a bad idea to import the files this way because there is no check that
# dom_ids are the same between all samples other than that they will have the
# number of rows. However, this is much faster and will probably be the
# solution that I use to because it will scale to hundreds or thousands of
# samples with large numbers of dominating sets.
acc <- snakemake@wildcards[['acc']]
file_prefix <- paste0("outputs/sgc_pangenome_catlases/", acc, "_k31_r10_abund/")

domset_abund_files <- snakemake@input[['dom_abund']]
domnames <- read_csv(domset_abund_files[1])[ , c("dom_id", "level")]     # read in domset id
dom_info <- read_csv(domset_abund_files[1]) %>%
  select(-abund) %>%
  mutate(dom_id = as.character(dom_id))

df <- do.call(cbind, lapply(domset_abund_files,
                            function(x) read_csv(x)[ , "abund"]))
df <- cbind(domnames, df)
colnames(df) <- c("dom_id", "level", domset_abund_files)
colnames(df) <- gsub(file_prefix, "", colnames(df))
colnames(df) <- gsub(".reads.gz.dom_abund.csv", "", colnames(df))

write_tsv(df, snakemake@output[['dom_abund']])
write_tsv(dom_info, snakemake@output[['dom_info']])

# prune to nodes present in > 100 samples at level 1 ---------------------------

df <- df %>%
  filter(level == 1) %>%
  select(-level)

tmp_pa <- df
rownames(tmp_pa) <- df$dom_id
tmp_pa <- tmp_pa[ , -c(1)]
tmp_pa[tmp_pa>0] <-1
row_sums_pa <- rowSums(tmp_pa)
row_sums_pa <- as.data.frame(row_sums_pa)
rownames(row_sums_pa) <- rownames(tmp_pa)

row_sums_pa <- row_sums_pa %>%
  rownames_to_column("dom_id") %>%
  filter(row_sums_pa >= 100) 

df_pruned <- row_sums_pa %>%
  select(dom_id) %>%
  mutate(dom_id = as.numeric(dom_id)) %>%
  left_join(df, by = "dom_id") %>%
  mutate(dom_id = as.character(dom_id))

write_tsv(df_pruned, snakemake@output[['dom_abund_pruned']])
