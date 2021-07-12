brary(readr)
library(dplyr)
library(tibble)
setwd("~/github/2020-ibd")

# use cbind approach ------------------------------------------------------

# It's a bad idea to import hte files this way because there is no check that
# dom_ids are the same between all samples other than that they will have the
# number of rows. However, this is much faster and will probably be the
# solution that I use to because it will scale to hundreds or thousands of
# samples with large numbers of dominating sets.

domset_abund_files <- Sys.glob("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10_abund/*dom_abund.csv")
domnames <- read_csv(domset_abund_files[1])[ , "dom_id"]     # read in domset id
dom_info <- read_csv(domset_abund_files[1]) %>%
  select(-abund) %>%
  mutate(dom_id = as.character(dom_id))

df <- do.call(cbind, lapply(domset_abund_files,
                            function(x) read_csv(x)[ , "abund"]))
df <- cbind(domnames, df)
colnames(df) <- c("dom_id", domset_abund_files)
colnames(df) <- gsub("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10_abund/", "", colnames(df))
colnames(df) <- gsub("_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz.dom_abund.csv", "", colnames(df))

write_tsv(df, "sandbox/test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_abund.tsv")
write_tsv(dom_info, "sandbox/test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_dom_info.tsv")

# prune to nodes present in > 100 samples ---------------------------------

tmp_pa <- df
rownames(tmp_pa) <- df$dom_id
tmp_pa <- tmp_pa[ , -c(1:3)]
tmp_pa[tmp_pa>0] <-1
row_sums_pa <- rowSums(tmp_pa)
row_sums_pa <- as.data.frame(row_sums_pa)
rownames(row_sums_pa) <- rownames(tmp_pa)
# filter to those that occur in 100 samples or more -- 25,910 pieces
row_sums_pa <- row_sums_pa %>%
  rownames_to_column("dom_id") %>%
  filter(row_sums_pa >= 100) 

df_pruned <- row_sums_pa %>%
  select(dom_id) %>%
  mutate(dom_id = as.numeric(dom_id)) %>%
  left_join(df, by = "dom_id") %>%
  mutate(dom_id = as.character(dom_id))

write_tsv(df_pruned, "sandbox/test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_abund_pruned.tsv")

