library(dplyr)
library(readr)
library(tidyr)

eggnog <- read_tsv(snakemake@input[['eggnog']], skip = 4, show_col_types = F) %>%
#eggnog <- read_tsv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/GCA_900036035.1_eggnog/GCA_900036035.1.emapper.annotations", skip = 3) %>%
  rename(query_name = `#query`) %>%
  mutate(query_name = gsub("_1$", "", query_name)) # remove trailling one added by eggnog or clustering

cdbg_annot <- read_csv(snakemake@input[['cdbg_annot']], show_col_types= F) %>%
#cdbg_annot <- read_csv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10_multifasta/multifasta.cdbg_annot.csv") %>%
  separate(record_name, into = c("record_id", "tmp"), remove = F, sep = " ") %>%
  select(-tmp)

sig_ccs <- read_tsv(snakemake@input[['sig_ccs']]) %>%
# sig_ccs <- read_tsv("test_corncob_dda/corncob_results_sig_ccs.tsv") %>%
  rename(dom_id = aa_seq)

cdbg_to_pieces <- read_csv(snakemake@input[['cdbg_to_pieces']])
#cdbg_to_pieces <- read_csv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/cdbg_to_pieces.csv")
dom_info <- read_tsv(snakemake@input[['dom_info']])
#dom_info <- read_tsv("test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_dom_info.tsv")

# filter to sig_ccs
sig_ccs <- left_join(sig_ccs, cdbg_to_pieces, by = c("dom_id" = "dominator"))
sig_ccs <- left_join(sig_ccs, cdbg_annot, by = c("cdbg_node" = "cdbg_id"))
sig_ccs <- left_join(sig_ccs, eggnog, by = c("record_id" = "query_name"))

#write_tsv(sig_ccs, "test_corncob_dda/corncob_results_sig_ccs_annotated.tsv")

# rm cdbg_node
sig_ccs <- sig_ccs %>%
  select(-cdbg_node) %>%
  distinct() 

sig_ccs <- left_join(sig_ccs, dom_info, by = "dom_id")

write_tsv(sig_ccs, snakemake@output[['sig_ccs_annot']])
# write_tsv(sig_ccs, "test_corncob_dda/corncob_results_sig_ccs_annotated_no_cdbg.tsv")
