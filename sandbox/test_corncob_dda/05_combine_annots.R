library(dplyr)
library(readr)
library(tidyr)
setwd("~/github/2020-ibd/sandbox")

eggnog <- read_tsv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/GCA_900036035.1_eggnog/GCA_900036035.1.emapper.annotations", skip = 3) %>%
  rename(query_name = `#query_name`)
cdbg_annot <- read_csv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10_multifasta/multifasta.cdbg_annot.csv") %>%
  separate(record_name, into = c("record_id", "tmp"), remove = F, sep = " ") %>%
  select(-tmp)

sig_ccs <- read_tsv("test_corncob_dda/corncob_results_sig_ccs.tsv") %>%
  rename(dom_id = aa_seq)
cdbg_to_pieces <- read_csv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/cdbg_to_pieces.csv")
dom_info <- read_tsv("test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_dom_info.tsv")

# filter to sig_ccs
sig_ccs <- left_join(sig_ccs, cdbg_to_pieces, by = c("dom_id" = "dominator"))
sig_ccs <- left_join(sig_ccs, cdbg_annot, by = c("cdbg_node" = "cdbg_id"))
sig_ccs <- left_join(sig_ccs, eggnog, by = c("record_id" = "query_name"))

write_tsv(sig_ccs, "test_corncob_dda/corncob_results_sig_ccs_annotated.tsv")

# rm cdbg_node
sig_ccs <- sig_ccs %>%
  select(-cdbg_node) %>%
  distinct() 

sig_ccs <- left_join(sig_ccs, dom_info, by = "dom_id")

write_tsv(sig_ccs, "test_corncob_dda/corncob_results_sig_ccs_annotated_no_cdbg.tsv")

