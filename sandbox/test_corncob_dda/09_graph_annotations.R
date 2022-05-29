library(dplyr)
library(readr)
library(tidyr)
setwd("~/github/2020-ibd/sandbox")


# read in metadata --------------------------------------------------------
eggnog <- read_tsv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/megahit_and_isolates_eggnog/pan_genome_reference.emapper.annotations", skip = 3) %>%
  rename(query_name = `#query_name`) %>%
  mutate(query_name = gsub("_1", "", query_name))
cdbg_annot <- read_csv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10_multifasta_roary_isolates_and_metagenomes/multifasta.cdbg_annot.csv") %>%
  separate(record_name, into = c("record_id", "tmp"), remove = F, sep = " ") %>%
  select(-tmp)

sig_ccs <- read_tsv("test_corncob_dda/corncob_results_sig_ccs.tsv") %>%
  rename(dom_id = aa_seq)
cdbg_to_pieces <- read_csv("test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/cdbg_to_pieces.csv")
dom_info <- read_tsv("test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_dom_info.tsv")


# QUESTIONS ---------------------------------------------------------------

# how many level 1 pieces have at least one annotation? -------------------
annots_level1 <- dom_info %>%
  filter(level == 1) %>%
  left_join(cdbg_to_pieces, by = c("dom_id" = "dominator")) %>%
  left_join(cdbg_annot, by = c("cdbg_node" = "cdbg_id"))

tmp <- annots_level1 %>%
  group_by(dom_id) %>%
  arrange(record_id) %>%
  slice(1) %>%
  ungroup()

table(is.na(tmp$record_id))
# FALSE   TRUE
# 34583 250759
# 34,583 domset pieces had an annotation

length(unique(cdbg_annot$cdbg_id))
# 1562541 cDBG nodes were annotated

# how many pieces are annotated as a marker gene? ------------------------
scg <- read_csv("../inputs/kegg_scg_13059_2015_610_MOESM1_ESM.csv", skip = 2)$KO
eggnog$KEGG_ko2 <- gsub("ko:", "", eggnog$KEGG_ko)
table(eggnog$KEGG_ko2 %in% scg)
# FALSE  TRUE
# 25459   458
# 458 pan genome sequences were annotated as single copy genes (e.g. in the multifasta query used for annotation)

annots_level1 <- annots_level1 %>%
  left_join(eggnog, by = c("record_id" = "query_name"))

annots_level1_scg <- annots_level1 %>%
  filter(KEGG_ko2 %in% scg)
length(unique(annots_level1_scg$cdbg_node))
# 63083 cDBG nodes were annotated as single copy marker genes

length(unique(annots_level1_scg$dom_id))
# 1097 domset level 1 pieces contained at least 1 single copy marker gene annotation

tmp <- annots_level1_scg %>%
  group_by(dom_id, KEGG_ko2) %>%
  tally() %>%
  group_by(dom_id) %>%
  tally() %>%
  arrange(desc(n))
table(tmp$n)
# 1    2
# 1047   50
# 50 pieces had 2 scg annotations, while 1047 contained one. 

mean(annots_level1_scg$size)
# [1] 2342.132

sd(annots_level1_scg$size)
# [1] 1107.838

tmp <- annots_level1_scg %>%
  select(dom_id, KEGG_ko2) %>%
  distinct() %>%
  group_by(KEGG_ko2) %>%
  tally() %>%
  arrange(desc(n))
mean(tmp$n)
# [1] 17.64615
# each of the 76 marker genes is annotated to 17.6 nodes (min 3, max 66)
# average number of annotations per piece? -------------------------------

# distribution of piece size ---------------------------------------------


