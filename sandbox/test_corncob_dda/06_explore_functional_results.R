library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(clusterProfiler)
setwd("~/github/2020-ibd")

# how many tested at each level? ------------------------------------------
all_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_all_ccs.tsv")
dom_info <- read_tsv("sandbox/test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_dom_info.tsv")
all_ccs <- left_join(all_ccs, dom_info, by = c("aa_seq" = "dom_id"))

all_ccs %>% 
  select(aa_seq, level) %>%
  distinct() %>%
  group_by(level) %>%
  tally()
# row level     n
# 1     1 25897
# 2     2 10095
# 3     3  4049
# 4     4  1664
# 5     5   724
# 6     6   315
# 7     7   123
# 8     8    41
# 9     9    29
# 10    10    29
# 11    11     1

dom_info %>%
  group_by(level) %>%
  tally()
# how many sig at each level? ---------------------------------------------
sig_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_roary_isolates_and_megahit_no_cdbg.tsv")

sig_ccs %>%
  select(mu, dom_id, level) %>%
  distinct() %>%
  group_by(mu, level) %>%
  tally()

# limit to level 1 --------------------------------------------------------


sig_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_roary_isolates_and_megahit_no_cdbg.tsv") %>%
  filter(level == 1) #%>%
  #filter(mu == "mu.diagnosisCD")

# how many differentially abundant dominators don't have an annotation?

# prokka query annotation:
sig_ccs %>%
  select(dom_id, record_id) %>%
  filter(is.na(record_id)) %>%
  distinct() %>%
  nrow()
# 1242 (down from 2446) don't have prokka annotations of 14,446 dom_ids
# 1239 don't have prokka annotations of 14,443 for CD
sig_ccs %>%
  select(dom_id) %>%
  distinct() %>%
  nrow()

# what is the average number of annotations a dominator has? -------------
sig_ccs %>% 
  group_by(dom_id) %>%
  tally() %>%
  summarize(mean = mean(n)) 
# 5.23 (up from 1.19) -- probably potential duplicates from using a pangenome reference

View(sig_ccs %>% 
  group_by(dom_id) %>%
  tally())
# many annotations for top dom id

# how many increased/decreased in UC/CD? ---------------------------------
sig_ccs$direction <-  ifelse(sig_ccs$estimate > 0, "increased", "decreased")
table(sig_ccs$mu, sig_ccs$direction)

#                  decreased increased
# mu.diagnosisCD      48647      24780
# mu.diagnosisUC       2100         0
# these numbers are inflated by the annottions 

tmp <- sig_ccs %>%
  select(dom_id, mu, direction) %>%
  distinct()
table(tmp$mu, tmp$direction)
#                 decreased increased
# mu.diagnosisCD      7955      6478
# mu.diagnosisUC       441         0


# just mu CD --------------------------------------------------------------

sig_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_roary_isolates_and_megahit_no_cdbg.tsv") %>%
  filter(level == 1) %>%
  filter(mu == "mu.diagnosisCD") %>%
  mutate(direction = ifelse(estimate > 0, "increased", "decreased"))

# prokka query annotation:
sig_ccs %>%
  select(dom_id, record_id, direction) %>%
  filter(!is.na(record_id)) %>%
  select(dom_id, direction) %>%
  distinct() %>%
  group_by(direction) %>%
  tally()
# distinct prokka annotations: 
sig_ccs %>%
  select(dom_id, record_id, direction) %>%
  filter(!is.na(record_id)) %>%
  select(record_id, direction) %>%
  distinct() %>%
  group_by(direction) %>%
  tally()

# KO query annotation:
sig_ccs %>%
  select(dom_id, KEGG_ko, direction) %>%
  arrange(KEGG_ko) %>%
  group_by(dom_id) %>%
  slice(1) %>%
  filter(!is.na(KEGG_ko)) %>%
  ungroup() %>%
  group_by(direction) %>%
  tally()

# distinct KO annotations:
sig_ccs %>%
  select(dom_id, KEGG_ko, direction) %>%
  filter(!is.na(KEGG_ko)) %>%
  select(KEGG_ko, direction) %>%
  distinct() %>%
  group_by(direction) %>%
  tally()

# eggnog query annotation:
sig_ccs %>%
  select(dom_id, `eggNOG OGs`, direction) %>%
  arrange(`eggNOG OGs`) %>%
  group_by(dom_id) %>%
  slice(1) %>%
  filter(!is.na(`eggNOG OGs`)) %>%
  ungroup() %>%
  group_by(direction) %>%
  tally()

# distinct eggnog annotations:
sig_ccs %>%
  select(dom_id, `eggNOG OGs`, direction) %>%
  filter(!is.na(`eggNOG OGs`)) %>%
  select(`eggNOG OGs`, direction) %>%
  distinct() %>%
  group_by(direction) %>%
  tally()

# what is the kegg enrichment for things that are annotated? -------------
pathways <- read_tsv("http://rest.kegg.jp/link/pathway/ko", col_names = c("KO", "path"))
pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

pathway_names <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))

# need to select the first KO that occurs
cd_up_kegg <- sig_ccs %>%
  filter(direction == "increased") %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(!is.na(KEGG_ko)) %>%
  mutate(KEGG_ko = gsub(",.*", "", KEGG_ko))

cd_up_kegg_res <- enricher(gene = cd_up_kegg$KEGG_ko, 
                           TERM2GENE = pathways, TERM2NAME = pathway_names)
cd_up_kegg_plt <-dotplot(cd_up_kegg_res, showCategory =100, font.size = 8,
                         title = "increased in CD") +
  geom_point(color = "black") +
  theme(legend.position = "bottom") +
  guides(color = FALSE)

cd_down_kegg <- sig_ccs %>%
  filter(direction == "decreased") %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(!is.na(KEGG_ko))  %>%
  mutate(KEGG_ko = gsub(",.*", "", KEGG_ko))

cd_down_kegg_res <- enricher(gene = cd_down_kegg$KEGG_ko, 
                             TERM2GENE = pathways, TERM2NAME = pathway_names)
cd_down_kegg_plt <- dotplot(cd_down_kegg_res, showCategory =100, font.size = 7.5, 
                            title = "decreased in CD") +
  geom_point(color = "black") +
  theme(legend.position = "bottom") +
  guides(color = FALSE)

cowplot::plot_grid(cd_up_kegg_plt, cd_down_kegg_plt, labels = c("A", "B"))


# hunting down KEGG overlaps ---------------------------------------------------
# some overlap -- is there overlap in the genes themselves?
table(cd_down_kegg$record_id %in% cd_up_kegg$record_id)
# FALSE  TRUE 
# 8630   184 
table(cd_up_kegg$record_id %in% cd_down_kegg$record_id)
# FALSE  TRUE 
# 12195   172

table(unique(cd_down_kegg$KEGG_ko) %in% unique(cd_up_kegg$KEGG_ko))
# FALSE  TRUE 
# 84     222
table(unique(cd_up_kegg$KEGG_ko) %in% unique(cd_down_kegg$KEGG_ko))
# FALSE  TRUE 
# 920    222 

library(KEGGREST)
ko <- keggList("ko")
ko <- ko %>%
  as.data.frame() %>%
  rownames_to_column("KO")

overlapping_kos <- cd_up_kegg[cd_up_kegg$KEGG_ko %in% cd_down_kegg$KEGG_ko, ]
overlapping_kos <- overlapping_kos %>%
  select(KEGG_ko, KEGG_Pathway, KEGG_Module, KEGG_Reaction, KEGG_rclass, BRITE,
         KEGG_TC, KEGG_Reaction) %>%
  distinct() %>%
  left_join(ko, by = c("KEGG_ko" = "KO"))

scg <- read_csv("inputs/kegg_scg_13059_2015_610_MOESM1_ESM.csv", skip =2)
table(scg$KO %in% gsub("ko:", "", overlapping_kos$KEGG_ko))
# FALSE  TRUE 
# 57    19 

overlapping_kos_res <- enricher(gene = overlapping_kos$KEGG_ko, 
                                TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(overlapping_kos_res, showCategory = 100, font.size = 7.5) +
  geom_point(color = "black") +
  theme(legend.position = "bottom") +
  guides(color = FALSE)

# hunting down all overlaps ---------------------------------------------------

cd_down <- sig_ccs %>%
  filter(direction == "decreased") %>%
  filter(mu == "mu.diagnosisCD")

cd_up <- sig_ccs %>%
  filter(direction == "increased") %>%
  filter(mu == "mu.diagnosisCD")

table(unique(cd_down$record_id) %in% unique(cd_up$record_id))
# FALSE  TRUE 
# 6361   119
table(unique(cd_up$record_id) %in% unique(cd_down$record_id))
# FALSE  TRUE 
# 10001   119 

# what are the different genes? -------------------------------------------

cd_up_only <- cd_up %>%
  filter(! record_id %in% cd_down$record_id)

cd_up_only_kegg_res <- enricher(gene = cd_up_only$KEGG_ko, 
                                TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(cd_up_only_kegg_res, showCategory = 100)

ox_ko<- c('K00362', 'K11717', 'K01926', 'K01069', 'K03671', 'K04565', 'K14155', 
          'K11717', 'K01069', 'K14155', 'K03386', 'K01758', 'K03564', 'K00362', 
          'K01926', 'K04565', 'K03676', 'K01759', 'K03671', 'K01758', 'K01759', 
          'K01758', 'K14155')
nox_ko <- c('K00362', 'K05601', 'K03386', 'K00362', 'K05601')

View(cd_up_only %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  filter(KEGG_ko %in% ox_ko))
