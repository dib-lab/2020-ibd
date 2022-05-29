library(dplyr)
library(readr)
library(clusterProfiler)
setwd("~/github/2020-ibd")

sig_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_no_cdbg.tsv")

table(sig_ccs$mu, sig_ccs$level)

# limit to level 1 --------------------------------------------------------


sig_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_no_cdbg.tsv") %>%
  filter(level == 1)

# how many differentially abundant dominators don't have an annotation?

# prokka query annotation:
sig_ccs %>%
  select(dom_id, record_id) %>%
  filter(is.na(record_id)) %>%
  distinct()
# 2,436 don't have prokka annotations of 14,446 dom_ids
sig_ccs %>%
  select(dom_id) %>%
  distinct() %>%
  nrow()

# what is the average number of annotations a dominator has? -------------
sig_ccs %>% 
  group_by(dom_id) %>%
  tally() %>%
  summarize(mean = mean(n)) 
# 1.19

View(sig_ccs %>% 
  group_by(dom_id) %>%
  tally())
# the dominator with 8 matches a tRNA

# how many increased/decreased in UC/CD? ---------------------------------
sig_ccs$direction <-  ifelse(sig_ccs$estimate > 0, "increased", "decreased")
table(sig_ccs$mu, sig_ccs$direction)

#                  decreased increased
# mu.diagnosisCD      8832      7877
# mu.diagnosisUC       462         0

# what is the kegg enrichment for things that are annotated? -------------

pathways <- "http://rest.kegg.jp/link/pathway/ko"
download.file(url = pathways, destfile = "inputs/kegg_pathways.tsv")
pathways <- read_tsv("inputs/kegg_pathways.tsv", col_names = c("KO", "path"))
pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

pathway_names <- "http://rest.kegg.jp/list/pathway"
download.file(url = pathway_names, destfile = "inputs/kegg_pathway_names.tsv")
pathway_names <- read_tsv("inputs/kegg_pathway_names.tsv", col_names = c("path", "name"))

cd_up_kegg <- sig_ccs %>%
  filter(direction == "increased") %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(!is.na(KEGG_ko))

cd_up_kegg_res <- enricher(gene = cd_up_kegg$KEGG_ko, 
                           TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(cd_up_kegg_res)

cd_down_kegg <- sig_ccs %>%
  filter(direction == "decreased") %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(!is.na(KEGG_ko))

cd_down_kegg_res <- enricher(gene = cd_down_kegg$KEGG_ko, 
                             TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(cd_down_kegg_res)

# a lot of overlap -- is there overlap in the genes themselves?
table(cd_down_kegg$record_id %in% cd_up_kegg$record_id)

# check for not just kegg, but all of the records: 

cd_down <- sig_ccs %>%
  filter(direction == "decreased") %>%
  filter(mu == "mu.diagnosisCD")

cd_up <- sig_ccs %>%
  filter(direction == "increased") %>%
  filter(mu == "mu.diagnosisCD")

table(unique(cd_down$record_id) %in% unique(cd_up$record_id))
# FALSE  TRUE 
# 487    27 
table(unique(cd_up$record_id) %in% unique(cd_down$record_id))
# FALSE  TRUE 
# 2255    27

# what are the different genes? -------------------------------------------

cd_up_only <- cd_up %>%
  filter(! record_id %in% cd_down$record_id)

cd_up_only_kegg_res <- enricher(gene = cd_up_only$KEGG_ko, 
                                TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(cd_up_only_kegg_res)
View(cd_up_only)

ox_ko<- c('K00362', 'K11717', 'K01926', 'K01069', 'K03671', 'K04565', 'K14155', 
          'K11717', 'K01069', 'K14155', 'K03386', 'K01758', 'K03564', 'K00362', 
          'K01926', 'K04565', 'K03676', 'K01759', 'K03671', 'K01758', 'K01759', 
          'K01758', 'K14155')
nox_ko <- c('K00362', 'K05601', 'K03386', 'K00362', 'K05601')

View(cd_up_only %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  filter(KEGG_ko %in% ox_ko))
