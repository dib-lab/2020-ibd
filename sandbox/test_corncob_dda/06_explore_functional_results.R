library(dplyr)
library(readr)
library(clusterProfiler)
setwd("~/github/2020-ibd")

sig_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_roary_isolates_and_megahit_no_cdbg.tsv")

table(sig_ccs$mu, sig_ccs$level)

# limit to level 1 --------------------------------------------------------


sig_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_roary_isolates_and_megahit_no_cdbg.tsv") %>%
  filter(level == 1)

# how many differentially abundant dominators don't have an annotation?

# prokka query annotation:
sig_ccs %>%
  select(dom_id, record_id) %>%
  filter(is.na(record_id)) %>%
  distinct() %>%
  nrow()
# 1242 (down from 2446) don't have prokka annotations of 14,446 dom_ids
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

# what is the kegg enrichment for things that are annotated? -------------

pathways <- read_tsv("http://rest.kegg.jp/link/pathway/ko", col_names = c("KO", "path"))
pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

pathway_names <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))

cd_up_kegg <- sig_ccs %>%
  filter(direction == "increased") %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(!is.na(KEGG_ko))

cd_up_kegg_res <- enricher(gene = cd_up_kegg$KEGG_ko, 
                           TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(cd_up_kegg_res, showCategory =100)

cd_down_kegg <- sig_ccs %>%
  filter(direction == "decreased") %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(!is.na(KEGG_ko))

cd_down_kegg_res <- enricher(gene = cd_down_kegg$KEGG_ko, 
                             TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(cd_down_kegg_res, showCategory =100)

# some overlap -- is there overlap in the genes themselves?
table(cd_down_kegg$record_id %in% cd_up_kegg$record_id)
# FALSE  TRUE 
# 8630   184 
# check for not just kegg, but all of the records: 

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
