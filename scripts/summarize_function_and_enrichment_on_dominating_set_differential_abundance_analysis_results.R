library(dplyr)
library(readr)
library(clusterProfiler)

# read in and separate sig diff abund dom set pieces ---------------------

sig_ccs <- read_tsv(snakemake@input[['sig_ccs_annot']]) %>%
  filter(mu == "mu.diagnosisCD") # focus only on CD

sig_ccs_up <- sig_ccs %>%
  filter(estimate > 0) %>%
  mutate(KEGG_ko = gsub(",.*", "", KEGG_ko)) %>%
  select(dom_id, KEGG_ko) %>%
  filter(!is.na(KEGG_ko)) %>%
  filter(KEGG_ko != "-") %>%
  distinct()

sig_ccs_down <- sig_ccs %>%
  filter(estimate < 0) %>%
  mutate(KEGG_ko = gsub(",.*", "", KEGG_ko)) %>%
  select(dom_id, KEGG_ko) %>%
  filter(!is.na(KEGG_ko)) %>%
  filter(KEGG_ko != "-") %>%
  distinct()

sig_ccs_up_only <- sig_ccs_up %>%
  filter(!KEGG_ko %in% sig_ccs_down$KEGG_ko) %>%
  mutate(abundance = "increased",
         species = snakemake@wildcards[['gtdb_species']])

sig_ccs_down_only <- sig_ccs_down %>%
  filter(!KEGG_ko %in% sig_ccs_up$KEGG_ko) %>%
  mutate(abundance = "decreased",
         species = snakemake@wildcards[['gtdb_species']])

sig_ccs_distinct <- bind_rows(sig_ccs_up_only, sig_ccs_down_only)
write_tsv(sig_ccs_distinct, snakemake@output[['distinct_df']])

# read gene sets of interest ----------------------------------------------

ox_ko <- read_tsv(snakemake@input[['ox_kos']])
marker_ko <- read_csv(snakemake@input[['marker_kos']])

# definitions <- read_csv('~/github/2020-ibd/inputs/kegg_ortholog_definitions.csv')
# kegg_scg <- read_csv("~/github/2020-ibd/inputs/kegg_scg_13059_2015_610_MOESM1_ESM.csv") %>%
#   janitor::clean_names()%>%
#   left_join(definitions) %>%
#   filter(average_copy_number == 1) %>%
#   distinct()
# write_csv(kegg_scg, "~/github/2020-ibd/inputs/kegg_scg_13059_2015_610_MOESM1_ESM_filtered.csv")

# read in KEGG pathway metadata for enrichment ----------------------------

## download ortholog to pathway table
pathways <- "http://rest.kegg.jp/link/pathway/ko"
if(!file.exists("inputs/kegg_pathways.tsv")){
  download.file(url = pathways, destfile = "inputs/kegg_pathways.tsv")
}
pathways <- read_tsv("inputs/kegg_pathways.tsv", col_names = c("KO", "path"))
pathways <- pathways %>%
  filter(grepl(pattern = "map", path)) %>%
  dplyr::select(path, KO)

pathway_names <- "http://rest.kegg.jp/list/pathway"
if(!file.exists("inputs/kegg_pathway_names.tsv")){
  download.file(url = pathway_names, destfile = "inputs/kegg_pathway_names.tsv")
}
pathway_names <- read_tsv("inputs/kegg_pathway_names.tsv", col_names = c("path", "name"))

# distinct kegg enrichment analysis ----------------------------------------

enriched_up <- enricher(gene = sig_ccs_up_only$KEGG_ko, 
                        TERM2GENE = pathways, TERM2NAME = pathway_names, maxGSSize = 500)
enriched_up_df <- enriched_up@result %>%
  mutate(species = snakemake@wildcards[['gtdb_species']],
         abundance = "increased")

enriched_down <- enricher(gene = sig_ccs_down_only$KEGG_ko, 
                        TERM2GENE = pathways, TERM2NAME = pathway_names, maxGSSize = 500)
enriched_down_df <- enriched_down@result %>%
  mutate(species = snakemake@wildcards[['gtdb_species']],
         abundance = "decreased")

enriched_df <- bind_rows(enriched_up_df, enriched_down_df)

write_tsv(enriched_df, snakemake@output[['enriched_distinct']])

# distinct ox -------------------------------------------------------------

ox_up_only <- sig_ccs_up_only %>%
  mutate(ko = gsub("ko:", "", KEGG_ko),
         abundance = "increased",
         species = snakemake@wildcards[['gtdb_species']]) %>%
  inner_join(ox_ko, by = "ko") 

ox_down_only <- sig_ccs_down_only %>%
  mutate(ko = gsub("ko:", "", KEGG_ko),
         abundance = "decreased",
         species = snakemake@wildcards[['gtdb_species']]) %>%
  inner_join(ox_ko, by = "ko") 

ox_only <- bind_rows(ox_up_only, ox_down_only)
write_tsv(ox_only, snakemake@output[['distinct_ox']])

# overlap analysis --------------------------------------------------------

overlapping_kos <- sig_ccs_up %>%
  filter(KEGG_ko %in% sig_ccs_down$KEGG_ko) %>%
  select(KEGG_ko)

overlapping_kos_df <- sig_ccs %>%
  filter(KEGG_ko %in% overlapping_kos$KEGG_ko) %>%
  mutate(species = snakemake@wildcards[['gtdb_species']]) 

write_tsv(overlapping_kos_df, snakemake@output[['overlap_df']])

enriched_overlapping <- enricher(gene = overlapping_kos$KEGG_ko, 
                                 TERM2GENE = pathways, TERM2NAME = pathway_names, maxGSSize = 500)

enriched_overlapping_df <- enriched_overlapping@result %>%
  mutate(species = snakemake@wildcards[['gtdb_species']],
         abundance = "overlap")

write_tsv(enriched_overlapping_df, snakemake@output[['enriched_overlap']])

overlapping_marker_kos <-  sig_ccs %>%
  mutate(KEGG_ko = gsub(",.*", "", KEGG_ko),
         filename = basename(filename)) %>%
  filter(gsub("ko:", "", KEGG_ko) %in% marker_ko$ko) %>%
  filter(KEGG_ko %in% overlapping_kos$KEGG_ko) %>%
  arrange(desc(estimate)) %>%
  arrange(record_id) 


write_tsv(overlapping_marker_kos, snakemake@output[['overlap_marker_kos']])
