library(clusterProfiler)
library(readr)

setwd("~/Downloads")

## download ortholog to pathway table

pathways <- read_tsv("http://rest.kegg.jp/link/pathway/ko", col_names = c("KO", "path"))
pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

pathway_names <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))


eggnog <- read_tsv("outputs/gather_matches_loso_multifasta/all-multifasta-query-results.emapper.annotations", 
                   comment = "#", 
                   col_names = c('query_name', 'seed_eggNOG_ortholog',	
                                 'seed_ortholog_evalue', 'seed_ortholog_score',
                                 'best_tax_level', 'Preferred_name', 'GOs', 'EC',
                                 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',	'KEGG_TC',
                                 'CAZy', 'BiGG_Reaction', 'taxonomic_scope', 
                                 'eggNOG_OGs', 'best_eggNOG_OG', 'COG_functional_cat',	
                                 'eggNOG', 'free_text_desc'))
tmp <- gsub(",.*", "", eggnog$KEGG_ko) # select first KO annotation
#tmp <- gsub("ko:", "", tmp)

# enrichment --------------------------------------------------------------


enriched<- enricher(gene = tmp[!is.na(tmp)], 
                    TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(enriched)

