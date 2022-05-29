library(dplyr)
library(tidyr)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, diagnosis) %>%
  mutate(library_name = gsub("-", "\\.", library_name)) %>%
  distinct()
table(info$library_name %in% all$X1)

all <- read_tsv("sandbox/at_least_5_studies_vita_vars.tsv")
all <- left_join(info, all, by = c("library_name" = "X1"))

head(colnames(all))

ggplot(all, aes(x = diagnosis, y = `4121156555200726`*100)) +
  geom_violin() +
  theme_minimal()

all_long_summarized <- all %>%
  pivot_longer(cols = `1119237125513`:`9141524607182576`, names_to = "hash",
               values_to = "norm_abund") %>%
  group_by(diagnosis, hash) %>%
  mutate(norm_abund = norm_abund * 100) %>%
  summarize(mean = mean(norm_abund),
            sd = sd(norm_abund)) 
View(all_long_summarized)

all_wide_summarized <- all_long_summarized %>%
  select(-sd) %>%
  pivot_wider(names_from = diagnosis, values_from = mean)

cd_up <- all_wide_summarized %>%
  filter(CD > nonIBD)%>%
  select(-UC) %>%
  mutate(diff = nonIBD - CD)

cd_down <- all_wide_summarized %>%
  filter(CD < nonIBD) %>%
  select(-UC) %>%
  mutate(diff = nonIBD - CD)


# enrichment --------------------------------------------------------------

# download kegg dbs
pathways <- read_tsv("http://rest.kegg.jp/link/pathway/ko", col_names = c("KO", "path"))
pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

pathway_names <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))

# join with multifasta query results

sgc <- read_csv("outputs/gather_matches_loso_multifasta/all-multifasta-query-results.csv") %>%
  filter(hashval != "hashval") %>%
  separate(record_name, sep = " ", into = c("sequence_name", "prokka_annotation"), extra = "merge") %>%
  mutate(hashval = as.character(hashval))

# read in eggnog results
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
eggnog$KEGG_ko <- gsub(",.*", "", eggnog$KEGG_ko) # select first KO annotation

# CD DOWN
cd_down <- left_join(cd_down, sgc, by = c("hash" = "hashval"))
cd_down <- left_join(cd_down, eggnog, by = c("sequence_name" = "query_name"))

cd_down2 <- cd_down %>%
  filter(!is.na(KEGG_ko)) %>%
  mutate(gather_genome = gsub(".fna.*", "", sequence_name)) %>%
  mutate(gather_genome = gsub(".fa.*", "", gather_genome))%>%
  mutate(gather_genome = gsub("_232_.*", "", gather_genome)) %>%
  select(-catlas_base) %>%
  distinct()

cd_down3 <- cd_down2 %>%
  select(hash, gather_genome) %>%
  distinct() %>%
  group_by(gather_genome) %>%
  tally() %>%
  mutate(percent = n/sum(n)) %>%
  mutate(set = "down")


View(cd_down3)
cd_down <- cd_down %>%
  filter(!is.na(KEGG_ko)) %>%
  select(hash, CD, nonIBD, diff, KEGG_ko) %>%
  distinct() 

cd_down_enriched<- enricher(gene = cd_down$KEGG_ko, 
                    TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(cd_down_enriched)

# CD UP
cd_up <- left_join(cd_up, sgc, by = c("hash" = "hashval"))
cd_up <- left_join(cd_up, eggnog, by = c("sequence_name" = "query_name"))

cd_up2 <- cd_up %>%
  filter(!is.na(KEGG_ko)) %>%
  mutate(gather_genome = gsub(".fna.*", "", sequence_name)) %>%
  mutate(gather_genome = gsub(".fa.*", "", gather_genome))%>%
  mutate(gather_genome = gsub("_232_.*", "", gather_genome)) %>%
  select(-catlas_base) %>%
  distinct()

cd_up3 <- cd_up2 %>%
  select(hash, gather_genome) %>%
  distinct() %>%
  group_by(gather_genome) %>%
  tally() %>%
  mutate(percent = n/sum(n)) %>%
  mutate(set = "up")

tmp <- rbind(cd_up3, cd_down3)
ggplot(tmp, aes(x = gather_genome, y = n, fill = set)) +
  geom_col() +
  theme_minimal()

cd_up_small <- cd_up %>%
  filter(!is.na(KEGG_ko)) %>%
  select(hash, CD, nonIBD, diff, KEGG_ko) %>%
  distinct()

cd_up_enriched<- enricher(gene = cd_up_small$KEGG_ko, 
                            TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(cd_up_enriched)
