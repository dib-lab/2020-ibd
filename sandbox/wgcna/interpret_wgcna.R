library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(clusterProfiler)
mods <- read_csv("~/Downloads/wgcna_protein_modules.csv") %>%
  mutate(gene = gsub(".fa", ".fna", gene))
mods <- separate(mods, col = gene, sep = ".fna_", into = c("genome", "protein"))

# how many modules do sig diff prots from each genome separate out into?
num_mods_per_genome <- mods %>% 
  group_by(genome, moduleColor) %>%
  count() %>% 
  group_by(genome) %>%
  tally()
ggplot(num_mods_per_genome, aes(x = n)) +
  geom_density() +
  theme_minimal()

# how many genomes are in each module?
num_genomes_per_mod <- mods %>% 
  group_by(genome, moduleColor) %>%
  count() %>% 
  group_by(moduleColor) %>%
  tally()
colnames(num_genomes_per_mod) <- c("moduleColor", "n_genomes_per_mod")
# how many genes are there per module?
num_prots_per_mod <- mods %>%
  group_by(moduleColor) %>%
  tally() 
colnames(num_prots_per_mod) <- c("moduleColor", "n_prots_per_mod")

mod_specs <- full_join(num_genomes_per_mod, num_prots_per_mod)

# Do clusters correspond to KO functions?
nog <- read_tsv("~/Downloads/job_MM_2y3lv21c_annotations2.txt", comment = "#", 
                col_names = c('query_name', 'seed_eggNOG_ortholog',	'seed_ortholog_evalue'	,
                              'seed_ortholog_score', 'best_tax_level', 'Preferred_name',
                              'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                              'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC',
                              'CAZy', 'BiGG_Reaction')) %>%
  mutate(query_name = gsub(".fa_", ".fna_", query_name))
nog <- separate(nog, col = query_name, into = c("genome", "protein"), 
                sep = ".fna_sig_ccs.faa.nostop.faa_", remove = T)

all <- full_join(mods, nog, by = c("genome", "protein"))

# how many KOs are there per module?
n_kegg_per_mod <- all %>% 
  filter(!is.na(KEGG_ko)) %>%
  group_by(moduleColor) %>%
  tally()
colnames(n_kegg_per_mod) <- c("moduleColor", "n_kegg_per_mod")

mod_specs <- full_join(mod_specs, n_kegg_per_mod)
mod_specs$perc_ko <- mod_specs$n_kegg_per_mod / mod_specs$n_prots_per_mod 
  
# are KOs enriched within a module?
pathways <- read_tsv("http://rest.kegg.jp/link/pathway/ko", col_names = c("KO", "path")) %>%
  mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

pathway_names <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))

pink <- all %>% 
  filter(moduleColor == "pink") %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))
enriched <- enricher(gene = pink$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(enriched)


pink <- all %>% 
  filter(moduleColor == "pink") %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))
enriched <- enricher(gene = pink$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(enriched)


thistle2<- all %>% 
  filter(moduleColor == "thistle2") %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))
enriched <- enricher(gene = thistle2$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(enriched)

grey60 <- all %>% 
  filter(moduleColor == "grey60") %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))
enriched <- enricher(gene = grey60$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(enriched)


grey60 <- all %>% 
  filter(moduleColor == "grey60") %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))
enriched <- enricher(gene = grey60$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(enriched)


black <- all %>% 
  filter(moduleColor == "black") %>%
  filter(!is.na(KEGG_ko)) %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))
enriched <- enricher(gene = black$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(enriched)


blue <- all %>% 
  filter(moduleColor == "blue") %>%
  filter(!is.na(KEGG_ko)) %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))
enriched <- enricher(gene = blue$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(enriched)


table(black$KEGG_ko %in% blue$KEGG_ko)
table(blue$KEGG_ko %in% black$KEGG_ko)
intersect(blue$KEGG_ko, black$KEGG_ko)

# are KOs that we observe multiple times more likely to be in the same cluster than in different clusters?
# need to know how many times we observe each KO, and how many modules we observe each KO in
num_mods_per_kegg_ko <- all %>%
  group_by(KEGG_ko, moduleColor) %>%
  tally() %>%
  group_by(KEGG_ko) %>%
  tally()
colnames(num_mods_per_kegg_ko) <- c("KEGG_ko", "num_mods")
num_each_kegg_ko <- all %>%
  group_by(KEGG_ko) %>%
  tally() 
colnames(num_each_kegg_ko) <- c("KEGG_ko", "num_times_observed")
ko_specs <- full_join(num_each_kegg_ko, num_mods_per_kegg_ko)
View(ko_specs)
