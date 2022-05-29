library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(rjson)
library(clusterProfiler)

# # read in library info
# info <- read_tsv("inputs/working_metadata.tsv") %>%
#   select(library_name, diagnosis) %>%
#   mutate(library_name = gsub("-", "\\.", library_name)) %>%
#   distinct()

eggnog_files <- list.files("sandbox/test_megahit_diginorm_nocat/eggnog", 
                           ".annotations$", full.names = T)



eggnog_pan <- eggnog_files %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, .id = "genome", comment = "#", 
          col_names = c('query_name', 'seed_eggNOG_ortholog',	
                        'seed_ortholog_evalue', 'seed_ortholog_score',
                        'best_tax_level', 'Preferred_name', 'GOs', 'EC',
                        'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                        'KEGG_Reaction', 'KEGG_rclass', 'BRITE',	'KEGG_TC',
                        'CAZy', 'BiGG_Reaction', 'taxonomic_scope', 
                        'eggNOG_OGs', 'best_eggNOG_OG', 'COG_functional_cat',	
                        'eggNOG')) %>%
  mutate(genome = gsub(".emapper.annotations", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/eggnog/", "", genome))

eggnog_pan$KEGG_ko <- gsub(",.*", "", eggnog_pan$KEGG_ko) # select first KO annotation


corncob_files <- list.files("sandbox/test_megahit_diginorm_nocat/corncob",
                            ".sig_ccs.tsv$", full.names = T) 
corncob_pan <- corncob_files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genome") %>%
  mutate(genome = gsub("_sig_ccs.tsv", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/corncob/", "", genome)) 

pan_results <- left_join(corncob_pan, eggnog_pan, by = c("genome", "aa_seq" = "query_name")) %>%
  mutate(genome = gsub(".fa", "", genome)) %>%
  mutate(genome = gsub(".fna", "", genome))

species <- read_tsv("sandbox/test_megahit_diginorm_nocat/gtdbtk/41genomes_out/gtdbtk.bac120.summary.tsv") %>%
  select(user_genome, classification) %>%
  mutate(classification = gsub("[dkpcofgs]__", "", classification)) %>%
  separate(data = ., col = classification, sep = ";",
           into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))

pan_results <- left_join(pan_results, species, by = c("genome" = "user_genome"))


# check up down -----------------------------------------------------------

pan_results_up_down <- pan_results %>%
  mutate(direction = ifelse(estimate < 0, "down", "up")) %>%
  group_by(genome, direction, mu) %>%
  tally()
View(pan_results_up_down)

ggplot(pan_results_up_down, aes(x = genome, y = n, fill = direction)) +
  geom_col() +
  theme_minimal() +
  facet_wrap(~mu, scales = "free_x") +
  coord_flip()

# pathway -----------------------------------------------------------------
kegg_def <- read_csv("inputs/kegg_ortholog_definitions.csv")

pathways <- read_delim("http://rest.kegg.jp/link/pathway/ko", delim = "\t", col_names = c("KO", "path"))

pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)

pathway_names <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))



#### R gnavus
rgnv_cd_up_gene <- pan_results %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate > 0) %>%
  filter(genome == "GCF_900036035.1_RGNV35913_genomic")

rgnv_cd_down_gene <- pan_results %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate < 0) %>%
  filter(genome == "GCF_900036035.1_RGNV35913_genomic")


rgnv_cd_up_not_in_cd_down <- rgnv_cd_up_gene %>%
  filter(!KEGG_ko %in% rgnv_cd_down_gene$KEGG_ko)
rgnv_cd_up_not_in_cd_down_enriched<- enricher(gene = rgnv_cd_up_not_in_cd_down$KEGG_ko, 
                                               TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(rgnv_cd_up_not_in_cd_down_enriched)
View(rgnv_cd_up_not_in_cd_down_enriched@result %>%
       filter(p.adjust < .05))

### C bolt
cbolt_cd_up_gene <- pan_results %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate > 0) %>%
  filter(genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic")

cbolt_cd_down_gene <- pan_results %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate < 0) %>%
  filter(genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic")

cbolt_cd_up_not_in_cd_down <- cbolt_cd_up_gene %>%
  filter(!KEGG_ko %in% cbolt_cd_down_gene$KEGG_ko)
cbolt_cd_up_not_in_cd_down_enriched<- enricher(gene = cbolt_cd_up_not_in_cd_down$KEGG_ko, 
                                       TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(cbolt_cd_up_not_in_cd_down_enriched)
View(cbolt_cd_up_not_in_cd_down_enriched@result %>%
       filter(p.adjust < .05))

### Flavronifractor plautii

fplaut_cd_up_gene <- pan_results %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate > 0) %>%
  filter(genome == "ERS537353_12")

fplaut_cd_down_gene <- pan_results %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate < 0) %>%
  filter(genome == "ERS537353_12")

fplaut_cd_up_not_in_cd_down <- fplaut_cd_up_gene %>%
  filter(!KEGG_ko %in% fplaut_cd_down_gene$KEGG_ko)
fplaut_cd_up_not_in_cd_down_enriched<- enricher(gene = fplaut_cd_up_not_in_cd_down$KEGG_ko, 
                                               TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(fplaut_cd_up_not_in_cd_down_enriched)
View(fplaut_cd_up_not_in_cd_down_enriched@result %>%
       filter(p.adjust < .05))


### Flavronifractor sp

fsp_cd_up_gene <- pan_results %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate > 0) %>%
  filter(genome == "GCF_000508885.1_ASM50888v1_genomic")

fsp_cd_down_gene <- pan_results %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate < 0) %>%
  filter(genome == "GCF_000508885.1_ASM50888v1_genomic")

fsp_cd_up_not_in_cd_down <- fsp_cd_up_gene %>%
  filter(!KEGG_ko %in% fsp_cd_down_gene$KEGG_ko)
fsp_cd_up_not_in_cd_down_enriched<- enricher(gene = fsp_cd_up_not_in_cd_down$KEGG_ko, 
                                                TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(fsp_cd_up_not_in_cd_down_enriched)
View(fsp_cd_up_not_in_cd_down_enriched@result %>%
       filter(p.adjust < .05))

# plot per species enrichment on single plot ------------------------

# crohn's disease

cd_up_enriched_orgs <- fsp_cd_up_not_in_cd_down_enriched@result %>%
  filter(qvalue < 0.05) %>%
  select(Description, GeneRatio, Count) %>%
  mutate(genome = "Flavonifractor sp000508885")

cd_up_enriched_orgs <- cbolt_cd_up_not_in_cd_down_enriched@result %>%
  filter(qvalue < 0.05) %>%
  select(Description, GeneRatio, Count) %>%
  mutate(genome = "Clostridium_M bolteae") %>%
  rbind(cd_up_enriched_orgs)

cd_up_enriched_orgs <- fplaut_cd_up_not_in_cd_down_enriched@result %>%
  filter(qvalue < 0.05) %>%
  select(Description, GeneRatio, Count) %>%
  mutate(genome = "Flavonifractor plautii")%>%
  rbind(cd_up_enriched_orgs)


cd_up_enriched_orgs <- rgnv_cd_up_not_in_cd_down_enriched@result %>%
  filter(qvalue < 0.05) %>%
  select(Description, GeneRatio, Count) %>%
  mutate(genome = "Faecalicatena gnavus")%>%
  rbind(cd_up_enriched_orgs)

common_descriptions <- cd_up_enriched_orgs %>% 
  group_by(Description) %>% 
  tally() %>% 
  filter(n >= 2) %>% 
  ungroup()

cd_up_enriched_orgs <- cd_up_enriched_orgs %>% 
  filter(Description %in% common_descriptions$Description) %>%
  separate(GeneRatio, into = c("Count", "Total"), sep = "\\/") %>%
  mutate(GeneRatio = as.numeric(Count)/as.numeric(Total)) %>%
  mutate(Description = gsub("biosynthesis", "biosynth", Description))
  
ggplot(cd_up_enriched_orgs, aes(x = reorder(Description, GeneRatio), y = GeneRatio, 
                                color = genome, size = as.numeric(Count))) +
  geom_point(alpha = .5) +
  theme_minimal() +
  coord_flip()+
  theme(legend.position = "right",
        axis.title.y = element_blank()) +
  labs(y = "gene ratio", size = "count") + 
  guides(colour = guide_legend(override.aes = list(size=4))) +
  scale_color_brewer(palette = "Paired", 
                     labels = c(expression(italic("Clostridium_M bolteae")), 
                                expression(italic("Faecalicatena gnavus")),
                                expression(italic("Flavonifractor plautii")),
                                expression(italic("Flavonifractor sp000508885"))))



# tally genes -------------------------------------------------------------


cbolt_cd_up_gene %>%
  select(aa_seq) %>%
  distinct() %>%
  nrow() 
# 7004  

cbolt_cd_down_gene %>%
  select(aa_seq) %>%
  distinct() %>%
  nrow() 
# 981

rgnv_cd_up_gene %>%
  nrow()
# 3041

rgnv_cd_down_gene %>%
  nrow()
# 2943

fplaut_cd_down_gene %>%
  nrow()
# 4446

fplaut_cd_up_gene %>%
  nrow()
# 423 

fsp_cd_down_gene %>%
  nrow()
# 1959

fsp_cd_up_gene %>%
  nrow()
# 494


# dig in ------------------------------------------------------------------
cbolt_cd_up_gene_enriched@result %>%
  filter(Description == "Quorum sensing")
cbolt_cd_up_gene_enriched@result %>%
  filter(Description == "Vancomycin resistance")

rgnv_cd_up_gene_enriched@result %>%
  filter(Description == "Vancomycin resistance")

# write a parser to go from enriched pathway to KO definitions

parse_kegg_enriched <- function(enriched_result, kegg_definitions, pathway){
  enriched_kos <- enriched_result@result %>%
    filter(Description == pathway)
  enriched_kos <- read.delim(textConnection(enriched_kos$geneID), sep = "/",
                             header = F) %>%
    t() %>%
    as.data.frame() %>%
    mutate(ko = gsub("ko:", "", V1)) %>%
    select(-V1) %>%
    left_join(kegg_definitions, by = "ko") %>%
    distinct
  return(enriched_kos) 
}


parse_kegg_enriched(enriched_result = rgnv_cd_up_not_in_cd_down_enriched, 
                    kegg_definitions = kegg_def, 
                    pathway = "Vancomycin resistance")

parse_kegg_enriched(enriched_result = cbolt_cd_up_not_in_cd_down_enriched, 
                    kegg_definitions = kegg_def, 
                    pathway = "Vancomycin resistance")

parse_kegg_enriched(enriched_result = fplaut_cd_up_not_in_cd_down_enriched, 
                    kegg_definitions = kegg_def, 
                    pathway = "Flagellar assembly")

parse_kegg_enriched(enriched_result = fsp_cd_up_not_in_cd_down_enriched, 
                    kegg_definitions = kegg_def, 
                    pathway = "Flagellar assembly")

cbolt_cd_up_not_in_cd_down_annot <- cbolt_cd_up_not_in_cd_down %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  left_join(kegg_def, by = c("KEGG_ko" = "ko")) %>%
  distinct()
View(cbolt_cd_up_not_in_cd_down_annot)

rgnv_cd_up_not_in_cd_down_annot <- rgnv_cd_up_not_in_cd_down %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  left_join(kegg_def, by = c("KEGG_ko" = "ko")) %>%
  distinct()
View(rgnv_cd_up_not_in_cd_down_annot)


# oxidative stress --------------------------------------------------------
rnf <- c("K03612", "K03616", "K03615")
ox <- c("K04565", "K03386", "K03564", "K00362", "K01759", "K01758", "K14155", "K11717", "K03671",
         "K01926", "K01069", "K03676", "K03564")

ox1 <- rgnv_cd_up_not_in_cd_down_annot %>%
  filter(KEGG_ko %in% ox) %>%
  select(species, KEGG_ko, definition) %>%
  distinct()

ox2 <- cbolt_cd_up_not_in_cd_down_annot %>%
  filter(KEGG_ko %in% ox) %>%
  select(species, KEGG_ko, definition) %>%
  distinct()

ox3 <- fsp_cd_up_not_in_cd_down %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  left_join(kegg_def, by = c("KEGG_ko" = "ko")) %>%
  distinct() %>%
  filter(KEGG_ko %in% ox) %>%
  select(species, KEGG_ko, definition) %>%
  distinct()

ox4 <- fplaut_cd_up_not_in_cd_down %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  left_join(kegg_def, by = c("KEGG_ko" = "ko")) %>%
  distinct() %>%
  filter(KEGG_ko %in% ox) %>%
  select(species, KEGG_ko, definition) %>%
  distinct()

ox_all <- rbind(ox1, ox2, ox3, ox4)
write_delim(ox_all, "tmp_ox.tsv", delim = "|")


nitro <- c("K03385", "K00362", "K05601", "K05916", "K12264", "K12265",
           "K03386", "K24119", "K24126", "K04561")

nitro1 <- rgnv_cd_up_not_in_cd_down_annot %>%
  filter(KEGG_ko %in% nitro) %>%
  select(species, KEGG_ko, definition) %>%
  distinct()

nitro2 <- cbolt_cd_up_not_in_cd_down_annot %>%
  filter(KEGG_ko %in% nitro) %>%
  select(species, KEGG_ko, definition) %>%
  distinct()

nitro3 <- fsp_cd_up_not_in_cd_down %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  left_join(kegg_def, by = c("KEGG_ko" = "ko")) %>%
  distinct() %>%
  filter(KEGG_ko %in% nitro) %>%
  select(species, KEGG_ko, definition) %>%
  distinct()

nitro4 <- fplaut_cd_up_not_in_cd_down %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  left_join(kegg_def, by = c("KEGG_ko" = "ko")) %>%
  distinct() %>%
  filter(KEGG_ko %in% nitro) %>%
  select(species, KEGG_ko, definition) %>%
  distinct()

nitro_all <- rbind(nitro1, nitro2, nitro3, nitro4)
write_delim(nitro_all, "tmp_nitro.tsv", delim = "|")
# ulcerative colitis ------------------------------------------------------

### Flavronifractor plautii

fplaut_uc_up_gene <- pan_results %>%
  filter(mu == "mu.diagnosisUC") %>%
  filter(estimate > 0) %>%
  filter(genome == "ERS537353_12")

fplaut_uc_down_gene <- pan_results %>%
  filter(mu == "mu.diagnosisUC") %>%
  filter(estimate < 0) %>%
  filter(genome == "ERS537353_12")

fplaut_uc_up_not_in_uc_down <- fplaut_uc_up_gene %>%
  filter(!KEGG_ko %in% fplaut_uc_down_gene$KEGG_ko)
fplaut_uc_up_not_in_uc_down_enriched<- enricher(gene = fplaut_uc_up_not_in_uc_down$KEGG_ko, 
                                                TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(fplaut_uc_up_not_in_uc_down_enriched)
View(fplaut_uc_up_not_in_uc_down_enriched@result %>%
       filter(p.adjust < .05))


### Flavronifractor sp

fsp_uc_up_gene <- pan_results %>%
  filter(mu == "mu.diagnosisUC") %>%
  filter(estimate > 0) %>%
  filter(genome == "GCF_000508885.1_ASM50888v1_genomic")

fsp_uc_down_gene <- pan_results %>%
  filter(mu == "mu.diagnosisUC") %>%
  filter(estimate < 0) %>%
  filter(genome == "GCF_000508885.1_ASM50888v1_genomic")

fsp_uc_up_not_in_uc_down <- fsp_uc_up_gene %>%
  filter(!KEGG_ko %in% fsp_uc_down_gene$KEGG_ko)
fsp_uc_up_not_in_uc_down_enriched<- enricher(gene = fsp_uc_up_not_in_uc_down$KEGG_ko, 
                                             TERM2GENE = pathways, TERM2NAME = pathway_names)
clusterProfiler::dotplot(fsp_uc_up_not_in_uc_down_enriched)
View(fsp_uc_up_not_in_uc_down_enriched@result %>%
       filter(p.adjust < .05))

fsp_uc_up_gene %>%
  nrow()
# 37

fsp_uc_down_gene %>%
  nrow()
# 51

fplaut_uc_up_gene %>%
  nrow()
# 26

fplaut_uc_down_gene %>%
  nrow()
# 125