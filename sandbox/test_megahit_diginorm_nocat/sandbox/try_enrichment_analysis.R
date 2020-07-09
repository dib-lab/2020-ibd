library(dplyr)
library(readr)
library(purrr)
library(clusterProfiler)


# read in files -----------------------------------------------------------
info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct

eggnog_files <- list.files("sandbox/test_megahit_diginorm_nocat/eggnog", 
                           ".annotations$", full.names = T)

eggnog <- eggnog_files %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, .id = "genome", skip = 3) %>%
  mutate(genome = gsub(".emapper.annotations", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/eggnog/", "", genome)) %>%
  rename(query_name = "#query_name")

corncob_files <- list.files("sandbox/test_megahit_diginorm_nocat/corncob",
                            ".sig_ccs.tsv$", full.names = T) 
corncob <- corncob_files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genome") %>%
  mutate(genome = gsub("_sig_ccs.tsv", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/corncob/", "", genome))
View(corncob)  

corncob <- left_join(corncob, eggnog, by = c("genome", "aa_seq" = "query_name"))

species <- read_tsv("~/Desktop/41pangenomes_gather_matches_ranked.tsv") %>%
  select(genome, GTDB, NCBI, rank)
corncob <- left_join(corncob, species, by = "genome")


# get kegg universe -------------------------------------------------------

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


# all by disease and direction --------------------------------------------
corncob_up_cd <- corncob %>%
  filter(estimate > 0) %>%
  filter(mu == "mu.diagnosisCD")
up_cd_en <- enricher(gene = corncob_up_cd$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(up_cd_en)

corncob_down_cd <- corncob %>%
  filter(estimate < 0)%>%
  filter(mu == "mu.diagnosisCD")
down_cd_en <- enricher(gene = corncob_down_cd$KEGG_ko, 
                     TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(down_cd_en)

tmp_down <- corncob_down_cd$KEGG_ko[!is.na(corncob_down_cd$KEGG_ko)]
tmp_up <- corncob_up_cd$KEGG_ko[!is.na(corncob_up_cd$KEGG_ko)]
table(tmp_down %in% tmp_up)
# many of the orthologs overlap. This is probably a problem if they overlap *within* a species. 
corncob_up_uc <- corncob %>%
  filter(estimate > 0) %>%
  filter(mu == "mu.diagnosisUC")

up_uc_en <- enricher(gene = corncob_up_uc$KEGG_ko, 
                       TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(up_uc_en)

corncob_down_uc <- corncob %>%
  filter(estimate < 0) %>%
  filter(mu == "mu.diagnosisUC")

down_uc_en <- enricher(gene = corncob_down_uc$KEGG_ko, 
                             TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(down_uc_en)


# try with top-ranked genome ----------------------------------------------

# SRS1719577_6.fna

rank2 <- corncob %>%
  filter(genome == "SRS1719577_6.fna")

rank2_up_uc <- rank2 %>%
  filter(estimate > 0 ) %>%
  filter(mu == "mu.diagnosisUC")
rank2_down_uc <- rank2 %>%
  filter(estimate < 0) %>%
  filter(mu == "mu.diagnosisUC")
table(is.na(rank2_down_uc$KEGG_ko))
rank2_down_uc_en <- enricher(gene = rank2_down_uc$KEGG_ko, 
                              TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank2_down_uc_en)


rank2_up_cd <- rank2 %>%
  filter(estimate > 0 ) %>%
  filter(mu == "mu.diagnosisCD")
rank2_up_cd_en <- enricher(gene = rank2_up_cd$KEGG_ko, 
                            TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank2_up_cd_en)
View(rank2_up_cd_en@result)

rank2_down_cd <- rank2 %>%
  filter(estimate < 0) %>%
  filter(mu == "mu.diagnosisCD")

rank2_down_cd_en <- enricher(gene = rank2_down_cd$KEGG_ko, 
                              TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank2_down_cd_en)
View(rank2_down_cd_en@result)

# try with r gnavus, which has a known mechanism --------

rank13 <- corncob %>%
  filter(genome == "GCF_900036035.1_RGNV35913_genomic.fna")

rank13_up_uc <- rank13 %>%
  filter(estimate > 0 ) %>%
  filter(mu == "mu.diagnosisUC")
rank13_down_uc <- rank13 %>%
  filter(estimate < 0) %>%
  filter(mu == "mu.diagnosisUC")
table(is.na(rank13_down_uc$KEGG_ko))
rank13_down_uc_en <- enricher(gene = rank13_down_uc$KEGG_ko, 
                         TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank13_down_uc_en)


rank13_up_cd <- rank13 %>%
  filter(estimate > 0 ) %>%
  filter(mu == "mu.diagnosisCD")
rank13_up_cd_en <- enricher(gene = rank13_up_cd$KEGG_ko, 
                              TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank13_up_cd_en)
View(rank13_up_cd_en@result)

rank13_down_cd <- rank13 %>%
  filter(estimate < 0) %>%
  filter(mu == "mu.diagnosisCD")

rank13_down_cd_en <- enricher(gene = rank13_down_cd$KEGG_ko, 
                              TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank13_down_cd_en)
View(rank13_down_cd_en@result)

rank13_down_cd$KEGG_ko %in% rank13_up_cd$KEGG_ko
tmp_down <- rank13_down_cd$KEGG_ko[!is.na(rank13_down_cd$KEGG_ko)]
tmp_up <- rank13_up_cd$KEGG_ko[!is.na(rank13_up_cd$KEGG_ko)]
table(unique(tmp_down) %in% unique(tmp_up))
table(unique(tmp_up) %in% unique(tmp_down))

tmp_up_only<- tmp_up[!tmp_up %in% tmp_down]
# select only first KO
tmp_up_only <- gsub(",.*", "", tmp_up_only)
rank13_up_cd_only_en <- enricher(gene = tmp_up_only, 
                              TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank13_up_cd_only_en)
View(rank13_up_cd_only_en@result)
rank13_up_cd_only_en@result %>%
  filter(Description == "Galactose metabolism")

# look at R gnavus specific mechanism of action

read_blast_csv <- function(path){
  blast <- read_csv(path,
                    col_names = c("qseqid", "sseqid", "pident", "length",
                                  "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                                  "evalue", "bitscore"))
}
blast <- read_blast_csv("~/Downloads/G7T9G12Z114-Alignment-HitTable.csv")
# select first blast match for each protein
best_hits <- blast %>%
  group_by(qseqid) %>%
  top_n(1) %>%
  filter(length > 100) %>%
  filter(pident > 60)

View(rank13_up_cd[rank13_up_cd$aa_seq %in% best_hits$sseqid, ])
rank13_down_cd[rank13_down_cd$aa_seq %in% best_hits$sseqid, ]

rank13_up_uc[rank13_up_uc$aa_seq %in% best_hits$sseqid, ]
rank13_down_uc[rank13_down_uc$aa_seq %in% best_hits$sseqid, ]

# look at r gnavus complete query nbhd assemblies
# bowers et al: (>90% complete, <5% contamination), medium-quality draft (>50% complete, <10% contamination) 
checkm <- read_tsv("sandbox/test_megahit_diginorm_nocat/sandbox/checkm_GCF_900036035.1_RGNV35913_genomic.fna/completeness.tsv") %>%
  filter(Completeness > 50) %>%
  filter(Contamination < 10) %>%
  mutate(library_name = gsub("_.*", "", `Bin Id`)) %>%
  left_join(info, by = "library_name")
View(checkm)
table(checkm$diagnosis)
write.table(checkm$library_name, 
            "sandbox/test_megahit_diginorm_nocat/sandbox/checkm_GCF_900036035.1_RGNV35913_genomic.fna/quality_mags.txt",
            quote = F, row.names = F)

## see if biosynthetic cluster appears in only CD metagenomes
rg_counts <- read_tsv("sandbox/tximport/GCF_900036035.1_RGNV35913_genomic.fna_counts_raw.tsv")
rg_biosynth_counts <- rg_counts %>%
  filter(protein %in% best_hits$sseqid)
tmp <- rg_biosynth_counts %>%
  select(-protein) %>%
  colSums() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "library_name") %>%
  mutate(mean_biosynth = ./18) %>%
  left_join(info, by = "library_name")

ggplot(tmp, aes(x = diagnosis, y = mean_biosynth)) +
  geom_violin() +
  scale_y_log10() +
  theme_minimal()

# c bolt ------------------------------------------------------------------

cbolt_cd <- corncob %>%
  filter(genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna") %>%
  filter(mu == "mu.diagnosisCD")

lozupone <- read_csv("sandbox/test_megahit_diginorm_nocat/sandbox/doi_10.1101_gr.138198.12_S4.csv")
cbolt_cd$KEGG2 <- gsub("ko:", "", cbolt_cd$KEGG_ko)

table(lozupone$KEGG %in% cbolt_cd$KEGG2)
View(cbolt_cd[cbolt_cd$KEGG2 %in% lozupone$KEGG, ])
View(lozupone[lozupone$KEGG %in% cbolt_cd$KEGG2, ])

cbolt_up_cd <- cbolt_cd %>%
  filter(estimate > 0 ) 
cbolt_up_cd_en <- enricher(gene = cbolt_up_cd$KEGG_ko, 
                            TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(cbolt_up_cd_en)
View(cbolt_up_cd_en@result)
cbolt_up_cd_en@result %>%
  filter(Description == 'Flagellar assembly')

cbolt_down_cd <- cbolt_cd %>%
  filter(estimate < 0)
cbolt_down_cd_en <- enricher(gene = cbolt_down_cd$KEGG_ko, 
                           TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(cbolt_down_cd_en)
View(cbolt_down_cd_en@result)


cbolt_uc <- corncob %>%
  filter(genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna") %>%
  filter(mu == "mu.diagnosisUC")


# rtim --------------------------------------------------------------------

rank9 <-  corncob %>%
  filter(genome == "SRS294916_20.fna")

rank9_up_uc <- rank9 %>%
  filter(estimate > 0 ) %>%
  filter(mu == "mu.diagnosisUC")
rank9_down_uc <- rank9 %>%
  filter(estimate < 0) %>%
  filter(mu == "mu.diagnosisUC")
table(is.na(rank9_down_uc$KEGG_ko))
rank9_down_uc_en <- enricher(gene = rank9_down_uc$KEGG_ko, 
                             TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank9_down_uc_en)
View(rank9_down_uc_en@result)
View(rank9_down_uc)

tmp <- enricher(gene = c("ko:K02032", "ko:K02931","ko:K08303", "ko:K07814", "ko:K02074"),
         TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(tmp)


# actulibacter ---------------------------------------------------------------------

rank19 <-  corncob %>%
  filter(genome == "XieH_2016__YSZC12003_37172__bin.63.fa")

rank19_up_uc <- rank19 %>%
  filter(estimate > 0 ) %>%
  filter(mu == "mu.diagnosisUC")
rank19_down_uc <- rank19 %>%
  filter(estimate < 0) %>%
  filter(mu == "mu.diagnosisUC")

rank19_up_uc_en <- enricher(gene = rank19_up_uc$KEGG_ko, 
                             TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank19_up_uc_en)
View(rank19_down_uc_en@result)
View(rank19_down_uc)

# anan --------



rank25 <-  corncob %>%
  filter(genome == "LoombaR_2017__SID1050_bax__bin.11.fa")

rank25_up_uc <- rank25 %>%
  filter(estimate > 0 ) %>%
  filter(mu == "mu.diagnosisUC")
rank25_down_uc <- rank25 %>%
  filter(estimate < 0) %>%
  filter(mu == "mu.diagnosisUC")

rank25_down_uc_en <- enricher(gene = rank25_down_uc$KEGG_ko, 
                              TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank25_down_uc_en)
View(rank25_down_uc_en@result)
View(rank25_down_uc)

rank25_up_uc_en <- enricher(gene = rank25_up_uc$KEGG_ko, 
                              TERM2GENE = pathways, TERM2NAME = pathway_names)
dotplot(rank25_up_uc_en)

