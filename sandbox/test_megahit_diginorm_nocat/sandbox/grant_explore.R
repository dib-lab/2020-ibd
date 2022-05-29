library(dplyr)
library(readr)
library(purrr)
library(tidyr)
eggnog_files <- list.files("sandbox/test_megahit_diginorm_nocat/eggnog", 
                           ".annotations$", full.names = T)

eggnog <- eggnog_files %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, .id = "genome", skip = 3) %>%
  mutate(genome = gsub(".emapper.annotations", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/eggnog/", "", genome)) %>%
  rename(query_name = "#query_name")


kegg_tally <- eggnog %>% 
  filter(!is.na(KEGG_ko)) %>%
  group_by(genome, KEGG_ko) %>%
  tally() %>%
  group_by(KEGG_ko) %>%
  tally()

ggplot(kegg_tally, aes(x = n)) +
  geom_histogram(bins = 13) +
  geom_text(stat='count', aes(label=..count..), vjust = -.5) +
  theme_minimal()


# limit to differentially abundant KOs ------------------------------------

# CD --------

sig_cd_files <- list.files("sandbox/test_megahit_diginorm_nocat/corncob", 
                           ".cd_sig_ccs_names.txt$", full.names = T)
sig_cd <- sig_cd_files %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, .id = "genome", col_names = "query_name") %>%
  mutate(genome = gsub("_cd_sig_ccs_names.txt", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/corncob/", "", genome)) 

sig_cd_eggnog <- left_join(sig_cd, eggnog, by = c('genome', 'query_name'))


kegg_tally_cd <- sig_cd_eggnog %>% 
  filter(!is.na(KEGG_ko)) %>%
  group_by(genome, KEGG_ko) %>%
  tally() %>%
  group_by(KEGG_ko) %>%
  tally()

ggplot(kegg_tally_cd, aes(x = n)) +
  geom_histogram(bins = 13) +
  geom_text(stat='count', aes(label=..count..), vjust = -.5) +
  theme_minimal()

K03205_cd <- sig_cd_eggnog %>%
  filter(KEGG_ko == "ko:K03205")

tmp <- K03205_cd %>%
  group_by(genome) %>%
  tally()


tmp2 <- eggnog %>%
  filter(KEGG_ko == "ko:K03205") %>%
  group_by(genome) %>%
  tally()


sig_ec_cd <- sig_cd_eggnog %>%
  filter(!is.na(EC))
View(sig_ec_cd)

length(unique(sig_ec_cd$EC))

franzosa_ec <- readxl::read_excel("~/Downloads/41564_2018_306_MOESM9_ESM.xlsx", skip = 1)
franzosa_ec<- separate(franzosa_ec, col = "# Feature / Statistic", into = c("EC", "description"), 
                       sep = ":", remove = T)
franzosa_ec_sig <- franzosa_ec %>%
  filter(`q-val_CD-Control` < .05 | `q-val_UC-Control` < .05)
View(franzosa_ec_sig)
table(franzosa_ec_sig$EC %in% sig_ec_cd$EC)

# UC --------

sig_uc_files <- list.files("sandbox/test_megahit_diginorm_nocat/corncob", 
                           ".uc_sig_ccs_names.txt$", full.names = T)
sig_uc <- sig_uc_files %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, .id = "genome", col_names = "query_name") %>%
  mutate(genome = gsub("_uc_sig_ccs_names.txt", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/corncob/", "", genome)) 

sig_uc_eggnog <- left_join(sig_uc, eggnog, by = c('genome', 'query_name'))

kegg_tally_uc <- sig_uc_eggnog %>% 
  filter(!is.na(KEGG_ko)) %>%
  group_by(genome, KEGG_ko) %>%
  tally() %>%
  group_by(KEGG_ko) %>%
  tally()

ggplot(kegg_tally_uc, aes(x = n)) +
  geom_histogram(bins = 13) +
  geom_text(stat='count', aes(label=..count..), vjust = -.5) +
  theme_minimal()

K03205_uc <- sig_uc_eggnog %>%
  filter(KEGG_ko == "ko:K03205")

tmp <- K03205_uc %>%
  group_by(genome) %>%
  tally()

sig_ec_uc <- sig_uc_eggnog %>%
  filter(!is.na(EC))

table(unique(sig_ec_uc$EC) %in% unique(sig_ec_cd$EC))
table(franzosa_ec_sig$EC %in% c(sig_ec_uc$EC, sig_ec_cd$EC))
