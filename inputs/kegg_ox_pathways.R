library(readr)
setwd("~/github/2020-ibd")

# Rnf complex is a respiratory enzyme that catalyzes the oxidation of reduced
# ferredoxin to the reduction of NAD+, and the negative free energy change of 
# this reaction is used to generate a transmembrane ion gradient 
# (https://journals.asm.org/doi/full/10.1128/JB.00357-18).
rnf_ko <- c("K03612", "K03616", "K03615")

# oxidative stress kos
ox_ko <- c('K00362', 'K11717', 'K01926', 'K01069', 'K03671', 'K04565', 'K14155', 
           'K11717', 'K01069', 'K14155', 'K03386', 'K01758', 'K03564', 'K00362', 
           'K01926', 'K04565', 'K03676', 'K01759', 'K03671', 'K01758', 'K01759', 
           'K01758', 'K14155')
# 
nox_ko <- c("K03385", "K00362", "K05601", "K05916", "K12264", "K12265",
         "K03386", "K24119", "K24126", "K04561")


# kegg files ------------------------------------------------------------
# kegg_def <- read_csv("inputs/kegg_ortholog_definitions.csv")
# pathways <- read_delim("http://rest.kegg.jp/link/pathway/ko", delim = "\t", col_names = c("KO", "path"))
# pathways <- pathways %>%
#   #mutate(KO = gsub("ko:", "", KO)) %>%
#   filter(grepl(pattern = "map", path)) %>%
#   #mutate(path = gsub("path:", "", path)) %>%
#   dplyr::select(path, KO)
# pathway_names <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))