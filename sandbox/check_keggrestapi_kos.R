library(KEGGREST)

pathway_names2 <- read_tsv("http://rest.kegg.jp/list/pathway", col_names = c("path", "name"))


pathways_names <- keggList("pathway")
pathways_names <- pathways %>% 
  as.data.frame() %>%
  rownames_to_column("pathway")

pathways_names <- keggList("pathway")
pathways_names <- pathways %>% 
  as.data.frame() %>%
  rownames_to_column("pathway")

# link all pathways to their kos
kegg_map_pathways <- character()
for(i in 1:nrow(pathway_names)){
  kegg_map_pathway <- keggLink("ko", pathway_names[i, ]$path)
  kegg_map_pathways <- c(kegg_map_pathway, kegg_map_pathways)
  Sys.sleep(3)
}

pathways <- read_tsv("http://rest.kegg.jp/link/pathway/ko", col_names = c("KO", "path"))
pathways <- pathways %>%
  #mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "map", path)) %>%
  #mutate(path = gsub("path:", "", path)) %>%
  dplyr::select(path, KO)
