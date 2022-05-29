library(ggtree)



tree <- read.tree("sandbox/test_megahit_diginorm_nocat/sandbox/try_strainplan/output_lax/RAxML_bestTree.s__Ruminococcus_gnavus.StrainPhlAn3.tre")
tree$tip.label <- gsub("_GCF_900036035.1_RGNV35913_genomic.fna.cdgb_ids.reads.diginorm", "", tree$tip.label)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, diagnosis, study_accession) %>%
  distinct() %>%
  filter(library_name %in% tree$tip.label)

df <- data.frame(library_name = tree$tip.label) %>%
  left_join(info) %>%
  mutate(diagnosis = ifelse(is.na(diagnosis), "RefSeq", diagnosis))
View(df)

ggtree(tree) %<+% df + 
  geom_tiplab(size = 1) 
  #geom_tippoint(aes(color = diagnosis)) +
  #scale_color_manual(values = c(CD = "darkorchid", nonIBD = "cadetblue", UC = "pink", 
  #                              RefSeq = "grey"))
  
