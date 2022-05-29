setwd("~/github/ibd")

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ggthemes)
library(devtools)
#devtools::install_github("GuangchuangYu/ggtree", force = T)
library(ggtree)
library(dplyr)
library(tidytree) 
library(ggnewscale)
library(ggstance)
library(ggtree)
library(cowplot)
library(treeio)

# read in and format taxonomy ---------------------------------------------

# Read in NCBI and GTDB taxonomies. Use GTDB taxonomy unless none is recorded,
# then use NCBI taxonomy. 

# gtdb contains 128 records; one was recorded as archaea due to contamination
gtdb <- read_tsv("sandbox/gather_vita_hashes/gtdbtk.bac120.summary.tsv") %>%
  mutate(accession = user_genome) %>%
  select(accession, classification) %>%
  mutate(classification = gsub("[dkpcofgs]__", "", classification)) %>%
  separate(data = ., col = classification, sep = ";",
           into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
  mutate(strain = NA) %>%
  mutate(source = "gtdbtk") 

# some sourmash gather results annotated to a single chromosome. 
# genbank_names records the mapping between the gather name and the
# gtdb accession name.
genbank_names<- read_tsv("genbank_gather_matches_full_genomes.tsv") %>%
  select(accession, id2)


# read in gather res ------------------------------------------------------

# read in gather results, which contain abundance information
gather <- read_csv("sandbox/vita_var_vs_sgc_merged.csv")
# edit the gather accession names so they'll match the gtdb accessions
# gsub everything after first space to get genbank accession
# gsub out other database paths
gather$accession <- gather$name
gather$accession <- gsub("\\s.*", "", gather$accession)
gather$accession <- gsub("HGM_v1\\.0_all_60664_fna\\/", "", gather$accession)
gather$accession <- gsub("\\.\\/[0-9]{4,5}\\/", "", gather$accession)
gather$accession <- gsub("\\.[12]$", "", gather$accession)
gather$accession <- gsub("\\.fa\\.gz$", "", gather$accession)
gather$accession <- gsub("\\.fa$", "", gather$accession)
gather$accession <- gsub("\\.fna$", "", gather$accession)

# join gather results wtih the gtdb taxonomy assignments
gather <- left_join(gather, genbank_names, by = "accession")
gather$accession <- ifelse(is.na(gather$id2), gather$accession, gather$id2)
res <- left_join(gather, gtdb, by = "accession")

res <- res %>% 
  filter(superkingdom == "Bacteria") %>%
  mutate(id = accession) %>%
  select(id, f_match, f_unique_weighted, superkingdom, phylum, class, 
         order, family, genus, species)



# my data -----------------------------------------------------------------

## read the phylogenetic tree produced by gtdbtk
tree <- read.newick("sandbox/gather_vita_hashes/gtdbtk.bac120.classify.tree")


# drop tips from GTDB. This tree contains all GTDB species. 
# drop these so we only visualize 129 genomes.
drop_tips <- tree$tip.label[grepl(pattern = "^GB_", tree$tip.label)]
tree <- drop.tip(tree, drop_tips)
drop_tips <- tree$tip.label[grepl(pattern = "^RS_", tree$tip.label)]
tree <- drop.tip(tree, drop_tips)

# order res by tip labels -- otherwise, tip labels will be pulled as 
# order or results, which will not accurately map the info.
res <- res[match(tree$tip.label, res$id), ]
# # make rownames tip labels
# rownames(res) <- res$id

## visualize the tree 
p <- ggtree(tree) 

p <- p + geom_facet(panel = "Fraction of predictive k-mers", 
                    data = res, geom = ggstance::geom_barh, 
                    aes(x = as.numeric(f_unique_weighted)), 
                    stat = "identity", width = .6) +
  geom_facet(panel = "Fraction of pangenome", 
             data = res, geom = ggstance::geom_barh, 
             aes(x = as.numeric(f_match)), 
             stat = "identity", width = .6) +
  theme_tree2(legend.position=c(.05, .85))

res$name <- ifelse(res$species == "", res$genus, res$species)
p$data$label <- ifelse(p$data$label %in% res$id, res$name, p$data$label)

p <- p + geom_tiplab(size = 2)

pdf("figures/gather_species_viz_sgc.pdf", height = 10, width = 35)
p
dev.off()

# by order ----------------------------------------------------------------

## visualize the tree 
p <- ggtree(tree) 
p$data$label <- ifelse(p$data$label %in% res$id, res$species, p$data$label)

p + 
  geom_cladelabel(node = 137, fontsize = 3, angle =-20, label = expression(italic("Oscillospirales"))) +
  geom_cladelabel(node = 222, fontsize = 3, angle =-20, label = expression(italic("Lactobacillales"))) +
  geom_cladelabel(node = 226, fontsize = 3, angle =-20, label = expression(italic("Erysipelotrichales"))) +
  geom_cladelabel(node = 234, fontsize = 3, angle =-20, label = expression(italic("Bacteroidales"))) +
  geom_cladelabel(node = 171, fontsize = 3, angle =-20, label = expression(italic("Lachnospirales")))
# geom_tiplab(color = "gray", size = 1)


p <- ggtree(tree) 
p <- p + 
  geom_cladelabel(node = 137, offset.text= .1, fontsize = 3, hjust = .5, angle =90, label = expression(italic("Oscillospirales"))) +
  geom_cladelabel(node = 222, offset.text= .1, fontsize = 3, hjust = .5, angle =90, label = expression(italic("Lactobacillales"))) +
  geom_cladelabel(node = 226, offset.text= .1, fontsize = 3, hjust = .5, angle =90, label = expression(italic("Erysipelotrichales"))) +
  geom_cladelabel(node = 234, offset.text= .1, fontsize = 3, hjust = .5, angle =90, label = expression(italic("Bacteroidales"))) +
  geom_cladelabel(node = 171, offset.text= .1, fontsize = 3, hjust = .5, angle =90, label = expression(italic("Lachnospirales")))
p <- p + geom_facet(panel = "Fraction of predictive k-mers", 
                    data = res, geom = ggstance::geom_barh, 
                    aes(x = as.numeric(f_unique_weighted)), 
                    stat = "identity", width = .6) +
  geom_facet(panel = "Fraction of pangenome", 
             data = res, geom = ggstance::geom_barh, 
             aes(x = as.numeric(f_match)), 
             stat = "identity", width = .6) +
  theme_tree2(legend.position=c(.05, .85))

pdf("figures/gather_order_viz_sgc.pdf", width = 10, height = 7)
p
dev.off()

res$name <- ifelse(res$species == "", res$genus, res$species)
p$data$label <- ifelse(p$data$label %in% res$id, res$name, p$data$label)


# check res with gtotree --------------------------------------------------

# mapping <- gtdb %>% 
#   select(accession, genus, species) %>%
#   mutate(species = ifelse(species == "", genus, species)) %>%
#   select(accession, species) %>%
#   mutate(accession = paste0(accession, ".fna"))
# 
# write_tsv(mapping, "~/Downloads/gather_genomes/genome_to_id_map.tsv")

