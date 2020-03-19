setwd("~/github/2020-ibd")

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
library(ggrepel)

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
genbank_names <- read_tsv("genbank_gather_matches_full_genomes.tsv") %>%
  select(accession, id2)


# read in gather res ------------------------------------------------------

# read in gather results for the original genomes
gather_db <- read_csv("sandbox/gather_vita_hashes/rf_vita_hashes.csv") 

# edit the gather accession names so they'll match the gtdb accessions
# gsub everything after first space to get genbank accession
# gsub out other database paths
gather_db$accession <- gather_db$name
gather_db$accession <- gsub("\\s.*", "", gather_db$accession)
gather_db$accession <- gsub("HGM_v1\\.0_all_60664_fna\\/", "", gather_db$accession)
gather_db$accession <- gsub("\\.\\/[0-9]{4,5}\\/", "", gather_db$accession)
gather_db$accession <- gsub("\\.[12]$", "", gather_db$accession)
gather_db$accession <- gsub("\\.fa\\.gz$", "", gather_db$accession)
gather_db$accession <- gsub("\\.fa$", "", gather_db$accession)
gather_db$accession <- gsub("\\.fna$", "", gather_db$accession)

# join to the genbank names
gather_db <- left_join(gather_db, genbank_names, by = "accession")
# replace accessions that were different between gather and genbank
# (i.e. some gather results have an accession that matches the first
# sequence in the fasta. This seq may be a plasmid, not the actual
# genome that a signature was calculated for. I hand identified the
# real genome accessions for these orgs.)
gather_db$accession <- ifelse(is.na(gather_db$id2), gather_db$accession, gather_db$id2)

# filter db results to relevant columns
gather_db <- gather_db %>% 
  mutate(f_unique_weighted_db = f_unique_weighted) %>%
  mutate(id = accession) %>%
  mutate(f_orig_query_db = f_orig_query) %>%
  select(id, f_orig_query_db, f_unique_weighted_db)

# filter sgc results to relevant columns
gather_sgc <- read_csv("sandbox/vita_var_vs_sgc_merged.csv") %>%
  mutate(f_unique_weighted_sgc = f_unique_weighted) %>%
  mutate(id = name) %>%
  mutate(f_orig_query_sgc = f_orig_query) %>%
  select(id, f_orig_query_sgc, f_unique_weighted_sgc)

# join sgc and db results
gather_all <- full_join(gather_db, gather_sgc, by = "id")
# join gather results wtih the gtdb taxonomy assignments
res <- left_join(gather_all, gtdb, by = c("id" = "accession"))

# calculate fold change
res$fc_f_orig_query <- res$f_orig_query_sgc / res$f_orig_query_db
res$fc_f_unique_weighted <- res$f_unique_weighted_sgc / res$f_unique_weighted_db

# my data -----------------------------------------------------------------

# read the phylogenetic tree produced by gtdbtk
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

# visualize the tree 
p <- ggtree(tree) 

# add fold change bar plots 
p <- p + geom_facet(panel = "Fold change in predictive k-mers, all", 
                    data = res, geom = ggstance::geom_barh, 
                    aes(x = as.numeric(fc_f_orig_query)), 
                    stat = "identity", width = .6) +
  geom_facet(panel = "Fold change in predictive k-mers, unique", 
             data = res, geom = ggstance::geom_barh, 
             aes(x = as.numeric(fc_f_unique_weighted)), 
             stat = "identity", width = .6) +
  theme_tree2(legend.position=c(.05, .85))

# label the tips
res$name <- ifelse(res$species == "", res$genus, res$species)
p$data$label <- ifelse(p$data$label %in% res$id, res$name, p$data$label)
p <- p + geom_tiplab(size = 2)

# write pdf 
pdf("figures/gather_species_foldchange.pdf", height = 10, width = 35)
p
dev.off()

# by order ----------------------------------------------------------------
# labels the tree by order, instead of the tips having species names
# looks cleaner, but probably less information bc theirs no 
# order-level pattern.

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

pdf("figures/gather_order_viz_foldchange.pdf", width = 10, height = 7)
p
dev.off()

res$name <- ifelse(res$species == "", res$genus, res$species)
p$data$label <- ifelse(p$data$label %in% res$id, res$name, p$data$label)

# plot scatter plot -------------------------------------------------------

head(res)
colnames(res)

ggplot(res, aes(x = f_orig_query_db, y = f_orig_query_sgc, color = order)) +
  geom_point(alpha = .5) + 
  xlim(0, .11) + ylim(0, .11) +
  theme_minimal() +
  labs(x = "Fraction in databases", y = "Fraction in databases + \npangenome nbhd",
       title = "Fraction of predictive hashes")

ggplot(res, aes(x = f_orig_query_db, y = f_orig_query_sgc)) +
  geom_point(alpha = .5) + 
  xlim(0, .11) + ylim(0, .11) +
  theme_minimal() +
  labs(x = "Fraction in databases", y = "Fraction in databases + \npangenome nbhd",
       title = "Fraction of predictive hashes") +
  geom_text_repel(data=subset(res, f_orig_query_sgc > .092),
                  aes(f_orig_query_db, y = f_orig_query_sgc, 
                      label= ifelse(species == "", genus, species)),
                  nudge_x = .05)


ggplot(res, aes(x = f_orig_query_db, y = f_orig_query_sgc)) +
  geom_point(alpha = .5) + 
  xlim(0, .11) + ylim(0, .11) +
  theme_minimal() +
  labs(x = "Fraction in databases", y = "Fraction in databases + \npangenome nbhd",
       title = "Fraction of predictive hashes") +
  geom_text_repel(data=subset(res, f_orig_query_sgc > .092),
                  aes(f_orig_query_db, y = f_orig_query_sgc, 
                      label= ifelse(species == "", genus, species)),
                  nudge_x = .05) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey")

ggplot(res, aes(x = f_unique_weighted_db, y = f_unique_weighted_sgc)) +
  geom_point(alpha = .5) + 
  xlim(0, .11) + ylim(0, .11) +
  theme_minimal() +
  labs(x = "Unique fraction in databases", y = "Unique fraction in databases + \npangenome nbhd",
       title = "Fraction of predictive hashes") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  geom_text_repel(data=subset(res, f_unique_weighted_sgc > .04),
                  aes(f_unique_weighted_db, y = f_unique_weighted_sgc, 
                      label= ifelse(species == "", genus, species)),
                  nudge_x = .05) 

ggplot(res, aes(x = f_unique_weighted_db, y = f_unique_weighted_sgc)) +
  geom_point(alpha = .5) + 
  xlim(0, .11) + ylim(0, .11) +
  theme_minimal() +
  labs(x = "Unique fraction in databases", y = "Unique fraction in databases + \npangenome nbhd",
       title = "Fraction of predictive hashes") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey") +
  geom_text_repel(data=subset(res, f_unique_weighted_sgc > .04),
                  aes(f_unique_weighted_db, y = f_unique_weighted_sgc, 
                      label= ifelse(species == "", genus, species)),
                  nudge_x = .05)
