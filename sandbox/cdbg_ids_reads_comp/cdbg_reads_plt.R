setwd("~/github/2020-ibd/")
library(readr)
library(Rtsne)
library(dplyr)
library(ggplot2)

comp <- read.csv("sandbox/cdbg_ids_reads_comp/GCF_000371685.1_Clos_bolt_90B3_V1_genomic.cdbg_ids.reads.comp.csv")
colnames(comp) <- gsub("_GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna", "", colnames(comp))

comp <- read.csv("sandbox/cdbg_ids_reads_comp/SRS1719577_6.cdbg_ids.reads.comp.csv")
colnames(comp) <- gsub("_SRS1719577_6.fna", "", colnames(comp))

comp <- read.csv("sandbox/cdbg_ids_reads_comp/SRS1719498_9.cdbg_ids.reads.comp.csv")
colnames(comp) <- gsub("_SRS1719498_9.fna", "", colnames(comp))

rownames(comp) <- colnames(comp)

# read in metadata
info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, diagnosis) %>%
  distinct() %>%
  mutate(sample =  gsub("-", "\\.", library_name))

# Transform for plotting
comp <- as.matrix(comp)

# Make an MDS plot

fit <- dist(comp)
fit <- cmdscale(fit)
x <- fit[, 1]
y <- fit[, 2]
tmp <- data.frame(x = x, y = y, sample = colnames(comp))
tmp$sample <- gsub("^X", "", tmp$sample)
tmp <- left_join(tmp, info, by = "sample")

ggplot(tmp, aes(x = x, y = y, color = diagnosis)) +
  geom_point() +
  theme_minimal() +
  scale_color_viridis_d()

tsne_model <- Rtsne(comp, check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
d_tsne <- as.data.frame(tsne_model$Y) 
d_tsne <- data.frame(x = d_tsne$V1, y = d_tsne$V2, sample = colnames(comp))
d_tsne$sample <- gsub("^X", "", d_tsne$sample)
d_tsne <- left_join(d_tsne, info, by = "sample")
plt <- ggplot(d_tsne, aes(x = x, y = y, color = diagnosis)) +
  geom_point() +
  theme_minimal() + 
  theme(plot.title.position = "plot") +
  scale_color_viridis_d() +
  labs(x = "", y = "", title = "SRS1719498_9")

ggExtra::ggMarginal(plt, type = "density",
                    groupColour = T, groupFill = T)
