setwd("~/github/ibd")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ggrepel)
library(stringr)
library(plotly)
library(viridis)

## read in metadata
info <- read_tsv("inputs/working_metadata.tsv") 

## sig names were not lib names in comp mat, so read in 1:1 lib name:sig name
## to replace in comp mat
sig_names <- read_tsv("sandbox/greater_than_one_filt_sigs/sig_names.txt", 
                      col_names = F)
sig_names <- sig_names %>%
  mutate(ind = rep(c(1, 2),length.out = n())) %>%
  group_by(ind) %>%
  mutate(id = row_number()) %>%
  spread(ind, X1) %>%
  select(-id) 
colnames(sig_names) <- c("filename", "md5")
sig_names$filename <- gsub("signature filename: ", "", sig_names$filename) 
sig_names$filename <- gsub("_filt.sig", "", sig_names$filename) 
sig_names$md5 <- gsub("signature: ", "", sig_names$md5)

## Read in comp mat
comp <- read_csv("sandbox/filt_comp/filt_comp_all.csv") # read in mat
## order sig_names by comp colnames
sig_names <- sig_names[order(match(sig_names$md5, colnames(comp))), ]
## remove the last row
sig_names <- sig_names[1:954, ]

colnames(comp) <- sig_names$filename
rownames(comp) <- colnames(comp)                       # Label the rows

dist <- dist(comp)                                     # compute distances
fit_all <- cmdscale(dist, eig = T)                     # calculate MDS
fit <- as.data.frame(fit_all$points)
fit$sample <- rownames(fit)

colnames(fit) <- c("dim1", "dim2", "sample")       # set column names
fit <- left_join(fit, info, by = c("sample" = "library_name")) # join with metadata

var <- round(fit_all$eig*100/sum(fit_all$eig), 1)

comp_plt <- ggplot(fit, aes(x = dim1, y = dim2, label = read_count,
                            color = study_accession, shape = diagnosis)) +
  geom_point() +
  theme_minimal() +
  ggtitle("IBD studies") +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  scale_color_viridis_d()

ggExtra::ggMarginal(comp_plt, type = "density",
                    groupColour = T, groupFill = T)
comp_plt
ggplotly(comp_plt) # interactive html plot


comp_plt <- ggplot(fit, aes(x = dim1, y = dim2, label = read_count,
                            color = diagnosis, shape = study_accession)) +
  geom_point() +
  theme_minimal() +
  ggtitle("IBD studies") +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  scale_color_viridis_d()

ggExtra::ggMarginal(comp_plt, type = "density",
                    groupColour = T, groupFill = T)
comp_plt

