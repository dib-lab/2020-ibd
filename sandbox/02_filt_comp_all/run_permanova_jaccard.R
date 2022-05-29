setwd("~/github/ibd")
library(vegan)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ggrepel)
library(stringr)
library(plotly)
library(viridis)


# read in comp and format  ----------------------------------------------
comp <- read_csv("sandbox/02_filt_comp_all/comp_all_jaccard.csv") # read in mat
colnames(comp) <- gsub("_filt", "", colnames(comp))
rownames(comp) <-colnames(comp) 


# read in and formate metadata -------------------------------------------

info <- read_tsv("inputs/working_metadata.tsv") %>%
  filter(library_name %in% colnames(comp)) %>%
  group_by(library_name, study_accession, diagnosis, subject) %>%
  summarise(read_count = sum(read_count)) %>%
  as.data.frame()

info_hmp <- read_tsv("inputs/hmp2_mgx_metadata.tsv") %>%
  mutate(study_accession = "hmp") %>%
  select(External.ID, study_accession, diagnosis, Participant.ID, reads_raw) %>%
  distinct() %>%
  filter(External.ID %in% colnames(comp))
colnames(info_hmp) <- c("library_name", "study_accession", "diagnosis", "subject", "read_count")

info <- rbind(info, info_hmp)


# dist and permanova ------------------------------------------------------

dist <- dist(comp)                                     # compute distances
# sort info by colnames
info <- info[match(colnames(comp), info$library_name), ]
info$diagnosis <- as.factor(info$diagnosis)
info$study_accession <- as.factor(info$study_accession)
info$subject <- as.factor(info$subject)

perm <- adonis(dist ~ diagnosis + study_accession + read_count + subject, 
               data = info, 
               permutations = 1000)


perm
# Call:
#   adonis(formula = dist ~ diagnosis + study_accession + read_count +  subject, data = info, permutations = 1000) 
# 
# Permutation: free
# Number of permutations: 1000
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
#   diagnosis          2     319.7 159.832 128.348 0.03464 0.000999 ***
#   study_accession    5    1124.9 224.972 180.657 0.12191 0.000999 ***
#   read_count         1      62.0  62.029  49.810 0.00672 0.000999 ***
#   subject          689    5735.8   8.325   6.685 0.62161 0.000999 ***
#   Residuals       1594    1985.0   1.245         0.21512             
# Total           2291    9227.3                 1.00000             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# plot --------------------------------------------------------------------

fit_jacc_all <- cmdscale(dist, eig = T)                             # calculate MDS
fit_jacc <- as.data.frame(fit_jacc_all$points)
fit_jacc$sample <- rownames(fit_jacc)

colnames(fit_jacc) <- c("dim1", "dim2", "sample")                   # set column names
fit_jacc <- left_join(fit_jacc, info, by = c("sample" = "library_name")) # join with metadata

var <- round(fit_jacc_all$eig*100/sum(fit_jacc_all$eig), 1)

study_plt <- ggplot(fit_jacc, aes(x = dim1, y = dim2, label = read_count,
                             color = study_accession, shape = diagnosis)) +
  geom_point() +
  theme_minimal() +
  ggtitle("IBD studies") +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  scale_shape(solid = F) +
  scale_color_brewer(palette = "Paired")

study_plt <- ggExtra::ggMarginal(study_plt, type = "density",
                    groupColour = T, groupFill = T)


pdf("sandbox/02_filt_comp_all/study_comp_plt.pdf",
    height = 8, width = 9)
study_plt
dev.off()



diag_plt_jacc <- ggplot(fit_jacc, aes(x = dim1, y = dim2, label = read_count,
                            color = diagnosis, shape = study_accession)) +
  geom_point() +
  theme_minimal() +
  #ggtitle("Jaccard distance") +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  scale_shape(solid = F, name = "Study") +
  scale_color_viridis_d(name = "Diagnosis") + 
  theme(legend.position = "none")

diag_plt_jacc <- ggExtra::ggMarginal(diag_plt_jacc, type = "density",
                    groupColour = T, groupFill = T)

pdf("sandbox/02_filt_comp_all/diag_comp_plt.pdf",
    height = 8, width = 9)
diag_plt_jacc
dev.off()


# tsne --------------------------------------------------------------------

library(Rtsne)
tsne_model <- Rtsne(comp, check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
tsne_y <- as.data.frame(tsne_model$Y) 
tsne_y$sample <- rownames(comp)
tsne_y<- left_join(tsne_y, info, by = c("sample" = "library_name")) # join with metadata
tsne_plt <- ggplot(tsne_y, aes(x = V1, y = V2, color = diagnosis, shape = study_accession)) +
  geom_point() +
  theme_minimal()  +
  labs(x = paste0("t-SNE 1"), 
       y = paste0("t-SNE 2")) + 
  scale_color_viridis_d()
tsne_plt

