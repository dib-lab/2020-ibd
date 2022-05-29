setwd("~/github/ibd/")

library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(viridis)
library(readr)
library(tidyr)
library(ggrepel)
library(stringr)

# read in jaccard ---------------------------------------------------------

comp_jacc <- read_csv("sandbox/02_filt_comp_all/comp_all_jaccard.csv") # read in mat
colnames(comp_jacc) <- gsub("_filt", "", colnames(comp_jacc)) # remove filt from names
rownames(comp_jacc) <-colnames(comp_jacc)                     # set rownames
# read in and formate metadata -------------------------------------------

info <- read_tsv("inputs/working_metadata.tsv") %>%
  filter(library_name %in% colnames(comp_jacc)) %>%
  group_by(library_name, study_accession, diagnosis) %>%
  summarise(read_count = sum(read_count)) %>%
  as.data.frame()

info_hmp <- read_tsv("inputs/hmp2_mgx_metadata.tsv") %>%
  mutate(study_accession = "iHMP") %>%
  select(External.ID, study_accession, diagnosis, reads_raw) %>%
  distinct() %>%
  filter(External.ID %in% colnames(comp_jacc))
colnames(info_hmp) <- c("library_name", "study_accession", "diagnosis", "read_count")

info <- rbind(info, info_hmp)


# filter comp -------------------------------------------------------------

# filter the comparison matrix to remove the two samples that failed processing

keep <- colnames(comp_jacc)[colnames(comp_jacc) %in% info$library_name]
comp_jacc <- select(comp_jacc, keep)
comp_jacc <- comp_jacc[rownames(comp_jacc) %in% keep, ]
dim(comp_jacc)
rownames(comp_jacc) <- colnames(comp_jacc)
# plot jaccard ------------------------------------------------------------

dist_jacc <- dist(comp_jacc)                        # compute distances
fit_jacc_all <- cmdscale(dist_jacc, eig = T)        # calculate MDS
fit_jacc <- as.data.frame(fit_jacc_all$points)
fit_jacc$sample <- rownames(fit_jacc)

colnames(fit_jacc) <- c("dim1", "dim2", "sample")        # set column names
fit_jacc <- left_join(fit_jacc, info, 
                      by = c("sample" = "library_name")) # join with metadata

var_jacc <- round(fit_jacc_all$eig*100/sum(fit_jacc_all$eig), 1)

diag_plt_jacc <- ggplot(fit_jacc, aes(x = dim1, y = dim2, color = diagnosis, shape = study_accession)) +
  geom_point() +
  theme_minimal() +
  #ggtitle("Jaccard distance") +
  labs(x = paste0("PCo 1 (", var_jacc[1], "%)"),
       y = paste0("PCo 2 (", var_jacc[2], "%)")) +
  scale_shape(solid = F, name = "Study") +
  scale_color_viridis_d(name = "Diagnosis") + 
  theme(legend.position = "none")

diag_plt_jacc <- ggMarginal(diag_plt_jacc, type = "density",
                            groupColour = T, groupFill = T)

# pdf("sandbox/02_filt_comp_all/diag_comp_jacc_plt.pdf",
#     height = 8, width = 9)
# diag_plt_jacc
# dev.off()

# study_plt_jacc <- ggplot(fit_jacc, aes(x = dim1, y = dim2, color = study_accession, shape = diagnosis)) +
#   geom_point() +
#   theme_minimal() +
#   labs(x = paste0("PCo 1 (", var[1], "%)"),
#        y = paste0("PCo 2 (", var[2], "%)")) +
#   scale_shape(solid = F) +
#   scale_color_brewer(palette = "Paired")
# 
# study_plt_jacc <- ggExtra::ggMarginal(study_plt_jacc, type = "density",
#                                  groupColour = T, groupFill = T)
# 
# 
# pdf("sandbox/02_filt_comp_all/study_comp_jacc_plt.pdf",
#     height = 8, width = 9)
# study_plt_jacc
# dev.off()

# prep cosine -------------------------------------------------------------

comp_cos <- read_csv("outputs/comp/all_filt_comp.csv") # read in mat
colnames(comp_cos) <- gsub("_filt", "", colnames(comp_cos))
rownames(comp_cos) <-colnames(comp_cos) 
dist_cos <- dist(comp_cos)                                     # compute distances


fit_all_cos <- cmdscale(dist_cos, eig = T)                             # calculate MDS
fit_cos <- as.data.frame(fit_all_cos$points)
fit_cos$sample <- rownames(fit_cos)

colnames(fit_cos) <- c("dim1", "dim2", "sample")                   # set column names
fit_cos <- left_join(fit_cos, info, by = c("sample" = "library_name")) # join with metadata

var_cos <- round(fit_all_cos$eig*100/sum(fit_all_cos$eig), 1)

diag_plt_cosine <- ggplot(fit_cos, aes(x = dim1, y = dim2, color = diagnosis, shape = study_accession)) +
  geom_point() +
  theme_minimal() +
  #ggtitle("Cosine distance") +
  labs(x = paste0("PCo 1 (", var_cos[1], "%)"),
       y = paste0("PCo 2 (", var_cos[2], "%)")) +
  scale_shape(solid = F, name = "Study") +
  scale_color_viridis_d(name = "Diagnosis") +
  theme(legend.position = "none")

diag_plt_cosine <- ggMarginal(diag_plt_cosine, type = "density",
                              groupColour = T, groupFill = T)

# pdf("sandbox/02_filt_comp_all/diag_cosine_comp_plt.pdf", height = 8, width = 9)
# diag_plt_cosine
# dev.off()

# study_plt <- ggplot(fit_cos, aes(x = dim1, y = dim2, color = study_accession, shape = diagnosis)) +
#   geom_point() +
#   theme_minimal() +
#   #ggtitle("IBD studies") +
#   labs(x = paste0("PCo 1 (", var_cos[1], "%)"),
#        y = paste0("PCo 2 (", var_cos[2], "%)")) +
#   scale_shape(solid = F) +
#   scale_color_brewer(palette = "Paired")
# 
# study_plt <- ggExtra::ggMarginal(study_plt, type = "density",
#                                  groupColour = T, groupFill = T)
# 
# pdf("sandbox/02_filt_comp_all/study_cosine_comp_plt.pdf",
#     height = 8, width = 9)
# study_plt
# dev.off()

# make legend -------------------------------------------------------------

# make with cosine
legend <- ggplot(fit_cos, aes(x = dim1, y = dim2, color = diagnosis, shape = study_accession)) +
  geom_point() +
  scale_shape(solid = F, name = "Study") +
  scale_color_viridis_d(name = "Diagnosis")
legend <- as_ggplot(get_legend(legend))

# make with jacc, same as cos but doesn't require cosine section to be run
legend <- ggplot(fit_jacc, aes(x = dim1, y = dim2, color = diagnosis, shape = study_accession)) +
  geom_point() +
  scale_shape(solid = F, name = "Study") +
  scale_color_viridis_d(name = "Diagnosis") + 
  theme_minimal() +
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=13))
legend <- as_ggplot(get_legend(legend))
# plot all ----------------------------------------------------------------

plt_all <- ggarrange(as_ggplot(diag_plt_jacc), as_ggplot(diag_plt_cosine), legend,
                     labels = c("A", "B"),
                     ncol = 3, heights = 3, widths = c(2, 2, 1))

pdf("sandbox/02_filt_comp_all/minhash-comp-plts-all.pdf", height = 3.5, width = 8)
plt_all 
dev.off()

pdf("figures/minhash-comp-plts-all.pdf", height = 3.5, width = 8)
plt_all 
dev.off()


plt_jacc <- ggarrange(as_ggplot(diag_plt_jacc), legend, ncol = 2, heights = 3,
                      widths = c(2, 1))
plt_jacc

pdf("figures/minhash-comp-plt-jacc-all.pdf", height = 3.5, width = 6)
plt_jacc
dev.off()
