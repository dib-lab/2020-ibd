library(readr)
library(dplyr)
library(ggplot2)
library(tibble)
library(Rtsne)
library(ggpubr)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, diagnosis) %>%
  #mutate(library_name = gsub("-", "\\.", library_name)) %>%
  distinct()

# c bolteae ---------------------------------------------------------------

cbolt_k21 <- read_csv("outputs/tmp_comp/GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna.cdbg_ids.reads.sig_k21_ignore_abund_comp.csv")
colnames(cbolt_k21) <- gsub("_GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna", "", colnames(cbolt_k21))
rownames(cbolt_k21) <- colnames(cbolt_k21)
cbolt_k21 <- as.matrix(cbolt_k21)

fit <- dist(cbolt_k21)
fit <- cmdscale(fit)
fit <- fit %>%
  as.data.frame() %>%
  rownames_to_column("library_name") %>%
  left_join(info)
cbolt_mds <- ggplot(fit, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  ggtitle("C bolt nMDS") + 
  theme_minimal()


tsne_model <- Rtsne(cbolt_k21, check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
d_tsne = as.data.frame(tsne_model$Y) 
d_tsne$library_name <- rownames(cbolt_k21)
d_tsne <- d_tsne %>%
  left_join(info)
cbolt_tsne <- ggplot(d_tsne, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  ggtitle("C bolt tSNE") + 
  theme_minimal()


# r gnavus ----------------------------------------------------------------

rgnv_k21 <- read_csv("outputs/tmp_comp/GCF_900036035.1_RGNV35913_genomic.fna.cdbg_ids.reads.sig_k21_ignore_abund_comp.csv")
colnames(rgnv_k21) <- gsub("_GCF_900036035.1_RGNV35913_genomic.fna", "", colnames(rgnv_k21))
rownames(rgnv_k21) <- colnames(rgnv_k21)
rgnv_k21 <- as.matrix(rgnv_k21)

fit <- dist(rgnv_k21)
fit <- cmdscale(fit)
fit <- fit %>%
  as.data.frame() %>%
  rownames_to_column("library_name") %>%
  left_join(info)
rgnv_mds <- ggplot(fit, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  ggtitle("R gnvs nMDS") + 
  theme_minimal()

tsne_model <- Rtsne(rgnv_k21, check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
d_tsne = as.data.frame(tsne_model$Y) 
d_tsne$library_name <- rownames(rgnv_k21)
d_tsne <- d_tsne %>%
  left_join(info)
rgnv_tsne <-ggplot(d_tsne, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  ggtitle("R gnvs tSNE") + 
  theme_minimal()

# r gnavus ----------------------------------------------------------------

fsp_k21 <- read_csv("outputs/tmp_comp/GCF_000508885.1_ASM50888v1_genomic.fna.cdbg_ids.reads.sig_k21_ignore_abund_comp.csv")
colnames(fsp_k21) <- gsub("_GCF_000508885.1_ASM50888v1_genomic.fna", "", colnames(fsp_k21))
rownames(fsp_k21) <- colnames(fsp_k21)
fsp_k21 <- as.matrix(fsp_k21)

fit <- dist(fsp_k21)
fit <- cmdscale(fit)
fit <- fit %>%
  as.data.frame() %>%
  rownames_to_column("library_name") %>%
  left_join(info)
fsp_mds <- ggplot(fit, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  ggtitle("F sp nMDS") + 
  theme_minimal()


tsne_model <- Rtsne(fsp_k21, check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
d_tsne = as.data.frame(tsne_model$Y) 
d_tsne$library_name <- rownames(fsp_k21)
d_tsne <- d_tsne %>%
  left_join(info)
fsp_tsne <- ggplot(d_tsne, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  ggtitle("F sp tSNE") + 
  theme_minimal()

# f plaut ----------------------------------------------------------------

fpl_k21 <- read_csv("outputs/tmp_comp/ERS537353_12.fna.cdbg_ids.reads.sig_k21_ignore_abund_comp.csv")
colnames(fpl_k21) <- gsub("_ERS537353_12.fna", "", colnames(fpl_k21))
rownames(fpl_k21) <- colnames(fpl_k21)
fpl_k21 <- as.matrix(fpl_k21)

fit <- dist(fpl_k21)
fit <- cmdscale(fit)
fit <- fit %>%
  as.data.frame() %>%
  rownames_to_column("library_name") %>%
  left_join(info)
fpl_mds <- ggplot(fit, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  ggtitle("F plautii nMDS") + 
  theme_minimal()


tsne_model <- Rtsne(fpl_k21, check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
d_tsne = as.data.frame(tsne_model$Y) 
d_tsne$library_name <- rownames(fpl_k21)
d_tsne <- d_tsne %>%
  left_join(info)
fpl_tsne <- ggplot(d_tsne, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  ggtitle("F plautii tSNE") + 
  theme_minimal()

ggarrange(fpl_mds, fpl_tsne,
          fsp_mds, fsp_tsne,
          nrow = 2, ncol = 2, common.legend = T)

ggarrange(rgnv_mds, rgnv_tsne,
          cbolt_mds, cbolt_tsne,
          nrow = 2, ncol = 2, common.legend = T)
