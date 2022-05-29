library(readr)
library(dplyr)
library(ggplot2)
library(plotly)

comp <- read_csv("sandbox/hash_corncob/gather_containment_comp/GCF_000371685.1_Clos_bolt_90B3_V1_genomic_cd_up_hashes_matches_comp.csv")
containment <- read_csv("sandbox/hash_corncob/gather_containment_results/GCF_000371685.1_Clos_bolt_90B3_V1_genomic_cd_up_hashes_containment.txt")

comp_high <- comp %>% 
  transmute(sample = colnames(comp), across(everything())) %>%
  select(sample, "DS480726.1 Clostridium bolteae ATCC BAA-613 Scfld_02_67 genomic scaffold, whole genome shotgun sequence",
         "KQ235845.1 [Clostridium] bolteae WAL-14578 genomic scaffold supercont1.1, whole genome shotgun sequence",
         "KB851183.1 Clostridium bolteae 90A5 genomic scaffold acBRZ-supercont1.1, whole genome shotgun sequence",
         "KB851159.1 Clostridium bolteae 90B7 genomic scaffold acsNV-supercont1.1, whole genome shotgun sequence",
         "KB851138.1 Clostridium bolteae 90B8 genomic scaffold acBSf-supercont1.1, whole genome shotgun sequence") %>%
  left_join(containment, by = c("sample" = "name"))
  
View(comp_high)

dist_comp <- dist(comp)
fit_all <- cmdscale(dist_comp, eig = T)                         # calculate MDS
fit <- as.data.frame(fit_all$points)                            # abstract data
fit$sample <- colnames(comp)                                    # set rownames
colnames(fit) <- c("dim1", "dim2", "sample")                    # set column names
var <- round(fit_all$eig*100/sum(fit_all$eig), 1)               # calc percent var
fit$color <- ifelse(fit$sample %in% colnames(comp_high), "high", "lower")
fit$color <- ifelse(fit$sample %in% c("BBZM01000001.1 Clostridia bacterium UC5.1-1F7 DNA, contig: contig00001, whole genome shotgun sequence",
                                      "KB851045.1 Clostridium clostridioforme 90A7 genomic scaffold acsOa-supercont1.1, whole genome shotgun sequence",
                                      "KB851178.1 Clostridium bolteae 90B3 genomic scaffold acBRc-supercont1.1, whole genome shotgun sequence",
                                      "BAHX01000115.1 Clostridium sp. VE202-10 DNA, contig: VE202-10_contig00115, whole genome shotgun sequence",
                                      "HF988233.1 Clostridium bolteae CAG:59 genomic scaffold, scf3, whole genome shotgun sequence"),
                    "out_candidate", fit$color)
fit <- left_join(fit, comp_high, by = "sample")

fit_isolates <- fit %>%
  filter(!grepl(pattern = "^ERR", sample)) %>%
  filter(!grepl(pattern = "^SRR", sample)) %>%
  filter(!grepl(pattern = "^DRR", sample)) %>%
  filter(!grepl(pattern = "^HGM", sample)) %>%
  filter(!grepl(pattern = "^\\.\\/", sample)) 



View(fit_isolates)
plt <- ggplot(fit_isolates, aes(x = dim1, y = dim2, color = color, size = similarity, label = sample)) +
  geom_point(alpha = .1) +
  theme_minimal() +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5))
ggplotly(plt)  

# rgnv --------------------------------------------------------------------

comp <- read_csv("sandbox/hash_corncob/gather_containment_comp/GCF_900036035.1_RGNV35913_genomic_cd_up_hashes_matches_comp.csv")
containment <- read_csv("sandbox/hash_corncob/gather_containment_results/GCF_900036035.1_RGNV35913_genomic_cd_up_hashes_containment.txt")

comp_high <- comp %>% 
  transmute(sample = colnames(comp), across(everything())) %>%
  select(sample, "SRR3131756_bin.8.fa.gz", "ERR525797_bin.4.fa.gz", 
         "./4584/BackhedF_2015__SID258_4M__bin.6.fa",
         "FCFA01000001.1 [Ruminococcus] gnavus genome assembly, contig: contig1_1, whole genome shotgun sequence",
         "AAYG02000043.1 Ruminococcus gnavus ATCC 29149 R_gnavus-1.0.1_Cont380, whole genome shotgun sequence") %>%
  left_join(containment, by = c("sample" = "name"))

dist_comp <- dist(comp)
fit_all <- cmdscale(dist_comp, eig = T)                         # calculate MDS
fit <- as.data.frame(fit_all$points)                            # abstract data
fit$sample <- colnames(comp)                                    # set rownames
colnames(fit) <- c("dim1", "dim2", "sample")                    # set column names
var <- round(fit_all$eig*100/sum(fit_all$eig), 1)               # calc percent var
fit$color <- ifelse(fit$sample %in% colnames(comp_high), "high", "lower")
# fit$color <- ifelse(fit$sample %in% c("SRR3131756_bin.8.fa.gz", "ERR525797_bin.4.fa.gz", 
#                                       "./4584/BackhedF_2015__SID258_4M__bin.6.fa",
#                                       "FCFA01000001.1 [Ruminococcus] gnavus genome assembly, contig: contig1_1, whole genome shotgun sequence",
#                                       "AAYG02000043.1 Ruminococcus gnavus ATCC 29149 R_gnavus-1.0.1_Cont380, whole genome shotgun sequence"),
#                     "out_candidate", fit$color)
fit <- left_join(fit, comp_high, by = "sample")

# fit_isolates <- fit %>%
#   filter(!grepl(pattern = "^ERR", sample)) %>%
#   filter(!grepl(pattern = "^SRR", sample)) %>%
#   filter(!grepl(pattern = "^DRR", sample)) %>%
#   filter(!grepl(pattern = "^HGM", sample)) %>%
#   filter(!grepl(pattern = "^\\.\\/", sample)) 

fit_tmp <- fit %>%
  filter(similarity > .7)

plt <- ggplot(fit_tmp, aes(x = dim1, y = dim2, color = color, size = similarity, label = sample)) +
  geom_point(alpha = .05) +
  theme_minimal() +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5))
ggplotly(plt)  

