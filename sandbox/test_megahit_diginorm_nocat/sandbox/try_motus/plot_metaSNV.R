# pot metaSNV

library(dplyr)
library(ggplot)
library(readr)
library(tibble)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct()

allele_dist <- read_tsv("sandbox/test_megahit_diginorm_nocat/sandbox/try_motus/snv_call_out3/distances-m2-d1-b10-c1-p0.1/ref_mOTU_v25_01594.filtered.allele.dist") %>%
  mutate(X1 = gsub("\\.bam", "", X1)) %>%
  column_to_rownames("X1")
colnames(allele_dist) <- gsub("\\.bam", "", colnames(allele_dist))

info <- info %>%
  filter(library_name %in% colnames(allele_dist))
fit_all <- cmdscale(allele_dist, eig = T)                              # calculate MDS
fit <- as.data.frame(fit_all$points)                            # abstract data
fit$sample <- rownames(fit)                                     # set rownames
colnames(fit) <- c("dim1", "dim2", "sample")                    # set column names
fit <- left_join(fit, info, by = c("sample" = "library_name"))  # join with metadata
var <- round(fit_all$eig*100/sum(fit_all$eig), 1)               # calc percent var

plt <- ggplot(fit, aes(x = dim1, y = dim2, color = diagnosis, label = sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("CD" = "#a52a2a", "nonIBD" = "#6494ed", "UC" = "#ffa600", name = "Diagnosis"))

ggplotly(plt)  
