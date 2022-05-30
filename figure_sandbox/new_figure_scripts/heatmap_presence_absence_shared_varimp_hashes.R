# plot 934 shared hashes as a presence absence heatmap, annotated with diagnosis
library(readr)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(viridis)
# define the intersection of all hashes that had a variable importance across
# all models
prjeb2054 <- read_tsv("outputs/vita_rf/PRJEB2054_vita_vars.txt", col_names = "hash")
prjna237362 <- read_tsv("outputs/vita_rf/PRJNA237362_vita_vars.txt", col_names = "hash")
ihmp <- read_tsv("outputs/vita_rf/iHMP_vita_vars.txt", col_names = "hash")
srp057027 <- read_tsv("outputs/vita_rf/SRP057027_vita_vars.txt", col_names = "hash")
prjna385949 <- read_tsv("outputs/vita_rf/PRJNA385949_vita_vars.txt", col_names = "hash")
prjna400072 <- read_tsv("outputs/vita_rf/PRJNA400072_vita_vars.txt", col_names = "hash")

intersect_hashes_all <- Reduce(intersect, list(prjeb2054$hash, prjna237362$hash,
                                               ihmp$hash, srp057027$hash, 
                                               prjna400072$hash, prjna385949$hash))

# read in abund info
abund <- read_csv("outputs/vita_rf/SRP057027_ibd_filt.csv")

# filter to only hashes that are predicted to be important in all models
abund_small <- abund %>%
  select(X1, all_of(as.character(intersect_hashes_all)))
abund_small <- as.data.frame(abund_small)
rownames(abund_small) <- abund_small$X1
abund_small <- abund_small[ , -1]


# set to presence absence -------------------------------------------------
abund_small[abund_small > 0] <- 1

# read info 
info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct()%>%
  mutate(library_name = gsub("\\-", ".", library_name))

table(rownames(abund_small) %in% info$library_name)

info <- info[match(rownames(abund_small), info$library_name), ] # match order of classification to mtx
#tmp <- info %>%
#  select_if(function(x) !(all(is.na(x)) | all(x=="")))

mat_row <- data.frame(diagnosis = info$diagnosis)
rownames(mat_row) <- rownames(abund_small)
#tmp_pal <- brewer.pal(12, "Paired")
tmp_pal <- viridis(3)
mat_colors <- list(diagnosis = tmp_pal)
names(mat_colors$diagnosis) <- unique(mat_row$diagnosis)

pdf("tmp.pdf", height = 25, width = 10)
pheatmap(mat                  = abund_small, 
         show_rownames        = T,
         fontsize_row         = 5,
         show_colnames        = F, 
         cluster_rows         = T, 
         cluster_cols         = T,
         color                = brewer.pal(11, "BrBG")[c(5, 1)],
         annotation_row       = mat_row,
         annotation_colors    = mat_colors,
         legend_breaks        = c(0, 1),
         legend_labels        = c("Absent", "Present"),
         fontsize             = 14,
         annotation_names_row = F)
dev.off()


