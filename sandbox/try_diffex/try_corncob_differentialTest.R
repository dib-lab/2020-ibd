setwd("~/github/2020-ibd/")

library(dplyr)
library(readr)
library(tximport)
library(corncob)
library(phyloseq)
library(edgeR)

files <- list.files("sandbox/try_diffex", "quant.sf$", recursive = T, full.names = T, )
files_root <- gsub("\\/SRS294916_20\\.fna_quant\\/quant\\.sf", "", files)
files_root <- gsub("sandbox\\/try_diffex\\/", "", files_root)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct() %>%
  filter(library_name %in% files_root)

reads <- read_tsv("sandbox/try_diffex/mqc_fastp_filtered_reads_plot_1.txt",
                  col_names = c("library_name", "passed_filter", "low_quality", "too_many_N", "too_short"),
                  skip = 1) %>%
  mutate(library_name = gsub("\\.abundtrim", "", library_name)) %>%
  group_by(library_name) %>%
  mutate(mm = sum(passed_filter, low_quality, too_many_N, too_short)) %>%
  select(library_name, mm)


info <- left_join(info, reads, by = "library_name")

# match order of info to order of files_root
info <- info[order(match(info$library_name, files_root)), ]

# check that order matches
all.equal(info$library_name, files_root)

info$salmon <- files

counts <- tximport(files = info$salmon, type = "salmon", txOut = T)
colnames(counts$counts) <- info$library_name


### ADD
info <- select(info, study_accession, library_name, diagnosis)
colnames(counts$counts) <- info$library_name
count_info <- counts$counts

# filter to remove amino acid seqs that are 0 across an entire group
info <- info[order(match(info$library_name, colnames(count_info))), ]
# check that order matches
all.equal(info$library_name, colnames(count_info))
tmp <- filterByExpr(y = as.matrix(count_info), 
                    group = paste0(info$study_accession, "_", info$diagnosis),
                    min.count = 1, min.total.count = 2)
tmp2 <- count_info[tmp, ]
count_info <- tmp2
dim(count_info)

# make counts integers
count_info <- apply(count_info, 1:2, round) 
####



## Create phyloseq object
my_counts <- otu_table(count_info, taxa_are_rows = TRUE)
# Pull out the gene ID to use as "taxa" names
taxa_names(my_counts) <- rownames(count_info)

# Remove first three columns
my_samp <- sample_data(select(info, study_accession, diagnosis))
sample_names(my_samp) <- info %>% pull(library_name)
ibd <- phyloseq(my_counts, my_samp)


my_output <- differentialTest(formula = ~ study_accession + diagnosis,
                              formula_null = ~ study_accession,
                              phi.formula = ~ 1,
                              phi.formula_null = ~ 1,
                              data = ibd, 
                              test = "LRT", boot = FALSE)

length(my_output$significant_taxa)
## Plot!
my_gene_counts <- count_info["20",]
my_gene_counts <- data.frame(x20 = my_gene_counts, sample = names(my_gene_counts))
my_gene_counts <- left_join(my_gene_counts, info, by = c("sample" = "library_name"))
ggplot(my_gene_counts, aes(x = diagnosis, y = x20)) +
  geom_violin() +
  theme_minimal()

pdf("tmp_corncob.pdf", height = 300, width = 5)
plot(my_output)
dev.off()
