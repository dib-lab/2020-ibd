library(readr)
library(tibble)
library(dplyr)
library(edgeR)

# read in metadata ---------------------------------------
info <- read_tsv(snakemake@input[["info"]]) %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct() 

# read in abundtrim lib size info
libsizes <- read_tsv(snakemake@input[['mqc_fastp']], skip = 1,
                  col_names = c("library_name", "passed_filter", "low_quality", 
                                "too_many_N", "too_short")) %>%
  mutate(library_name = gsub("\\.abundtrim", "", library_name)) %>%
  group_by(library_name) %>%
  mutate(libsize = sum(passed_filter, low_quality, too_many_N, too_short)) %>%
  select(library_name, libsize)

# join metadata
info <- left_join(info, libsizes, by = "library_name")

# read in and filter counts ------------------------------
count_info <- read_tsv(snakemake@input[['counts']])
rownames(count_info) <- count_info$protein
count_info <- as.data.frame(count_info[ , -1])

info_salmon <- info %>%
  filter(library_name %in% colnames(count_info))

# filter to remove amino acid seqs that are 0 across an entire group
info_salmon <- info_salmon[order(match(info_salmon$library_name, colnames(count_info))), ]
# check that order matches
stopifnot(all.equal(info_salmon$library_name, colnames(count_info)))
keep <- filterByExpr(y = as.matrix(count_info), 
                     group = paste0(info_salmon$study_accession, "_", info_salmon$diagnosis),
                     lib.size = info_salmon$libsize,
                     min.count = 1, min.total.count = 2)
count_info <- count_info[keep, ]
write_tsv(as.data.frame(dim(count_info)), path = snakemake@output[["dim"]])

# make counts integers
count_info <- apply(count_info, 1:2, round) 

# re-add columns that had no mapping reads from salmon step
no_salmon <- info$library_name[!info$library_name  %in% info_salmon$library_name]
no_salmon_mat  <- matrix(ncol = length(no_salmon), nrow = nrow(count_info), data = 0)
colnames(no_salmon_mat) <- no_salmon
count_info <- cbind(count_info, no_salmon_mat)
count_info <- count_info %>%
  as.data.frame() %>%
  add_column("protein" = rownames(count_info), .before=TRUE)
# write full counts to file
write_tsv(count_info, path = snakemake@output[["filt"]])
