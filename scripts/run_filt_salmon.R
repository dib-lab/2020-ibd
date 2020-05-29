library(readr)
library(dplyr)
library(edgeR)


info <- read_tsv(snakemake@input[["info"]]) %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct() 

count_info <- read_tsv(snakemake@input[['counts']])

info_salmon <- info %>%
  select(colnames(count_info))

# filter to remove amino acid seqs that are 0 across an entire group
info_salmon <- info_salmon[order(match(info_salmon$library_name, colnames(count_info))), ]
# check that order matches
stopifnot(all.equal(info_salmon$library_name, colnames(count_info)))
keep <- filterByExpr(y = as.matrix(count_info), 
                     group = paste0(info_salmon$study_accession, "_", info_salmon$diagnosis),
                     min.count = 1, min.total.count = 2)
count_info <- count_info[keep, ]
write_tsv(dim(count_info), path = snakemake@output[["dim"]])

# make counts integers
count_info <- apply(count_info, 1:2, round) 

# re-add columns that had no mapping reads from salmon step
no_salmon <- info$library_name[!info$library_name  %in% info_salmon$library_name]
no_salmon_mat  <- matrix(ncol = length(no_salmon), nrow = nrow(count_info), data = 0)
colnames(no_salmon_mat) <- no_salmon
count_info <- cbind(count_info, no_salmon_mat)
write_tsv(count_info, path = snakemake@output[["filt"]])
