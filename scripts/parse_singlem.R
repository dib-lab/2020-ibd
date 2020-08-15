library(readr)

defaults <- read_tsv(snakemake@input[['default']], col_types = "ccccddcccl")
s16 <- read_tsv(snakemake@input[['s16']], col_types = "ccccddcccl")

all_res <- rbind(defaults, s16)
#write_tsv(all_res, "singlem_combined.tsv")
write_tsv(all_res, snakemake@output[['res']])
