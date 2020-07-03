library(readr)

sig_ccs <- read_tsv(snakemake@input[['sig']])
write.table(sig_ccs$aa_seq, snakemake@output[['names']],
            quote = F, row.names = F, col.names = F)
