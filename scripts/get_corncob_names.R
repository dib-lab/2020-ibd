library(readr)

sig_ccs <- read_tsv(snakemake@input)
write.table(sig_ccs$aa_seq, snakemake@output,
            quote = F, row.names = F, col.names = F)
