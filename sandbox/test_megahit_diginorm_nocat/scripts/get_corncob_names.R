library(readr)

sig_ccs <- read_tsv(snakemake@input[['sig']])
aa_seq <- gsub("^V", "", sig_ccs$aa_seq)
write.table(aa_seq, snakemake@output[['names']],
            quote = F, row.names = F, col.names = F)
