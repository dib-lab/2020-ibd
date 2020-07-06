library(readr)
library(dplyr)
sig_ccs <- read_tsv(snakemake@input[['sig']])
uc <- sig_ccs %>%
  filter(mu == "mu.diagnosisUC")

cd <- sig_ccs %>%
  filter(mu == "mu.diagnosisCD")

write.table(uc$aa_seq, snakemake@output[['uc']],
            quote = F, row.names = F, col.names = F)

write.table(cd$aa_seq, snakemake@output[['cd']],
            quote = F, row.names = F, col.names = F)
