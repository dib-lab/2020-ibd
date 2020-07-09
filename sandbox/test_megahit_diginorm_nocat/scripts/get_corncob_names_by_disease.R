library(readr)
library(dplyr)
sig_ccs <- read_tsv(snakemake@input[['sig']])
uc_up <- sig_ccs %>%
  filter(mu == "mu.diagnosisUC") %>%
  filter(estimate > 0)

uc_down <- sig_ccs %>%
  filter(mu == "mu.diagnosisUC") %>%
  filter(estimate < 0)

cd_up <- sig_ccs %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate > 0) 

cd_down <- sig_ccs %>%
  filter(mu == "mu.diagnosisCD") %>%
  filter(estimate < 0)

write.table(uc_up$aa_seq, snakemake@output[['uc_up']],
            quote = F, row.names = F, col.names = F)

write.table(uc_down$aa_seq, snakemake@output[['uc_down']],
            quote = F, row.names = F, col.names = F)

write.table(cd_up$aa_seq, snakemake@output[['cd_up']],
            quote = F, row.names = F, col.names = F)

write.table(cd_down$aa_seq, snakemake@output[['cd_down']],
            quote = F, row.names = F, col.names = F)
