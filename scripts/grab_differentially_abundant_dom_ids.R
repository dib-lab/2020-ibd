library(dplyr)
library(readr)

#sig_ccs <- read_tsv("outputs/sgc_pangenome_catlases_corncob/GCA_000210075.1_sig_ccs.tsv") %>%
sig_ccs <- read_tsv(snakemake@input[[sig_ccs]]) %>%
  select(mu, estimate, dom_id = aa_seq) %>%
  mutate(mu = gsub("mu.diagnosis", "", mu)) %>%
  mutate(abundance = ifelse(estimate > 0, "increased", "decreased"))

# write out dominators
sig_ccs %>%
  select(-estimate) %>%
  group_by(mu, abundance) %>%
  #group_walk(~ write_tsv(.x, paste0(.y$mu, "_", .y$abundance, "_nbhds.gz"), col_names = F))
  group_walk(~ write_tsv(.x, paste0(snakemake@params[['outdir']], 
                                    snakemake@wildcard[['acc']], "_",
                                    .y$mu, "_", .y$abundance, "_dom_ids.tsv.gz"), 
                         col_names = F))

# create empty file if no matching criteria
for(file in unlist(snakemake@output)){
  if (!file.exists(file)) {
    file.create(file)
}
