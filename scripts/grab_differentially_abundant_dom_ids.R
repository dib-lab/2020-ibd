library(dplyr)
library(readr)

if (!dir.exists(snakemake@params[['outdir']])) {dir.create(snakemake@params[['outdir']])}

cdbg_to_pieces <- read_csv(snakemake@input[['cdbg_to_pieces']])

#sig_ccs <- read_tsv("outputs/sgc_pangenome_catlases_corncob/GCA_000210075.1_sig_ccs.tsv") %>%
sig_ccs <- read_tsv(snakemake@input[['sig_ccs']]) %>%
  select(mu, estimate, dom_id = aa_seq) %>%
  mutate(mu = gsub("mu.diagnosis", "", mu)) %>%
  mutate(abundance = ifelse(estimate > 0, "increased", "decreased")) %>%
  left_join(cdbg_to_pieces, by = c("dom_id" = "dominator")) 

# write out dominators
sig_ccs %>%
  select(-estimate, -dom_id) %>%
  group_by(mu, abundance) %>%
  #group_walk(~ write_tsv(.x, paste0(.y$mu, "_", .y$abundance, "_nbhds.gz"), col_names = F))
  group_walk(~ write_tsv(.x, paste0(snakemake@params[['outdir']], 
                                    snakemake@wildcards[['acc']], "_",
                                    .y$mu, "_", .y$abundance, "_cdbg_ids.tsv.gz"), 
                         col_names = F))

# create empty file if no matching criteria
for(file in unlist(snakemake@output)){
  if (!file.exists(file)) {file.create(file)}
}

