library(dplyr)
library(rjson)
library(purrr)
library(tidyr)
library(ggplot2)

sigs <- Sys.glob("sandbox/test_charcoal/output.genbank_genomes/clean_sigs/*sig")
hash_to_assembly <- data.frame()
for(i in 1:length(sigs)){
  sig <- sigs[i]
  sig_json = fromJSON(file = sig)
  df <- as.data.frame(sig_json) 
  hash_to_assembly <- rbind(hash_to_assembly, df)
}
length(unique(hash_to_assembly$filename))

hash_to_assembly <- hash_to_assembly %>%
  mutate(name = gsub("_genomic.fna.gz.clean.fa.gz", "", filename))

# read in and normalize importances from models ---------------------------

read_varimp <- function(path_optimal_rf){
  study <- gsub("outputs\\/optimal_rf_seed\\/", "", path_optimal_rf)
  study <- gsub("_optimal_rf", "", study)
  study <- gsub(".RDS", "", study)
  optimal_rf <- readRDS(path_optimal_rf)
  varimp <- data.frame(hash = names(optimal_rf$variable.importance), 
                       importance = optimal_rf$variable.importance,
                       study = study)
  varimp <- separate(varimp, col = study, into = c("study", "seed"), sep = "_")
  rownames(varimp) <- seq(1:nrow(varimp))
  varimp$hash <- as.numeric(varimp$hash)
  # add a column where varimp is normalized by the total var imp
  # e.g., divide by the sum of all variable importances
  # this will make all variable importance measures sum to 1
  varimp$model_norm_imp <- varimp$importance / sum(varimp$importance)
  return(varimp)
}

varimp <- Sys.glob("outputs/optimal_rf_seed/*RDS") %>%
  map_dfr(read_varimp) %>%
  mutate(all_model_norm_imp = model_norm_imp / 36)

# join with assembly information ------------------------------------------

varimp <- left_join(varimp, hash_to_assembly, by = c("hash" ="signatures.mins"))
varimp <- varimp %>%
  mutate(accession = gsub("\\..*", "", name))

lineages <- read_csv("sandbox/test_gather_lineage_summarize/gtdb-rs202-genbank-protozoa-viral-fungi-lineage.csv")

varimp <- left_join(varimp, lineages, by = c("accession" = "ident"))


# label hashes that are contained in multiple genomes
# THIS DOESN'T DO WHAT IT'S SUPPOSED TO; this labels duplicates ALSO based on 
# presence in multiple models. A separate df needs to be created and grouped
# by hash & genome to ascertain duplicate status
varimp <- varimp %>%
  mutate(duplicate_hash = duplicated(hash)) 

# add a column that sums importance per genome

name_imp <- varimp %>%
  group_by(name) %>%
  summarize(total_all_model_norm_imp = sum(all_model_norm_imp))
varimp <- left_join(varimp, name_imp, by = "name")
# plot --------------------------------------------------------------------

varimp_tmp <- varimp %>% filter(!is.na(name))
pdf("tmp_varimp_clean.pdf", width = 10, height = 30)
ggplot(varimp_tmp, aes(x = reorder(species, total_all_model_norm_imp), 
                       y = all_model_norm_imp)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip()
dev.off()

# calculate percentage of total importance if 1% threshold is used --------

top_imp_genomes <- name_imp %>%
  filter(total_all_model_norm_imp >= 0.01) 

# plot genomes that hold >1% imp
varimp_tmp <- varimp %>% 
  filter(!is.na(name)) %>%
  filter(name %in% top_imp_genomes$name) 

# calc how many genomes each hash appears in
hash_dups <- varimp_tmp %>%
  select(hash, name) %>%
  distinct() %>%
  group_by(hash) %>% 
  mutate(duplicate_hash_times = seq(n()))

varimp_tmp <- left_join(varimp_tmp, hash_dups)
# create color vector for plotting
varimp_tmp <- varimp_tmp %>%
  mutate(duplicate_hash_times_color = ifelse(duplicate_hash_times >= 6, ">=6", duplicate_hash_times))

ggplot(varimp_tmp, aes(x = reorder(species, total_all_model_norm_imp), 
                       y = all_model_norm_imp, fill=as.character(duplicate_hash_times_color))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip()

ggplot(varimp_tmp, aes(x = reorder(species, total_all_model_norm_imp), 
                       y = all_model_norm_imp, fill=duplicate_hash)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip()


tmp <- varimp %>%
  filter(name %in% top_imp_genomes$name) %>%
  select(hash, all_model_norm_imp) %>%
  distinct()
sum(tmp$all_model_norm_imp)  
View(tmp)

tmp <- varimp %>%
  filter(!is.na(name)) %>%
  select(hash, all_model_norm_imp) %>%
  distinct()
sum(tmp$all_model_norm_imp)  

# lineage sheet ------------------------------------------------------------

charcoal_lineages <- varimp %>%
  select(name, superkingdom, phylum, class, order, family, genus, species) %>%
  distinct() %>%
  filter(!is.na(name)) %>%
  mutate(filename = gsub(" .*", "_genomic.fna.gz", name)) %>%
  mutate(superkingdom = gsub("Eukaryota", "d__Eukaryota", superkingdom)) %>%
  select(filename, superkingdom, phylum, class, order, family, genus, species)
write_csv(charcoal_lineages, "sandbox/test_charcoal/lineages.csv", col_names = F)

