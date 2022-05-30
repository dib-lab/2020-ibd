library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(tibble)
library(vegan)
files <- list.files(path = "outputs/nbhd_reads_sigs_csv", pattern = ".csv$",
                    full.names = TRUE, recursive = TRUE)

# gsub_pattern <- paste0("\\/", snakemake@params[['gather_genome']]), ".cdbg_ids.reads.csv")
# values_from_pattern <- paste0(snakemake@params[['gather_genome']]), ".cdbg_ids.reads.sig")

sigs_long <- files %>%
  set_names() %>% 
  map_dfr(read_csv, .id = "library_name", col_types = "cd") %>%
  mutate(library_name = gsub("outputs\\/nbhd_reads_sigs_csv\\/", "", library_name)) %>%
  mutate(library_name = gsub("\\/LoombaR_2017__SID1050_bax__bin.11.fa.cdbg_ids.reads.csv", "", library_name))

sigs <- pivot_wider(sigs_long, id_cols = library_name, 
                    values_from = "LoombaR_2017__SID1050_bax__bin.11.fa.cdbg_ids.reads.sig",
                    names_from = "minhash")
sigs[is.na(sigs)] <- 0

counts <- sigs %>%
  column_to_rownames(var = "library_name")


# write_tsv()


info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct() %>%
  filter(library_name %in% rownames(counts))


# plot_species_accumulation <- function(info, counts, pangenome_name){
  info <- info[order(match(info$library_name, rownames(counts))), ]
  stopifnot(all.equal(info$library_name, rownames(counts)))
  
  uc <- filter(info, diagnosis == "UC")
  sp_uc <- specaccum(counts[rownames(counts) %in% uc$library_name, ], "exact")
  sp_uc_df <- data.frame(num_samples = sp_uc$sites, proteins = sp_uc$richness, sd = sp_uc$sd) 
  
  cd <- filter(info, diagnosis == "CD")
  sp_cd <- specaccum(counts[rownames(counts) %in% cd$library_name, ], "exact")
  sp_cd_df <- data.frame(num_samples = sp_cd$sites, proteins = sp_cd$richness, sd = sp_cd$sd)
  
  nonIBD <- filter(info, diagnosis == "nonIBD")
  sp_nonibd <- specaccum(counts[rownames(counts) %in% nonIBD$library_name, ], "exact")
  sp_nonibd_df <- data.frame(num_samples = sp_nonibd$sites, proteins = sp_nonibd$richness, sd = sp_nonibd$sd) 
  
  # plot combined accumulation curve
  
  ggplot() +
    geom_ribbon(data = sp_cd_df, aes(x = num_samples, y = proteins, ymin = proteins - sd, ymax = proteins + sd), fill = '#440154FF', alpha = 1/3) +
    geom_point(data = sp_cd_df, aes(x = num_samples, y = proteins, color = "CD"), size = 1) + 
    geom_ribbon(data = sp_uc_df, aes(x = num_samples, y = proteins, ymin = proteins - sd, ymax = proteins + sd), fill = '#FDE725FF', alpha = 1/3) +
    geom_point(data = sp_uc_df, aes(x = num_samples, y = proteins, color = 'UC'), size = 1) + 
    geom_ribbon(data = sp_nonibd_df, aes(x = num_samples, y = proteins, ymin = proteins - sd, ymax = proteins + sd), fill = '#1F968BFF', alpha = 1/3) +
    geom_point(data = sp_nonibd_df, aes(x = num_samples, y = proteins, color = 'nonIBD'), size = 1) + 
    theme_minimal() +
    #theme(plot.title.position = "plot") +
    scale_color_manual(name = "Diagnosis", 
                       values = c("CD" = '#440154FF', "UC" = '#FDE725FF', "nonIBD" = '#1F968BFF'),
                       labels = c("CD", "UC", "nonIBD"),
                       guide = 'legend') +
    labs(y = "MinHash k-mers")
# }