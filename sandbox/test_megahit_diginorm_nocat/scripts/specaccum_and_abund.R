library(readr)
library(dplyr)
library(ggplot2)
library(vegan)

info <- read_tsv(snakemake@input[['info']]) %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct()

generate_counts <- function(path){
  # read in an format counts
  counts <- read_tsv(path)
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$protein
  counts <- counts[ , -1]
}

generate_presence_absence <- function(path, threshold = 0){
  counts <- generate_counts(path)
  # replace any count > 0 with 1
  counts_pa <- ifelse(counts > threshold, 1, 0)
}

plot_gene_abund <- function(counts_pa, pangenome_name){
  n_proteins_observed_per_sample <- data.frame(library_name = colnames(counts_pa), 
                                               n_prots = colSums(counts_pa))
  
  n_proteins_observed_per_sample <- full_join(n_proteins_observed_per_sample, 
                                              info)
  
  plt <- ggplot(n_proteins_observed_per_sample, aes(x = diagnosis, y = n_prots)) +
    geom_boxplot() +
    theme_minimal() +
    #theme(plot.title.position = "plot") +
    labs(title = pangenome_name, y = "Number of proteins per sample")
  
  tukey <- TukeyHSD(aov(n_prots~diagnosis, data = n_proteins_observed_per_sample))
  return(list(plt = plt, tukey = tukey))
}

plot_species_accumulation <- function(info, counts, pangenome_name){
  info <- info[order(match(info$library_name, colnames(counts))), ]
  stopifnot(all.equal(info$library_name, colnames(counts)))
  
  uc <- filter(info, diagnosis == "UC")
  sp_uc <- specaccum(t(counts[ , colnames(counts) %in% uc$library_name]), "exact")
  sp_uc_df <- data.frame(num_samples = sp_uc$sites, proteins = sp_uc$richness, sd = sp_uc$sd) 
  
  cd <- filter(info, diagnosis == "CD")
  sp_cd <- specaccum(t(counts[ , colnames(counts) %in% cd$library_name]), "exact")
  sp_cd_df <- data.frame(num_samples = sp_cd$sites, proteins = sp_cd$richness, sd = sp_cd$sd)
  
  nonIBD <- filter(info, diagnosis == "nonIBD")
  sp_nonibd <- specaccum(t(counts[ , colnames(counts) %in% nonIBD$library_name]), "exact")
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
    labs(title = pangenome_name)
}

do_all <- function(path, pangenome_name, threshold, info = info, 
                   out_gene_abund_plt,  out_gene_abund_tukey,
                   out_specaccum_plt){
  counts <- generate_counts(path)
  counts_pa <- generate_presence_absence(path, threshold = threshold)
  counts_plt <- plot_gene_abund(counts_pa, pangenome_name = pangenome_name)
  ggsave(out_gene_abund_plt, counts_plt$plt)
  write_tsv(as.data.frame(counts_plt$tukey$diagnosis), path = out_gene_abund_tukey)
  specaccum_plt <- plot_species_accumulation(info = info, counts = counts, 
                                         pangenome_name = pangenome_name)
  
  ggsave(out_specaccum_plt, specaccum_plt)
}


do_all(path = snakemake@input[['counts']], info = info,
       pangenome_name = snakemake@params[['gather_genome']], threshold = 0,
       out_gene_abund_plt = snakemake@output[['gene_abund_plt']],
       out_gene_abund_tukey = snakemake@output[['gene_abund_tukey']],
       out_specaccum_plt = snakemake@output[['specaccum_plt']])
