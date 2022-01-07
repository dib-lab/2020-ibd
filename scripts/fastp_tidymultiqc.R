library(purrr)
library(dplyr)
library(TidyMultiqc)

multiqc_fastp <- unlist(snakemake@input) %>%
  set_names() %>%
  map_dfr(load_multiqc, .id = "library_name") %>%
  mutate(library_name = gsub("outputs/sgc_genome_queries_fastp/", "", library_name),
         library_name = gsub("/multiqc_data/multiqc_data.json", "", library_name),
         accession = gsub("_genomic.fna", "", metadata.sample_id))

write_tsv(multiqc_fast, snakemake@output[['tsv']]
