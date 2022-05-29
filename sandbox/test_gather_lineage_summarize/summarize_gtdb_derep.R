library(dplyr)
library(readr)
library(tidyr)
library(purrr)

setwd("~/github/2020-ibd/")

# read in gtdb derep + genbank protozoa, fungi, viral gather results.
# subset to assemblies identified in all models. 
# summarize to species level to group queries of the same species together.

# read in gather results and subset to accessions shared between all models
gather_results <- Sys.glob("outputs/gather/*gtdb*seed*csv") %>%
  set_names() %>%
  map_dfr(read_csv, .id = "source")  %>%
  mutate(source = gsub("outputs/gather/", "", source)) %>%
  mutate(source = gsub("vita_vars_gtdb_", "", source)) %>%
  mutate(source = gsub(".csv", "", source)) %>%
  separate(col = name, into = c("accession"), remove = F, sep = " ", extra = "drop") %>%
  separate(col = source, into = c("study", "seed"), remove = T, sep = "_", extra = "drop") %>%
  mutate(accession = gsub("\\..*", "", accession))

shared_assemblies <- gather_results %>%
  group_by(name) %>%
  tally() %>%
  filter(n == 36)

dim(shared_assemblies) # 364 shared assemblies

gather_results_shared_assemblies <- gather_results %>%
  filter(name %in% shared_assemblies$name)


# read in summarized lineages from gather results ----------------------------
lineages <- read_csv("sandbox/test_gather_lineage_summarize/gtdb-rs202-genbank-protozoa-viral-fungi-lineage.csv")

gather_results_shared_assemblies <- gather_results_shared_assemblies %>%
  left_join(lineages, by = c("accession" ="ident")) %>%
  mutate(lineage = paste(superkingdom, phylum, class, order, family, genus, species, sep = ";"))


# replace the "unassigneds" with NAs to match formating with lineages DF
# gather_results_shared_assemblies[gather_results_shared_assemblies == "unassigned" ] <- NA

length(unique(gather_results_shared_assemblies$name))
length(unique(gather_results_shared_assemblies$lineage))
View(gather_results_shared_assemblies)

tmp <- gather_results_shared_assemblies %>%
  select(-intersect_bp, -f_orig_query, -f_match, -f_unique_to_query, 
         -f_unique_weighted, -average_abund, -median_abund, 
         -std_abund, -filename, -md5, -f_match_orig, -unique_intersect_bp, 
         -gather_result_rank, -remaining_bp, -X2) %>%
  group_by(lineage) %>%
  tally()

# what fraction of hashes from each model is classifiable with ncbi --------
gather_frac <- gather_results %>%
  group_by(study, seed) %>%
  summarize(total_original_query = sum(f_unique_to_query))


# investigate eukaryotic fraction -----------------------------------------


View(gather_results_shared_assemblies %>%
  filter(superkingdom == "Eukaryota"))

gather_results_shared_assemblies %>%
  filter(superkingdom == "Eukaryota") %>%
  group_by(name) %>%
  summarize(sum_f_unique_to_query = sum(f_unique_to_query))

# name                                                                       sum_f_unique_to_query
# <chr>                                                                                      <dbl>
# 1 GCA_000256725.2 Toxoplasma gondii TgCatPRC2 strain=TgCatPRC2, TGCATPRC2 v2               0.0101 
# 2 GCA_002754635.1 Plasmodium vivax strain=CMB-1, CMB-1_v2                                  0.0302 
# 3 GCA_002893335.1 Cyclospora cayetanensis strain=CDC:HCTX208:15, ASM289333v1               0.0376 
# 4 GCA_002893375.1 Cyclospora cayetanensis strain=CDC:HCVA02:15, ASM289337v1                0.638  
# 5 GCA_002893405.1 Cyclospora cayetanensis strain=CDC:HCFL47:13, ASM289340v1                0.0619 
# 6 GCA_003057635.1 Cyclospora cayetanensis strain=CDC:HCTX205:15, ASM305763v1               0.00739
# 7 GCA_900088545.1 Plasmodium ovale wallikeri, Pow2                                         0.00817

# Choose GCA_002893375.1 Cyclospora cayetanensis strain=CDC:HCVA02:15, ASM289337v1 as representative