library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(rjson)

# read in gtdb derep + genbank protozoa, fungi, viral gather results.
# subset to assemblies identified in all models. 
# summarize to species level to group queries of the same species together.

# the output of this script is the representative assemblies as a dummy 
# gather csv output file, to be used as input to genome-grist. Genome-grist will 
# orchestrate the download of these accessions and generate spacegraphcats 
# configuration files. 
# All shared assemblies are output as one csv file, but the assemblies that
# hold >1% of variable importance are the ones will be downloaded by genome-grist

# read in gather results and subset to accessions shared between all models
gather_files <- snakemake@input[["gather"]]
gather_files <- unlist(gather_files, use.names=FALSE)
# gather_files <- Sys.glob("outputs/gather/*gtdb*seed*csv")
gather_results <- gather_files %>%
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

gather_results_shared_assemblies <- gather_results %>%
  filter(name %in% shared_assemblies$name)

# use an arbitrary gather result file, filter to the shared assemblies, 
# and output the file as dummy input to genome-grist download.

gather_dummy_all <- gather_files %>%
  map_dfr(read_csv)  %>%
  filter(name %in% gather_results_shared_assemblies$name) %>%
  group_by(name) %>%
  slice(1) %>%
  ungroup() %>%
  filter(! name %in% c("GCA_002893335.1 Cyclospora cayetanensis strain=CDC:HCTX208:15, ASM289333v1",  
                       "GCA_002893405.1 Cyclospora cayetanensis strain=CDC:HCFL47:13, ASM289340v1",
                       "GCA_003057635.1 Cyclospora cayetanensis strain=CDC:HCTX205:15, ASM305763v1"))

write_csv(gather_dummy_all, snakemake@output[["gather_all_shared"]])


# filter based on variable importance -------------------------------------

# shared assemblies are present in all models; read in smallest set
# reads in signatures of matches and generates a hash:assembly dataframe
# sig_json <- fromJSON(file = "outputs/gather/SRP057027_vita_vars_gtdb_seed6.matches")
sig_json <- fromJSON(file = unlist(snakemake@input[["gather"]], use.names = F)[1])
hash_to_assembly <- lapply(sig_json, data.frame, stringsAsFactors = FALSE)
hash_to_assembly <- bind_rows(hash_to_assembly)

hash_to_assembly_shared_assemblies <- hash_to_assembly %>%
  filter(name %in% shared_assemblies$name) %>%
  select(name, signatures.mins)

# read in and normalize importances from models 
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

#varimp <- Sys.glob("outputs/optimal_rf_seed/*RDS") %>%
varimp <- unlist(snakemake@input[["varimp"]], use.names = F) %>%
  map_dfr(read_varimp) %>%
  mutate(all_model_norm_imp = model_norm_imp / 36)

# join with assembly information 
varimp <- left_join(varimp, hash_to_assembly_shared_assemblies, by = c("hash" ="signatures.mins"))
varimp <- varimp %>%
  mutate(accession = gsub("\\..*", "", name))
# retain only assemblies that have > 1% of variable importance
top_imp_genomes <- name_imp %>%
  filter(total_all_model_norm_imp >= 0.01) 

# this is sort of cheating, but remove the assemblies that I know become <1%
# after running charcoal to remove contamination
gather_dummy_top <- gather_dummy_all %>%
  filter(name %in% top_imp_genomes$name)
  filter(! name %in% c("GCA_002893375.1 Cyclospora cayetanensis strain=CDC:HCVA02:15, ASM289337v1"))

write_csv(gather_dummy_top, snakemake@output[["gather_grist"]])
