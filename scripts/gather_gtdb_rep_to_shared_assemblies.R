library(dplyr)
library(readr)
library(tidyr)
library(purrr)

# read in gtdb derep + genbank protozoa, fungi, viral gather results.
# subset to assemblies identified in all models. 
# summarize to species level to group queries of the same species together.

# the two output of this script are
# 1. representative assemblies summarized to lineage
# 2. representative assemblies as a dummy gather csv output file, to be used
#    as input to genome-grist. Genome-grist will orchestrate the download of 
#    these accessions and generate spacegraphcats configuration files. 

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

gather_dummy <- gather_files %>%
  map_dfr(read_csv)  %>%
  filter(name %in% gather_results_shared_assemblies$name) %>%
  group_by(name) %>%
  slice(1) %>%
  ungroup() %>%
  filter(! name %in% c("GCA_002893335.1 Cyclospora cayetanensis strain=CDC:HCTX208:15, ASM289333v1",  
                       "GCA_002893405.1 Cyclospora cayetanensis strain=CDC:HCFL47:13, ASM289340v1",
                       "GCA_003057635.1 Cyclospora cayetanensis strain=CDC:HCTX205:15, ASM305763v1"))

# select Eukaryotic representative species for 

write_csv(gather_dummy, snakemake@outputs[["gather_grist"]])