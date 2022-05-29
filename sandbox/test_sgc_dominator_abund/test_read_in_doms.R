library(readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
setwd("~/github/2020-ibd")

domset_abund_files <- list.files("sandbox/test_sgc_dominator_abund/sgc_nbhd_rgnv/rgnv_count_dom_out", 
                           "_domset_abund.txt$", full.names = T)

domset_abund <- domset_abund_files %>%
  set_names() %>%
  map_dfr(read_csv, .id = "sample")
  
domset_abund <- domset_abund %>%
  pivot_wider(id_cols = sample, names_from = dom_id, values_from = abund) 

domset_abund %>%
  mutate(sample = gsub("sandbox/test_sgc_dominator_abund/sgc_nbhd_rgnv/rgnv_count_dom_out/", "", sample)) %>%
  mutate(sample = gsub("_rgnv_nbhd_domset_abund.txt", "", sample))
head(domset_abund)



# try cbind approach ------------------------------------------------------

# It's a bad idea to import hte files this way because there is no check that 
# dom_ids are the same between all samples other than that they will have the 
# number of rows. However, this is much faster and will probably be the
# solution that I use to because it will scale to hundreds or thousands of 
# samples with large numbers of dominating sets.

domset_abund_files <- list.files("sandbox/test_sgc_dominator_abund/sgc_nbhd_rgnv/rgnv_count_dom_out", 
                                 "_domset_abund.txt$", full.names = T)
domnames <- read_csv(domset_abund_files[1])[ , "dom_id"]     # read in domset id
df <- do.call(cbind, lapply(domset_abund_files, 
                            function(x) read_csv(x)[ , "abund"]))
df <- cbind(domnames, df)
colnames(df) <- c("dom_id", domset_abund_files)
colnames(df) <- gsub("sandbox/test_sgc_dominator_abund/sgc_nbhd_rgnv/rgnv_count_dom_out/", "", colnames(df))
colnames(df) <- gsub("_rgnv_nbhd_domset_abund.txt", "", colnames(df))

tmp <- as.matrix(df[ , -1])
hist(tmp)



# on the cluster ----------------------------------------------------------

domset_abund_files <- list.files("rgnv_count_dom_out", "_domset_abund.txt$", full.names = T)
domnames <- read_csv(domset_abund_files[1])[ , "dom_id"]     # read in domset id
df <- do.call(cbind, lapply(domset_abund_files, 
                            function(x) read_csv(x)[ , "abund"]))
df <- cbind(domnames, df)
colnames(df) <- c("dom_id", domset_abund_files)
colnames(df) <- gsub("rgnv_count_dom_out/", "", colnames(df))
colnames(df) <- gsub("_rgnv_nbhd_domset_abund.txt", "", colnames(df))
# save the output
write_tsv(df, "rgnv_count_domset_abund.tsv")
saveRDS(df, "rgnv_count_domset_abund.RDS")

# read in when crashes
df <- read_tsv("rgnv_count_domset_abund.tsv")

# make a presence absence matrix
pa <- df %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

# determine how many samples have a domset abund of 1 or more
pa %>%
  mutate(total = rowSums(across(where(is.numeric))))
total <- rowSums(pa[ , -1])

sorted_total <- sort(total)