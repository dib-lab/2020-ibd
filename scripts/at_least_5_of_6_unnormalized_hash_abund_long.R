library(dplyr)
library(tidyr)
library(readr)

## Get file names
# files <- snakemake@input
# files <- unlist(files, use.names=FALSE)
files <- list.files("outputs/filt_sigs_named_csv", "_filt_named.csv$", full.names = T)

## read files into list with each row labelled by sample.
ibd_long <- list()
for(i in 1:length(files)){
  print(i)
  sig <- read.csv(files[i])                # read in signature csv as df
  sample <- colnames(sig)[2]               # set lib name using sample
  sample <- gsub("_filt_named.sig", "", sample)  # remove file suffix
  sample <- gsub("^X", "", sample)         # remove X appended by R on import
  sig$sample <- sample                     # set libname as col
  colnames(sig) <- c("minhash", "abund", "sample")
  ibd_long[[i]] <- sig
}

## bind into one dataframe
ibd_long <- do.call(rbind, ibd_long)

hash_lib_size <- ibd_long %>%
  group_by(sample) %>%
  summarise(hash_lib_size = sum(abund))

at_least_5_hashes <- read.table("outputs/vita_rf/at_least_5_studies_vita_vars.txt")

ibd_long_filt <- ibd_long %>%
  filter(minhash %in% at_least_5_hashes$V1)

tmp <- pivot_wider(ibd_long_filt, id_cols = "sample", names_from = "minhash", values_from = "abund")

write_tsv(tmp, "outputs/hash_tables/at_least_5_of_6_unnormalized_abund_hashes_wide.tsv")
write_tsv(hash_lib_size, "outputs/hash_tables/hash_lib_sizes.tsv")
