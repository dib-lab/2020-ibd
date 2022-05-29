setwd("~/github/ibd/")

library(dplyr)
library(data.table)
library(readr)

## Get file names and metadata
files <- list.files("sandbox/greater_than_one_filt_sig_csvs/", 
                    pattern = ".csv$", 
                    full.names = T)
info <- read_tsv("inputs/working_metadata.tsv")

## read files into list with each row labelled by sample. 
## normalize files by number of hashes when importing.
ibd_long <- list()
for(i in 1:length(files)){
  print(i)
  sig <- read.csv(files[i])                # read in signature csv as df
  num_hashes <- nrow(sig)                  # get number of hashes
  sig[ , 2] <- sig[ , 2] / num_hashes      # normalize by number of hashes
  sample <- colnames(sig)[2]               # set lib name using sample
  sample <- gsub("_filt.sig", "", sample)  # remove file suffix
  sample <- gsub("^X", "", sample)         # remove X appended by R on import
  sig$sample <- sample                     # set libname as col
  colnames(sig) <- c("minhash", "abund", "sample")
  ibd_long[[i]] <- sig
}

# bind into one dataframe
ibd_long <- do.call(rbind, ibd_long)
write.csv(ibd_long, "sandbox/greater_than_one_filt_sig_csvs/long_ibd.csv",
          quote = F, row.names = F)
