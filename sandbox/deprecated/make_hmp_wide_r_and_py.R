
# R -----------------------------------------------------------------------

## Get file names
files <- list.files(path = "outputs/filt_sigs_named_csv_hmp", 
                    pattern = "_filt_named.csv$", 
                    full.names = T)

## read files into list with each row labelled by sample.
## normalize files by number of hashes when importing.
ibd_long <- list()
for(i in 1:length(files)){
  print(i)
  sig <- read.csv(files[i])                # read in signature csv as df
  sample <- colnames(sig)[2]               # set lib name using sample
  sample <- gsub("_filt.sig", "", sample)  # remove file suffix
  sample <- gsub("^X", "", sample)         # remove X appended by R on import
  sig$sample <- sample                     # set libname as col
  colnames(sig) <- c("minhash", "abund", "sample")
  ibd_long[[i]] <- sig
}

## bind into one dataframe
ibd_long <- do.call(rbind, ibd_long)
write.csv(ibd_long, "outputs/hash_tables/hmp_unnormalized_abund_hashes_long.csv",
          quote = F, row.names = F)


# python ------------------------------------------------------------------

import pandas as pd
import feather

input = "outputs/hash_tables/hmp_unnormalized_abund_hashes_long.csv"
output = "outputs/hash_tables/hmp_unnormalized_abund_hashes_wide.feather"

ibd = pd.read_csv(str(input), dtype = {"minhash" : "int64", "abund" : "float64", "sample" : "object"})
ibd_wide=ibd.pivot(index='sample', columns='minhash', values='abund')
ibd_wide = ibd_wide.fillna(0)
ibd_wide['sample'] = ibd_wide.index
ibd_wide = ibd_wide.reset_index(drop=True)
ibd_wide.columns = ibd_wide.columns.astype(str)
ibd_wide.to_feather(str(output)) 
