## Get file names
files <- snakemake@input
files <- unlist(files, use.names=FALSE)

## read files into list with each row labelled by sample. 
## normalize files by number of hashes when importing.
ibd_long <- list()
for(i in 1:length(files)){
  print(i)
  if(file.size(files[i]) == 0) {next} else
  {
    sig <- read.csv(files[i])                # read in signature csv as df
    num_hashes <- nrow(sig)                  # get number of hashes
    if(num_hashes == 0) {next} else
    {
      sig[ , 2] <- sig[ , 2] / num_hashes      # normalize by number of hashes
      sample <- colnames(sig)[2]               # set lib name using sample
      sample <- gsub("_singlem_reads\\.fq", "", sample)  # remove file suffix
      sample <- gsub("^outputs\\/sgc_genome_queries_singlem_reads\\/", "", sample)         # remove X appended by R on import
      sig$sample <- sample                     # set libname as col
      colnames(sig) <- c("minhash", "abund", "sample")
      ibd_long[[i]] <- sig
    }
  }
}
## bind into one dataframe
ibd_long <- do.call(rbind, ibd_long)
write.csv(ibd_long, snakemake@output[['csv']],
          quote = F, row.names = F)

