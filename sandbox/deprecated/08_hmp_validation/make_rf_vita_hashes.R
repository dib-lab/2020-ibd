setwd("~/github/ibd")

ibd_novalidation_filt <- read.csv("sandbox/rf_filt/ibd_novalidation_filt.csv",
                                  row.names = 1)
hashes <- colnames(ibd_novalidation_filt)
hashes <- gsub("X", "", hashes)
write.table(hashes, "sandbox/hmp_validation/rf_vita_hashes.txt", quote = F, 
            row.names = F, col.names = F)
