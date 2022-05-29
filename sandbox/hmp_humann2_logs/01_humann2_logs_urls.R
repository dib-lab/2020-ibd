# Generate download links for HUMANn2 logs

setwd("~/github/ibd/")

library(dplyr)
library(readr)

info <- read_tsv("inputs/hmp2_mgx_metadata.tsv") 

urls <- paste0("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/",
               info$External.ID, 
               "_humann2.tar.bz2")

write.table(urls, file = "sandbox/hmp_humann2_logs//urls.txt",
            col.names = F, row.names = F, quote = F)

# wget -i urls.txt
# for infile in *bz2
# do
# tar xjf $infile
# done
# 
# for infile in *log
# do
# grep "overall alignmnet rate" $infile >> OAR.txt
# done
