library(dplyr)
hmp <- read_csv("~/github/cosmo-kmers/inputs/hmp2_metadata.csv")
colnames(hmp) <- make.names(colnames(hmp))
hmp <- hmp %>%
  filter(data_type == "metatranscriptomics") %>%
  filter(!grepl("_P", External.ID)) %>%
  filter(!grepl("PSM6XBTP_TR", External.ID))

write_tsv(hmp, "~/Downloads/hmp2_metadata_metatranscriptomics.tsv")
