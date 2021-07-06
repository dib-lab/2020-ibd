library(readr)
library(dplyr)
library(stringr)

singlem <- read_tsv(snakemake@input[['singlem']])
reads <- unlist(str_split(string = singlem$read_names, pattern = " "))
write.table(reads, snakemake@output[['reads']], quote = F, col.names = F, row.names = F)
