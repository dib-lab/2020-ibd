library(readr)

prjna237362 <- read_tsv(snakemake@input[[1]], col_names = "hash")
prjna400072 <- read_tsv(snakemake@input[[2]], col_names = "hash")
prjeb2054 <- read_tsv(snakemake@input[[3]], col_names = "hash")
srp057027 <- read_tsv(snakemake@input[[4]], col_names = "hash")
prjna385949 <- read_tsv(snakemake@input[[5]], col_names = "hash")
ihmp <- read_tsv(snakemake@input[[6]], col_names = "hash")

intersect_hashes_all <- Reduce(intersect, list(prjeb2054$hash, prjna237362$hash,
                                               ihmp$hash, srp057027$hash, 
                                               prjna400072$hash, prjna385949$hash))

intersect_hashes_no_prjeb2054 <- Reduce(intersect, list(prjna237362$hash,
                                                        ihmp$hash, srp057027$hash, 
                                                        prjna400072$hash, prjna385949$hash))

intersect_hashes_no_prjna237362 <- Reduce(intersect, list(prjeb2054$hash,
                                                          ihmp$hash, srp057027$hash, 
                                                          prjna400072$hash, prjna385949$hash))

intersect_hashes_no_ihmp <- Reduce(intersect, list(prjeb2054$hash, prjna237362$hash,
                                                   srp057027$hash, 
                                                   prjna400072$hash, prjna385949$hash))

intersect_hashes_no_srp057027 <- Reduce(intersect, list(prjeb2054$hash, prjna237362$hash,
                                                        ihmp$hash,
                                                        prjna400072$hash, prjna385949$hash))

intersect_hashes_no_prjna400072 <- Reduce(intersect, list(prjeb2054$hash, prjna237362$hash,
                                                          ihmp$hash, srp057027$hash, 
                                                          prjna385949$hash))

intersect_hashes_no_prjna385949 <- Reduce(intersect, list(prjeb2054$hash, prjna237362$hash,
                                                          ihmp$hash, srp057027$hash, 
                                                          prjna400072$hash))

at_least_5_hashes <- unique(c(intersect_hashes_all, intersect_hashes_no_ihmp, intersect_hashes_no_prjeb2054,
                              intersect_hashes_no_prjna237362, intersect_hashes_no_prjna385949, 
                              intersect_hashes_no_prjna400072, intersect_hashes_no_srp057027))

write.table(at_least_5_hashes, snakemake@output[['at_least_5']],
            quote = F, row.names = F, col.names = F)
