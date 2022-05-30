setwd("~/github/2020-ibd")
library(readr)
library(UpSetR)
prjeb2054 <- read_tsv("outputs/vita_rf/PRJEB2054_vita_vars.txt", col_names = "hash")
prjna237362 <- read_tsv("outputs/vita_rf/PRJNA237362_vita_vars.txt", col_names = "hash")
ihmp <- read_tsv("outputs/vita_rf/iHMP_vita_vars.txt", col_names = "hash")
srp057027 <- read_tsv("outputs/vita_rf/SRP057027_vita_vars.txt", col_names = "hash")
prjna385949 <- read_tsv("outputs/vita_rf/PRJNA385949_vita_vars.txt", col_names = "hash")
prjna400072 <- read_tsv("outputs/vita_rf/PRJNA400072_vita_vars.txt", col_names = "hash")


var_list <- list(PRJEB2054 = prjeb2054$hash, 
                 PRJNA237362 = prjna237362$hash,
                 iHMP = ihmp$hash,
                 SRP057027 = srp057027$hash,
                 PRJNA400072 = prjna400072$hash, 
                 PRJNA385949 = prjna385949$hash)
pdf(file = "new_figure_scripts/upset_var_imp_hashes.pdf", width = 20, height = 5)
upset(fromList(var_list), order.by = "degree", nsets = 6, nintersects = 100,
      number.angles = 30, text.scale = c(2, 1.8, 2, 1.8, 1.5, 1.2))
dev.off()


pdf(file = "new_figure_scripts/upset_var_imp_hashes_small.pdf", width = 8, height = 5)
upset(fromList(var_list), order.by = "degree", nsets = 6, nintersects = 7,
      number.angles = 30, text.scale = c(2, 1.8, 2, 1.35, 1.5, 1.6))
dev.off()

num_total_hashes <- length(unique(c(prjeb2054$hash, prjna237362$hash,
         ihmp$hash, srp057027$hash, 
         prjna400072$hash, prjna385949$hash)))

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

length(unique(c(intersect_hashes_all, intersect_hashes_no_ihmp, intersect_hashes_no_prjeb2054,
              intersect_hashes_no_prjna237362, intersect_hashes_no_prjna385949, 
              intersect_hashes_no_prjna400072, intersect_hashes_no_srp057027)))

at_least_5_hashes <- unique(c(intersect_hashes_all, intersect_hashes_no_ihmp, intersect_hashes_no_prjeb2054,
                              intersect_hashes_no_prjna237362, intersect_hashes_no_prjna385949, 
                              intersect_hashes_no_prjna400072, intersect_hashes_no_srp057027))

write.table(at_least_5_hashes, "outputs/vita_rf/at_least_5_studies_vita_vars.txt",
            quote = F, row.names = F, col.names = F)
