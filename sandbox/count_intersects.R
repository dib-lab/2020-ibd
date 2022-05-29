setwd("github/2020-ibd")
library(readr)
library(UpSetR)
prjeb2054 <- read_tsv("outputs/vita_rf/PRJEB2054_vita_vars.txt", col_names = "hash")
prjna237362 <- read_tsv("outputs/vita_rf/PRJNA237362_vita_vars.txt", col_names = "hash")
ihmp <- read_tsv("outputs/vita_rf/iHMP_vita_vars.txt", col_names = "hash")
srp057027 <- read_tsv("outputs/vita_rf/SRP057027_vita_vars.txt", col_names = "hash")
prjna385949 <- read_tsv("outputs/vita_rf/PRJNA385949_vita_vars.txt", col_names = "hash")
prjna400072 <- read_tsv("outputs/vita_rf/PRJNA400072_vita_vars.txt", col_names = "hash")

length(intersect(intersect(intersect(prjeb2054$hash, prjna237362$hash), ihmp$hash), srp057027$hash))

var_list <- list(prjeb2054 = prjeb2054$hash, 
                 prjna237362 = prjna237362$hash,
                 ihmp = ihmp$hash,
                 srp057027 = srp057027$hash,
                 prjna400072 = prjna400072$hash, 
                 prjna385949 = prjna385949$hash)
upset(fromList(var_list), order.by = "degree", nsets = 6, nintersects = 100)
