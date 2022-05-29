library(readr)
library(dplyr)

hashes <- read_tsv("outputs/vita_rf/at_least_5_studies_vita_vars.txt", col_names = "hash")

ihmp <- read_csv("outputs/vita_rf/iHMP_ibd_filt.csv")
ihmp_filt <- ihmp %>%
  select(X1, colnames(ihmp)[colnames(ihmp) %in% as.character(hashes$hash)])

prjeb2054 <- read_csv("outputs/vita_rf/PRJEB2054_ibd_filt.csv")
prjeb2054_filt <- prjeb2054 %>%
  select(X1, colnames(prjeb2054)[colnames(prjeb2054) %in% as.character(hashes$hash)])

PRJNA237362 <- read_csv("outputs/vita_rf/PRJNA237362_ibd_filt.csv")
PRJNA237362_filt <- PRJNA237362 %>%
  select(X1, colnames(PRJNA237362)[colnames(PRJNA237362) %in% as.character(hashes$hash)])

PRJNA385949<- read_csv("outputs/vita_rf/PRJNA385949_ibd_filt.csv")
PRJNA385949_filt <- PRJNA385949  %>%
  select(X1, colnames(PRJNA385949)[colnames(PRJNA385949) %in% as.character(hashes$hash)])

PRJNA400072 <- read_csv("outputs/vita_rf/PRJNA400072_ibd_filt.csv")
PRJNA400072_filt <- PRJNA400072  %>%
  select(X1, colnames(PRJNA400072)[colnames(PRJNA400072) %in% as.character(hashes$hash)])

SRP057027 <- read_csv("outputs/vita_rf/SRP057027_ibd_filt.csv")
SRP057027_filt <- SRP057027  %>%
  select(X1, colnames(SRP057027)[colnames(SRP057027) %in% as.character(hashes$hash)])

all <- full_join(ihmp_filt, prjeb2054_filt)
all <- full_join(all, PRJNA237362_filt)
all <- full_join(all, PRJNA385949_filt)
all <- full_join(all, PRJNA400072_filt)
all <- full_join(all, SRP057027_filt)
dim(all)

write_tsv(all, "sandbox/at_least_5_studies_vita_vars.tsv")
