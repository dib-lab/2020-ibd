setwd("~/github/ibd")

library(feather)
library(dplyr)
library(readr)
library(ranger)

# ibd <- read_feather("wide_ibd2.feather")
# ibd <- as.data.frame(ibd)
# rownames(ibd) <- ibd$sample
# ibd <- ibd[ , ncol(ibd)]

# info <- read_tsv("inputs/working_metadata.tsv") %>% 
#   select(study_accession, library_name, diagnosis) %>%
#   mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
#   filter(library_name %in% rownames(ibd)) %>%
#   distinct()
# 
# info_validation <- info %>% 
#   filter(study_accession %in% c("SRP057027", "PRJNA385949"))

# ibd_validation <- ibd[rownames(ibd) %in% info_validation$library_name, ]
# 
# ibd_vita <- readRDS(ibd_vita, "sandbox/greater_than_one_filt_sigs_rf/ibd_vita.RDS")
# 
# var <- ibd_vita$var                        # separate out selected predictive hashes
# var <- gsub("X", "", var)                  # remove the X from the beginning of hashes
# ibd_validation_filt <- ibd_validation[ , colnames(ibd_validation) %in% var] # subset ibd to hashes in ibd_vita


# or read in subsetted ibd validation filter file
ibd_validation_filt <- read.csv("sandbox/rf_filt_vita_output/ibd_validation_filt.csv",
                                row.names = 1)

info_validation <- read_tsv("inputs/working_metadata.tsv") %>% 
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(ibd_validation_filt)) %>%
  distinct()


# order diagnosis vector by order of ibd_validation_filt rownames
info_validation <- info_validation[match(rownames(ibd_validation_filt), info_validation$library_name), ]
# check order...
tail(info_validation)
tail(rownames(ibd_validation_filt))

## add diagnosis vector to validation set
ibd_validation_filt$diagnosis <- info_validation$diagnosis

optimal_ranger <- readRDS("sandbox/rf_filt/optimal_ranger.RDS")

pred_valid <- predict(optimal_ranger, ibd_validation_filt)
table(observed = ibd_validation_filt$diagnosis, 
      predicted = pred_valid$predictions)

pred_valid_df <- data.frame(sample = rownames(ibd_validation_filt),
                            diagnosis = ibd_validation_filt$diagnosis,
                            prediction = pred_valid$predictions)
# predicted
# observed  CD nonIBD  UC
# CD     288    105   6
# nonIBD   5     40   0
# UC       1      2  45
288 + 105 + 6 + 5 + 40 + 1 + 2 + 45

(288+40+45)/492
(105+6+5+1+2+0)/492

info_all <- read_tsv("inputs/working_metadata.tsv") %>% 
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(ibd)) %>%
  distinct()


# separate by cohort ------------------------------------------------------


