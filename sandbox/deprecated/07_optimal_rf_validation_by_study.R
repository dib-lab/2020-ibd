setwd("~/github/ibd")

library(dplyr)
library(readr)
library(ranger)

set.seed(1)

# validation by study -----------------------------------------------------

optimal_ranger <- readRDS("outputs/optimal_rf/optimal_ranger.RDS")

# or read in subsetted ibd validation filter file
ibd_validation_filt <- read.csv("sandbox/rf_filt_vita_output/ibd_validation_filt.csv",
                                row.names = 1)

info_validation <- read_tsv("inputs/working_metadata.tsv") %>% 
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(ibd_validation_filt)) %>%
  distinct()


# SRP057027 ----------------------------------------------------------------

info_srp <- filter(info_validation, study_accession == "SRP057027")
ibd_srp <- ibd_validation_filt[rownames(ibd_validation_filt) %in% info_srp$library_name, ]
info_srp <- info_srp[match(rownames(ibd_srp), info_srp$library_name), ]
ibd_srp$diagnosis <- info_srp$diagnosis

pred_srp <- predict(optimal_ranger, ibd_srp)
pred_srp_tab <- table(observed = ibd_srp$diagnosis, predicted = pred_srp$predictions)
write.table(pred_srp_tab, "outputs/rf_validation/pred_srp057027.txt")

# observed  CD nonIBD  UC
# CD       251     56   5
# nonIBD     5     20   0
# (251+20)/nrow(ibd_srp)
# (56+5+5) / nrow(ibd_srp)

pred_srp_df <- data.frame(sample = rownames(ibd_srp),
                            diagnosis = ibd_srp$diagnosis,
                            prediction = pred_srp$predictions)

write.csv(pred_srp_df, "outputs/rf_validation/pred_srp057027.csv",
          quote = F, row.names = F)
# PRJNA285949 -------------------------------------------------------------

info_prjna <- filter(info_validation, study_accession == "PRJNA385949")
ibd_prjna <- ibd_validation_filt[rownames(ibd_validation_filt) %in% info_prjna$library_name, ]
info_prjna <- info_prjna[match(rownames(ibd_prjna), info_prjna$library_name), ]
ibd_prjna$diagnosis <- info_prjna$diagnosis

pred_prjna <- predict(optimal_ranger, ibd_prjna)
pred_prjna_tab <- table(observed = ibd_prjna$diagnosis, predicted = pred_prjna$predictions)
write.table(pred_prjna_tab, "outputs/rf_validation/pred_prjna285949.txt")

# observed CD nonIBD UC
# CD       37     49  1
# nonIBD    0     20  0
# UC        1      2 45
# (37+20+45)/nrow(ibd_prjna)
# (49+1+2+1) / nrow(ibd_prjna)

pred_prjna_df <- data.frame(sample = rownames(ibd_prjna),
                            diagnosis = ibd_prjna$diagnosis,
                            prediction = pred_prjna$predictions)
write.csv(pred_prjna_df, "outputs/rf_validation/pred_prjna285949.csv",
          quote = F, row.names = F)

