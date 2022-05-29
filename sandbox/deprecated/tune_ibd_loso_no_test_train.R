setwd("~/github/2020-ibd/")

library(readr)
library(dplyr)
library(ggplot2)
library(ranger)
library(mlr)
library(tuneRanger)
set.seed(1)

ibd_filt <- read_csv("outputs/vita_rf/PRJNA385949_ibd_filt.csv")
ibd_filt <- as.data.frame(ibd_filt)
rownames(ibd_filt) <- ibd_filt$X1
ibd_filt <- ibd_filt[ , -1]

## read in study metadata
## collapse duplicate libraries so each sample only has one row
info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(ibd_filt)) %>%
  distinct()

## set validation cohort and remove it from variable selection
info_validation <- info %>%
  filter(study_accession == "PRJNA385949") %>%
  mutate(library_name = gsub("-", "\\.", library_name))
ibd_validation <- ibd_filt[rownames(ibd_filt) %in% info_validation$library_name, ]
# match order of to ibd_filt
info_validation <- info_validation[order(match(info_validation$library_name, rownames(ibd_validation))), ]
# check names
all.equal(info_validation$library_name, rownames(ibd_validation))
# make diagnosis var
diagnosis_validation <- info_validation$diagnosis


## remove validation cohort from training data
# using tuneRanger, we do not need to use a train/test/validation framework.
# Instead, tuneRanger does not need a test set because each tree is only trained 
# on a subset of the data (bag), so we can use the rest (out of bag) to obtain 
# an unbiased performance estimation of a single tree and therefore of all trees.
# see: https://github.com/PhilippPro/tuneRanger/issues/8

info_novalidation <- info %>%
  filter(study_accession != "PRJNA385949") %>%
  mutate(library_name = gsub("-", "\\.", library_name))
ibd_novalidation <- ibd_filt[rownames(ibd_filt) %in% info_novalidation$library_name, ]
# match order of to ibd_filt
info_novalidation <- info_novalidation[order(match(info_novalidation$library_name, rownames(ibd_novalidation))), ]
# check names
all.equal(info_novalidation$library_name, rownames(ibd_novalidation))
# make diagnosis var
diagnosis_novalidation <- info_novalidation$diagnosis

# Include classification vars as cols in df
ibd_novalidation$diagnosis <- diagnosis_novalidation
ibd_validation$diagnosis <- diagnosis_validation

# tune ranger -------------------------------------------------------------
# We make an mlr task with the ibd_train dataset here 
tmp <- ibd_novalidation
colnames(tmp) <-  make.names(colnames(tmp))
ibd_task <- makeClassifTask(data = tmp, target = "diagnosis")
# Rough Estimation of the tuning time
# estimateTimeTuneRanger(ibd_task)
# Tuning process (takes 5 hours)
res <- tuneRanger(ibd_task, num.threads = 20)

# # write model parameters to a file
# write_tsv(res$recommended.pars, "tmp.tsv")
# model -------------------------------------------------------------------

# extract model parameters 

# use model parameters to build optimized RF
# ibd_novalidation$diagnosis <- as.factor(ibd_novalidation$diagnosis)
# ibd_rf <- ranger(
#   dependent.variable.name = "diagnosis",
#   mtry            = res$recommended.pars$mtry,
#   num.trees       = 10000,
#   data            = ibd_novalidation,
#   sample.fraction = res$recommended.pars$sample.fraction,
#   min.node.size   = res$recommended.pars$min.node.size,
#   seed            = 1,
#   importance      = 'impurity'
# )

# # from luiz's computer; PREJEB2054
# ibd_novalidation$diagnosis <- as.factor(ibd_novalidation$diagnosis)
# ibd_rf <- ranger(
#   dependent.variable.name = "diagnosis",
#   mtry            = 7346,
#   num.trees       = 10000,
#   data            = ibd_novalidation,
#   sample.fraction = .8665791,
#   min.node.size   = 3,
#   seed            = 1,
#   importance      = 'impurity'
# )

# from luiz's computer; PRJNA237362
# ibd_novalidation$diagnosis <- as.factor(ibd_novalidation$diagnosis)
# ibd_rf <- ranger(
#   dependent.variable.name = "diagnosis",
#   mtry            = 5488,
#   num.trees       = 10000,
#   data            = ibd_novalidation,
#   sample.fraction = .8574862,
#   min.node.size   = 4,
#   seed            = 1,
#   importance      = 'impurity'
# )

# iHMP
# ibd_novalidation$diagnosis <- as.factor(ibd_novalidation$diagnosis)
# ibd_rf <- ranger(
#   dependent.variable.name = "diagnosis",
#   mtry            = 15828,
#   num.trees       = 10000,
#   data            = ibd_novalidation,
#   sample.fraction = 0.8861247,
#   min.node.size   = 2,
#   seed            = 1,
#   importance      = 'impurity'
# )
# saveRDS(ibd_rf, "sandbox/ihmp_as_validation_optimal_RF.RDS")
# SRP057027
# ibd_novalidation$diagnosis <- as.factor(ibd_novalidation$diagnosis)
# ibd_rf <- ranger(
#   dependent.variable.name = "diagnosis",
#   mtry            = 7400,
#   num.trees       = 10000,
#   data            = ibd_novalidation,
#   sample.fraction = 0.7533606,
#   min.node.size   = 2,
#   seed            = 1,
#   importance      = 'impurity'
# )

#PRJNA400072
# ibd_novalidation$diagnosis <- as.factor(ibd_novalidation$diagnosis)
# ibd_rf <- ranger(
#   dependent.variable.name = "diagnosis",
#   mtry            = 772,
#   num.trees       = 10000,
#   data            = ibd_novalidation,
#   sample.fraction = 0.8394355,
#   min.node.size   = 2,
#   seed            = 1,
#   importance      = 'impurity'
# )

# PRJNA385949
ibd_novalidation$diagnosis <- as.factor(ibd_novalidation$diagnosis)
ibd_rf <- ranger(
  dependent.variable.name = "diagnosis",
  mtry            = 8378,
  num.trees       = 10000,
  data            = ibd_novalidation,
  sample.fraction = 0.8267662,
  min.node.size   = 2,
  seed            = 1,
  importance      = 'impurity'
)
source("~/github/2020-ibd/sandbox/loo_rf/function_evaluate_model.R")

evaluate_model(optimal_ranger = ibd_rf, 
               data = ibd_novalidation, 
               reference_class = diagnosis_novalidation, 
               set = "novalidation", 
               study_as_validation = "PRJNA385949")

evaluate_model(optimal_ranger = ibd_rf, 
               data = ibd_validation, 
               reference_class = diagnosis_validation, 
               set = "validation", 
               study_as_validation = "PRJNA385949")


# assess CD as non-IBD assignment -----------------------------------------

# 17 CD predicted to be nonIBD, while 43 CD are predicted correctly. 
# The code below will run a chi-square test, but only 7 of the 17 
# that were misclassified have metadata associated with CD location,
# so there is not enough data to run statistics correctly.

# generate metadata table of CD location.
hmp<- read_tsv("inputs/hmp2_mgx_metadata.tsv") %>%
  filter(External.ID %in% info$library_name) %>%
  select(External.ID, baseline_montreal_location) %>%
  right_join(info, by = c( "External.ID" = "library_name")) %>%
  filter(!is.na(baseline_montreal_location))

# generate a dataframe of predictions
pred_validation <- predict(ibd_rf, ibd_validation)
cd_test <- info_validation
cd_test$prediction <- pred_validation$predictions

# join to montreal classification
cd_test <- left_join(hmp, cd_test, by = c("External.ID" = "library_name" )) %>%
  select(External.ID, baseline_montreal_location, diagnosis.x, prediction)


cd_test %>%
  group_by(diagnosis.x, prediction, baseline_montreal_location) %>%
  tally()


colon <- vector()
for(i in 1:length(cd_test$baseline_montreal_location)){
  iter = cd_test$baseline_montreal_location[i]
  if(iter == "L1") {
    colon[i] <- "non-colon"
  } else if(iter == "L4") {
    colon[i] <- "non-colon"
  } else if(iter == "L1+L4") {
    colon[i] <- "non-colon"
  } else {
    colon[i] <- "colon"
  }
}

cd_test$loc <- colon
cd_test$prediction  <- droplevels(cd_test$prediction)
tbl <- table(cd_test$prediction, cd_test$loc)
tbl
chisq.test(tbl)
