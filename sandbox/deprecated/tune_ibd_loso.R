setwd("~/github/2020-ibd/")

library(readr)
library(dplyr)
library(ggplot2)
library(ranger)
library(mlr)
library(tuneRanger)
set.seed(1)

ibd_filt <- read_csv("outputs/vita_rf/PRJEB2054_ibd_filt.csv")
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
  filter(study_accession == "PRJEB2054") %>%
  mutate(library_name = gsub("-", "\\.", library_name))
ibd_validation <- ibd_filt[rownames(ibd_filt) %in% info_validation$library_name, ]
# match order of to ibd_filt
info_validation <- info_validation[order(match(info_validation$library_name, rownames(ibd_validation))), ]
# check names
all.equal(info_validation$library_name, rownames(ibd_validation))
# make calprotectin var
diagnosis_validation <- info_validation$diagnosis


## remove validation cohort from test/train
info_novalidation <- info %>%
  filter(study_accession != "PRJEB2054") %>%
  mutate(library_name = gsub("-", "\\.", library_name))
ibd_novalidation <- ibd_filt[rownames(ibd_filt) %in% info_novalidation$library_name, ]
# match order of to ibd_filt
info_novalidation <- info_novalidation[order(match(info_novalidation$library_name, rownames(ibd_novalidation))), ]


## set train/test samples
info_train <- info_novalidation


# Random sample indexes
train_index <- sample(1:nrow(info_novalidation), 0.8 * nrow(info_novalidation))
test_index <- setdiff(1:nrow(info_novalidation), train_index)

# Build train set
info_train <- info_novalidation[train_index, ]
# build ibd_train
ibd_train <- ibd_filt[rownames(ibd_filt) %in% info_train$library_name, ]
# match order of to ibd_filt
info_train <- info_train[order(match(info_train$library_name, rownames(ibd_train))), ]
# check names
all.equal(info_train$library_name, rownames(ibd_train))
# make calprotectin var
diagnosis_train <- info_train$diagnosis

# Build test set
info_test <- info_novalidation[test_index, ]
# build ibd_test
ibd_test<- ibd_filt[rownames(ibd_filt) %in% info_test$library_name, ]
# match order of to ibd_filt
info_test <- info_test[order(match(info_test$library_name, rownames(ibd_test))), ]
# check names
all.equal(info_test$library_name, rownames(ibd_test))
# make calprotectin var
diagnosis_test <- info_test$diagnosis


# Include classification vars as cols in df
ibd_train$diagnosis <- diagnosis_train
ibd_test$diagnosis <- diagnosis_test
ibd_validation$diagnosis <- diagnosis_validation

# tune ranger -------------------------------------------------------------
# We make an mlr task with the ibd_train dataset here 
tmp <- ibd_train
colnames(tmp) <-  make.names(colnames(tmp))
ibd_task <- makeClassifTask(data = tmp, target = "diagnosis")
# Rough Estimation of the tuning time
# estimateTimeTuneRanger(ibd_task)
# Tuning process (takes 5 hours)
res <- tuneRanger(ibd_task, num.threads = 3)

# write model parameters to a file
write_tsv(res$recommended.pars, "tmp.tsv")
# model -------------------------------------------------------------------

# extract model parameters 

# use model parameters to build optimized RF
ibd_train$diagnosis <- as.factor(ibd_train$diagnosis)
ibd_rf <- ranger(
  dependent.variable.name = "diagnosis",
  mtry            = res$recommended.pars$mtry,
  num.trees       = 10000,
  data            = ibd_train,
  sample.fraction = res$recommended.pars$sample.fraction,
  min.node.size   = res$recommended.pars$min.node.size,
  seed            = 1,
  importance      = 'impurity'
)

# from luiz's computer
# ibd_train$diagnosis <- as.factor(ibd_train$diagnosis)
# ibd_rf <- ranger(
#   dependent.variable.name = "diagnosis", 
#   mtry            = 4238, 
#   num.trees       = 10000,
#   data            = ibd_train,
#   sample.fraction = .82412,
#   min.node.size   = 6,
#   seed            = 1,
#   importance      = 'impurity'
# )


evaluate_model(optimal_ranger = ibd_rf, 
               data = ibd_train, 
               reference_class = diagnosis_train, 
               set = "train", 
               study_as_validation = "PRJEB2054")

evaluate_model(optimal_ranger = ibd_rf, 
               data = ibd_test, 
               reference_class = diagnosis_test, 
               set = "test", 
               study_as_validation = "PRJEB2054")

evaluate_model(optimal_ranger = ibd_rf, 
               data = ibd_validation, 
               reference_class = diagnosis_validation, 
               set = "validation", 
               study_as_validation = "PRJEB2054")
