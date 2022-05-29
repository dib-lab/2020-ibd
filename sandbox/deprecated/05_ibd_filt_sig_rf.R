setwd("~/github/ibd/")

# remotes::install_github("silkeszy/Pomona")
library(Pomona)
library(ranger)
library(dplyr)
library(randomForest)
library(caret)
library(ggplot2)
library(scales)
library(tidyr)

library(data.table)
# library(feather)
# library(readr)

set.seed(1)

# perform variable selection within random forests on ibd minhashes (k 31, 
# scaled 2000, with hashes that were present only once across all samples
# removed). Uses the Pomona package vita implementation, which wraps the ranger 
# package. Saves output as RDS for faster loading of output objects into 
# subsequent R sessions. 


## read in data ------------------------------------------------------------

## format hash table (samples x features)
ibd <- fread("sandbox/greater_than_one_filt_sig_csvs/wide_ibd.csv.gz")
rownames(ibd) <- ibd[ , 1] # set rownames as sample names
ibd <- ibd[ , -1]          # remove sample names from df

## format classification vector
info <- read_csv("inputs/working_metadata.tsv") %>% 
  filter(library_accession %in% rownames(ibd))  %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct() %>%
  

## remove validation set from info
info <- info[match(rownames(ibd), info$library_name), ] # match order of classification to ibd
diagnosis <- info$diagnosis # make diagnosis var

## remove validation data set (SRP...)

## split to test and train
train <- sample(nrow(ibd), 0.7*nrow(ibd), replace = FALSE)
train_set <- ibd[train, ]
test_set <- ibd[-train, ]

diagnosis_train <- diagnosis[train]
diagnosis_test <- diagnosis[-train]

# run vita ----------------------------------------------------------------

## perform variant selection
## var.sel.vita calculates p-values based on the empirical null distribution 
## from non-positive VIMs as described in Janitza et al. (2015). 
ibd_vita <- var.sel.vita(x = train_set, y = diagnosis_train, p.t = 0.05, 
                         ntree = 5000, mtry.prop = 0.2, nodesize.prop = 0.1, 
                         no.threads = 10, method = "ranger", 
                         type = "classification")
saveRDS(ibd_vita, "sandbox/greater_than_one_filt_sigs_rf/ibd_vita.RDS")


# build rf with sel feat ---------------------------------------------------

ibd_vita <- readRDS(ibd_vita, "sandbox/greater_than_one_filt_sigs_rf/ibd_vita.RDS")

var <- ibd_vita$var                        # separate out selected predictive hashes
var <- gsub("X", "", var)                  # remove the X from the beginning of hashes
ibd_filt <- ibd[ , colnames(ibd) %in% var] # subset ibd to hashes in ibd_vita


## Run random forest classification using the hashes from variable selection RF
## split the dataset into train and validation set in the ratio 70:30 
## Because seed was set at the beginning of the script, the same samples should 
## be in this training dataset as were in the var sel rf dataset. 

train <- sample(nrow(ibd_filt), 0.7*nrow(ibd_filt), replace = FALSE)
train_set <- ibd_filt[train, ]
test_set <- ibd_filt[-train, ]
diagnosis_train <- diagnosis[train]
diagnosis_test <- diagnosis[-train]

ibd_filt_rf <- randomForest(x = train_set, y = as.factor(diagnosis_train), 
                            importance = TRUE)
ibd_filt_rf
saveRDS(ibd_filt_rf, "sandbox/greater_than_one_filt_sigs_rf/ibd_filt_rf.RDS")

# work with model ---------------------------------------------------------

ibd_filt_rf <- readRDS("sandbox/ibd_sig_csv_10k_rf/ibd_filt_rf.RDS")

pred_train <- predict(ibd_filt_rf, train_set)
table(observed = diagnosis_train, predicted = pred_train)

# predicted
# observed  CD nonIBD  UC
# CD     238      0   0
# nonIBD   1    138   0
# UC       0      0 152

pred_test <- predict(ibd_filt_rf, test_set)
table(observed = diagnosis_test, predicted = pred_test) 
pred_test <- factor(pred_test, levels = c("nonIBD", "CD", "UC"))
diagnosis_test <- factor(diagnosis_test, levels = c("nonIBD", "CD", "UC"))
cm_test <- caret::confusionMatrix(data = pred_test, 
                                  reference = as.factor(diagnosis_test))
# predicted
# observed  CD nonIBD  UC
# CD     109      0   0
# nonIBD   3     49   1
# UC       7      1  57

varImpPlot(ibd_filt_rf)

# subset ibd_filt to top identifiers; spread to long format; make boxplots
# by disease state

imp <- ibd_filt_rf$importance %>%
  as.data.frame() %>%
  mutate(hash = rownames(ibd_filt_rf$importance)) %>%
  arrange(desc(MeanDecreaseGini)) 

top_hash <- imp$hash[1:10]
ibd_top_hash <- ibd_filt[ , colnames(ibd_filt) %in% top_hash]

ibd_top_hash$sample <- rownames(ibd_top_hash)
top_long <- gather(ibd_top_hash, key = "hash", value = "abund", -sample)
top_long <- left_join(top_long, info, by = c("sample" = "External.ID"))
top_long$diagnosis <- factor(top_long$diagnosis, levels = c("nonIBD", "CD", "UC"))
ggplot(top_long, aes(x = diagnosis, y = abund)) +
  #geom_boxplot() +
  geom_jitter(alpha = 1/10) +
  facet_wrap(~hash, scales = "free") +
  theme_minimal()

ggplotConfusionMatrix <- function(m){
  mycaption <- paste("Accuracy", percent_format()(m$overall[1]),
                     "Kappa", percent_format()(m$overall[2]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
    theme_minimal() +
    theme(legend.position = "none", 
          text = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    labs(caption = mycaption, title = "Random Forest Accuracy")
  return(p)
}

ggplotConfusionMatrix(cm_test)

ggplotConfusionMatrix(cfm)

# do vita hashes cluster samples by disease better than full hashes? NO! ------

library(vegan)
dist <- vegdist(ibd_filt, method="bray")
pcoa <- cmdscale(dist, eig=TRUE)

ibd_dist_all <- as.data.frame(pcoa$points)
ibd_dist_all$sample <- rownames(ibd)
ibd_dist_all <- left_join(ibd_dist_all, info, by = c("sample" = "External.ID"))

# percent variance explained:
var <- round(pcoa$eig*100/sum(pcoa$eig), 1)

ibd_pcoa_plt <- ggplot(ibd_dist_all, 
                       aes(x = V1, y = V2, color = diagnosis)) +
  geom_point() +
  theme_minimal() +
  labs(x = paste0("PCo 1 (", var[1], "%)"), 
       y = paste0("PCo 2 (", var[2], "%)"),
       title = "Metatranscriptome MinHashes")  +
  scale_color_viridis_d() +
  theme(legend.position = "none")

ibd_pcoa_plt <- ggExtra::ggMarginal(ibd_pcoa_plt, type = "density",
                                    groupColour = T, groupFill = T)
ibd_pcoa_plt



# make a two level classifier ---------------------------------------------

diagnosis_train <- diagnosis[train]
diagnosis_train2 <- gsub("CD", "IBD", diagnosis_train)
diagnosis_train2 <- gsub("UC", "IBD", diagnosis_train2)

diagnosis_test <- diagnosis[-train]
diagnosis_test2 <- gsub("CD", "IBD", diagnosis_test)
diagnosis_test2 <- gsub("UC", "IBD", diagnosis_test2)

ibd_filt_rf_2lvl <- randomForest(x = train_set, y = as.factor(diagnosis_train2), 
                                 importance = TRUE)

ibd_filt_rf_2lvl # look at the model
saveRDS(ibd_filt_rf_2lvl, "sandbox/ibd_sig_csv_10k_rf/ibd_filt_rf_2lvl.RDS")

pred_train2 <- predict(ibd_filt_rf_2lvl, train_set)
table(observed = diagnosis_train2, predicted = pred_train2)

# predicted
# predicted
# observed IBD nonIBD
# IBD    390      0
# nonIBD   1    138

pred_test2 <- predict(ibd_filt_rf_2lvl, test_set)
table(observed = diagnosis_test2, predicted = pred_test2) 

# predicted
# observed IBD nonIBD
# IBD    174      0
# nonIBD   8     45

varImpPlot(ibd_filt_rf_2lvl)

