setwd("~/github/ibd/")
library(ranger)
#library(permutations)
set.seed(1)

# https://uc-r.github.io/random_forests#tune

# read in data ------------------------------------------------------------

ibd_novalidation_filt <- read.csv("sandbox/rf_filt_vita_output/ibd_novalidation_filt.csv",
                                  row.names = 1)
diagnosis <- read.table("sandbox/rf_filt_vita_output/ibd_novalidation_filt_diagnosis.txt",
                        header = T)
diagnosis <- diagnosis$x
ibd_novalidation_filt$diagnosis <- diagnosis


# make test and train -----------------------------------------------------

train <- sample(nrow(ibd_novalidation_filt), 
                0.7*nrow(ibd_novalidation_filt), 
                replace = FALSE)
train_set <- ibd_novalidation_filt[train, ]
test_set <- ibd_novalidation_filt[-train, ]
diagnosis_train <- diagnosis[train]
diagnosis_test <- diagnosis[-train]

# tune rf -----------------------------------------------------------------

# hyperparameter grid search
hyper_grid2 <- expand.grid(
  mtry       = seq(sqrt(ncol(ibd_novalidation_filt))/2, 
                   sqrt(ncol(ibd_novalidation_filt))*8, 
                   by = 20),
  node_size  = c(5, 7, 10),
  sampe_size = c(.70, .80),
  OOB_RMSE   = 0
)

for(i in 1:nrow(hyper_grid2)) {
  
  # train model
  model <- ranger(
    formula         = diagnosis_train ~ ., 
    data            = train_set, 
    num.trees       = 10000,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sampe_size[i],
    seed            = 1
  )
  
  # add OOB error to grid
  hyper_grid2$OOB_RMSE[i] <- sqrt(model$prediction.error)
}

# hyper_grid2 %>% 
#   dplyr::arrange(OOB_RMSE) %>%
#   head(10)
# 
# cnt <- hyper_grid2 %>% 
#   dplyr::arrange(OOB_RMSE) %>%
#   group_by(OOB_RMSE) %>%
#   tally() %>%
#   arrange(n)
# 
# min(hyper_grid2$OOB_RMSE)
# hist(hyper_grid2$OOB_RMSE)


optimal_ranger <- ranger(
  formula         = diagnosis_train ~ ., 
  data            = train_set, 
  num.trees       = 10000,
  mtry            = 960,
  min.node.size   = 5,
  sample.fraction = .7,
  seed            = 1,
  importance      = 'impurity'
)

var_imp <- as.data.frame(optimal_ranger$variable.importance) 
var_imp$names <- rownames(var_imp)
colnames(var_imp) <- c("x", "names")
  
var_imp %>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  theme_minimal() +
  geom_col() +
  coord_flip() +
  ggtitle("Top 25 important variables")

pred_test <- predict(optimal_ranger, test_set)
pred_test_tab <- table(observed = diagnosis_test, predicted = pred_test$predictions)
write.table(pred_test_tab, 'outputs/optimal_rf/pred_test_tab.txt')
# predicted
# observed CD nonIBD UC
# CD     21      7  0
# nonIBD  1     74  0
# UC      3      3 30
# 14/139 = 10%

pred_train <- predict(optimal_ranger, train_set)
pred_train_tab <- table(observed = diagnosis_train, predicted = pred_train$predictions)
write.table(pred_train_tab, 'outputs/optimal_rf/pred_train_tab.txt')
# predicted
# observed  CD nonIBD  UC
# CD      70      0   0
# nonIBD   0    171   0
# UC       0      0  82

saveRDS(optimal_ranger, "outputs/optimal_rf/optimal_ranger.RDS")
