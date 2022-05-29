
iris$genus <- c(rep("A", nrow(iris)/2), rep("B", nrow(iris)/2))
# to avoid using cbind() and typing out all column names in manova() call, 
# pass a matrix where variables (abundances) are specified by subsetting brackets
# as.matrix(iris[, 1:4])

res.man <- manova(cbind(Sepal.Length, Petal.Length, Sepal.Width, Petal.Width) ~ Species*genus, data = iris)
summary(res.man)
res.man <- manova(as.matrix(iris[, 1:4]) ~ Species*genus, data = iris)
summary(res.man)
res.man

# manova on ibd -----------------------------------------------------------
setwd("~/github/ibd")
library(feather)
library(dplyr)
library(readr)

# read in abundances
ibd <- read_feather("wide_ibd2.feather")
ibd <- read_csv("sandbox/rf_filt/")
# read in info
info <- read_tsv("../inputs/working_metadata.tsv") %>% 
  select(study_accession, library_name, diagnosis) %>%
  #mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% ibd$sample) %>%
  distinct()


# sort info by ibd$sample
info <- info[match(ibd$sample, info$library_name), ]
info <- select(info, -library_name)
# remove sample vector
rownames(ibd) <- ibd$sample
ibd <- ibd[ , -ncol(ibd)]
# join info to ibd
ibd <- cbind(ibd, info)

# run manova
res_man <- manova(as.matrix(ibd[ , 1:(ncol(ibd)-2)]) ~ study_accession*diagnosis, 
                  data = ibd)
