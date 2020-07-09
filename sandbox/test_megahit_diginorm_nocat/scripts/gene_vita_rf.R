library(ranger)
library(dplyr)
library(readr)
library(feather)
library(Pomona)

set.seed(1)

# perform variable selection within random forests on pangenome gene abundances. 
# Uses the Pomona package vita implementation, which wraps the ranger
# package. Saves output as RDS for faster loading of output objects into
# subsequent R sessions.

## read in data ------------------------------------------------------------

## format hash table (samples x features)

vsd <- read_tsv(snakemake@input[['vsd']]) # read in hash abund table
vsd <- as.data.frame(vsd)                 # transform to dataframe
rownames(vsd) <- vsd$protein              # set rownames to gene name
vsd <- vsd[ , -ncol(vsd)]                 # remove "protein" column
vsd <- t(vsd)                             # transpose so samples are rows, gense are cols

## read in study metadata
## collapse duplicate libraries so each sample only has one row
info <- read_tsv(snakemake@input[['info']]) %>%
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(vsd)) %>%
  distinct()

## remove validation cohort from variable selection
info_novalidation <- info %>%
  filter(study_accession != snakemake@params[["validation_study"]]) %>%
  mutate(library_name = gsub("-", "\\.", library_name))
vsd_novalidation <- vsd[rownames(vsd) %in% info_novalidation$library_name, ]

## make classification vector
## match order of to vsd
info_novalidation <- info_novalidation[match(rownames(vsd_novalidation), info_novalidation$library_name), ]
## make diagnosis var
diagnosis_novalidation <- info_novalidation$diagnosis

# run vita ----------------------------------------------------------------

## perform variant selection
## var.sel.vita calculates p-values based on the empirical null distribution
## from non-positive VIMs as described in Janitza et al. (2015).
vsd_vita <- var.sel.vita(x = vsd_novalidation, y = diagnosis_novalidation, p.t = 0.05,
                         ntree = 10000, mtry.prop = 0.2, nodesize.prop = 0.1,
                         no.threads = snakemake@params[["threads"]], 
                         method = "ranger", type = "classification")
saveRDS(vsd_vita, snakemake@output[["vita_rf"]])

# write files -------------------------------------------------------------

## write predictive hashes
var <- vsd_vita$var                 # separate out selected predictive hashes
var <- gsub("X", "", var)           # remove the X from the beginning of hashes
write.table(var, snakemake@output[['vita_vars']],
            quote = F, col.names = F, row.names = F)

## filter to predictive hashes and write to fie
vsd_filt <- vsd[ , colnames(vsd) %in% var] # subset vsd to hashes in vsd_vita
write.csv(vsd_filt, snakemake@output[['vsd_filt']], quote = F)
