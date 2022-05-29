library(magrittr)
library(dplyr)
library(tibble)
library(readr)
library(corncob)

# run corncob on hash abundances for hashes in at least 5 of 6 models to
# determine differentially abundant hashes (e.g. direction of abund difference). 
# files outputs to sandbox/hash_corncob
# read in metadata --------------------------------------------------------

# read in sample metadata
# info <- read_tsv(snakemake@input[["info"]]) %>%
info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("-", ".", library_name)) %>%
  distinct() 

# read in abundtrim lib size info
# libsizes <- read_tsv(snakemake@input[['hash_lib_sizes']])
libsizes <- read_tsv("outputs/hash_tables/hash_lib_sizes.tsv")

# join metadata
info <- left_join(info, libsizes, by = c("library_name" = "sample"))

# import counts -----------------------------------------------------------

# count_info <- read_tsv(snakemake@input[['filt']])
count_info <- read_tsv("outputs/hash_tables/at_least_5_of_6_unnormalized_abund_hashes_wide.tsv") 
count_info$sample <- gsub("_filt_named.sig", "", count_info$sample)
count_info[is.na(count_info)] <- 0
# head(count_info[1:5, 1:5])

## Join up the data -- observations are rows and variables are columns
# join takes too long, use a different method
info <- info[order(match(info$library_name, count_info$sample)), ]
# check that order matches; fail if not
stopifnot(all.equal(info$library_name, count_info$sample))
df <- as.data.frame(cbind(info, count_info[ , -1])) 
# change levels of diagnosis so nonIBD is default
df$diagnosis <- factor(df$diagnosis, levels = c("nonIBD", "CD", "UC"))

# Run corncob -------------------------------------------------------------

## function to test for differences in all AA seqs
## note inherits df from global env
fit_corncob <- function(col_num) {
  # record protein being queried
  aa_name <- df %>%
    select(all_of(col_num)) %>%
    colnames()
  # capture the warning message output by corncob when a group contains zero counts
  tc <- textConnection("messages","w")
  sink(tc, type="message")
  # run corncob with LRT
  corncob_out <- df %>%
    select(study_accession, hash_lib_size, diagnosis, all_of(col_num)) %>%
    rename(ww = 4) %>%
    corncob::bbdml(formula = cbind(ww, hash_lib_size - ww) ~ study_accession + diagnosis,
                   formula_null = cbind(ww, hash_lib_size - ww) ~ study_accession,
                   phi.formula = ~ 1,
                   phi.formula_null = ~ 1,
                   data = .,
                   test = "LRT", 
                   boot = FALSE) %>%
    summary()
  # close message text connection
  sink(NULL, type="message")
  close(tc)
  # make a dataframe for the final results
  corncob_coeff <- data.frame(mu = rownames(corncob_out$coefficients), 
                              estimate = corncob_out$coefficients[ , 'Estimate'],
                              standard_error = corncob_out$coefficients[ , 'Std. Error'],
                              t_value = corncob_out$coefficients[ , "t value"],
                              p_value = corncob_out$coefficients[ , "Pr(>|t|)"],
                              aa_seq = aa_name,
                              separation_in_abund_model = ifelse(length(messages) == 0, "none", "separation"))
  return(corncob_coeff)
}

# run function on all AAs
all_ccs <- sapply(5:ncol(df), fit_corncob, simplify=F)
all_ccs <- do.call(rbind, all_ccs)
# write_tsv(all_ccs, path = snakemake@output[["all_ccs"]])
write_tsv(all_ccs, path = "sandbox/hash_corncob/all_ccs_hashes.tsv")

# perfom pvalue adustment
sig <- all_ccs %>%
  filter(mu %in% c("mu.diagnosisCD", "mu.diagnosisUC")) %>%
  group_by(mu) %>%
  mutate(bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
  filter(bonferroni < .05)

# write_tsv(sig, path = snakemake@output[["sig_ccs"]])
write_tsv(sig, path = "sandbox/hash_corncob/sig_ccs_hashes.tsv")
