library(magrittr)
library(dplyr)
library(readr)
library(tximport)
library(corncob)
library(edgeR)

# run corncob on pangenome read count files to determine which genes are
# differentially abundant in IBD subtypes. Note that this script back-adds
# columns for samples that had no reads map (and therefore have empty files)
# during salmon quantification. 

# produce quant file list -------------------------------------------------

files <- unlist(snakemake@input[['quant']]) # read in list of quant files
# files <- list.files("sandbox/try_diffex", "quant.sf$", recursive = T, full.names = T)
files <- files[file.size(files) > 0]        # remove empty files

gather_genome <- snakemake@params[['gather_genome']]
files_root <- gsub("outputs\\/nbhd_reads_salmon\\/", "", files) # derive sample name from file name
gsub_string <- paste0(gather_genome, "_quant\\/quant\\.sf")
files_root <- gsub(gsub_string, "", files_root) # derive sample name from file name

# read in metadata --------------------------------------------------------

# read in sample metadata
info <- read_tsv(snakemake@input[["info"]]) %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct() 

# read in abundtrim lib size info
libsizes <- read_tsv(snakame@input[['mqc_fastp']], skip = 1,
                  col_names = c("library_name", "passed_filter", "low_quality", 
                                "too_many_N", "too_short")) %>%
  mutate(library_name = gsub("\\.abundtrim", "", library_name)) %>%
  group_by(library_name) %>%
  mutate(libsize = sum(passed_filter, low_quality, too_many_N, too_short)) %>%
  select(library_name, libsize)

# join metadata
info <- left_join(info, libsizes, by = "library_name")

# import counts -----------------------------------------------------------

# generate a column in the metadata info with the path to the quant.sf salmon file
info_salmon <- info %>%
  filter(library_name %in% files_root)
info_salmon <- info_salmon[order(match(info_salmon$library_name, files_root)), ] # match order of info to order of files_root
stopifnot(all.equal(info_salmon$library_name, files_root)) # check that order matches
info_salmon$salmon <- files # add col

counts <- tximport(files = info_salmon$salmon, type = "salmon", txOut = T)
colnames(counts$counts) <- info_salmon$library_name
count_info <- counts$counts

# filter to remove amino acid seqs that are 0 across an entire group
info_salmon <- info_salmon[order(match(info_salmon$library_name, colnames(count_info))), ]
# check that order matches
stopifnot(all.equal(info_salmon$library_name, colnames(count_info)))
keep <- filterByExpr(y = as.matrix(count_info), 
                     group = paste0(info_salmon$study_accession, "_", info_salmon$diagnosis),
                     min.count = 1, min.total.count = 2)
count_info <- count_info[keep, ]
write_tsv(dim(count_info), path = snakemake@output[["dim"]])

# make counts integers
count_info <- apply(count_info, 1:2, round) 

# re-add columns that had no mapping reads from salmon step
no_salmon <- info$library_name[!info$library_name  %in% info_salmon$library_name]
no_salmon_mat  <- matrix(ncol = length(no_salmon), nrow = nrow(count_info), data = 0)
colnames(no_salmon_mat) <- no_salmon
count_info <- cbind(count_info, no_salmon_mat)
write_tsv(count_info, path = snakemake@output[["counts"]])

# format counts for bdml ----------------------------------------------------

# transpose counts
count_info_t <- count_info %>% 
  as.data.frame %>% 
  t %>% 
  as_tibble

# re-add sample names
count_info_t %<>%
  add_column("sample" = colnames(count_info), .before=TRUE)

## Join up the data -- observations are rows and variables are columns
# join takes too long, use a different method
info <- info[order(match(info$library_name, count_info_t$sample)), ]
# check that order matches; fail if not
stopifnot(all.equal(info$library_name, count_info_t$sample))
df <- as.data.frame(cbind(info, count_info_t[ , -1])) 

# Run corncob -------------------------------------------------------------

## function to test for differences in all AA seqs
## note inherits df from global env
fit_corncob <- function(col_num) {
  # record protein being queried
  aa_name <- df %>%
    select(col_num) %>%
    colnames()
  # capture the warning message output by corncob when a group contains zero counts
  tc <- textConnection("messages","w")
  sink(tc, type="message")
  # run corncob with LRT
  corncob_out <- df %>%
    select(study_accession, libsize, diagnosis, all_of(col_num)) %>%
    rename(ww = 4) %>%
    corncob::bbdml(formula = cbind(ww, libsize - ww) ~ study_accession + diagnosis,
                   formula_null = cbind(ww, libsize - ww) ~ study_accession,
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
}

# run function on all AAs
all_ccs <- sapply(5:ncol(df), fit_corncob, simplify=F)
all_ccs <- do.call(rbind, all_ccs)
write_tsv(all_ccs, path = snakemake@output[["all_ccs"]])

# perfom pvalue adustment
sig <- all_ccs %>%
  filter(mu %in% c("mu.diagnosisnonIBD", "mu.diagnosisUC")) %>%
  group_by(mu) %>%
  mutate(bonferroni = p.adjust(p_vale, method = "bonferroni")) %>%
  filter(bonferroni < .05)
write_tsv(sig_ccs, path = snakemake@output[["sig_ccs"]])

