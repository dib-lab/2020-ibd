### Read in the data
setwd("~/github/2020-ibd")
library(tidyverse)
library(magrittr)
library(dplyr)
library(readr)
library(tximport)
library(corncob)
library(edgeR)

files <- list.files("sandbox/try_diffex", "quant.sf$", recursive = T, full.names = T)
files_root <- gsub("\\/SRS294916_20\\.fna_quant\\/quant\\.sf", "", files)
files_root <- gsub("sandbox\\/try_diffex\\/", "", files_root)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct() %>%
  filter(library_name %in% files_root)

reads <- read_tsv("sandbox/try_diffex/mqc_fastp_filtered_reads_plot_1.txt",
                  col_names = c("library_name", "passed_filter", "low_quality", "too_many_N", "too_short"),
                  skip = 1) %>%
  mutate(library_name = gsub("\\.abundtrim", "", library_name)) %>%
  group_by(library_name) %>%
  mutate(mm = sum(passed_filter, low_quality, too_many_N, too_short)) %>%
  select(library_name, mm)
  

info <- left_join(info, reads, by = "library_name")

# match order of info to order of files_root
info <- info[order(match(info$library_name, files_root)), ]

# check that order matches
all.equal(info$library_name, files_root)

info$salmon <- files

counts <- tximport(files = info$salmon, type = "salmon", txOut = T)
info <- select(info, study_accession, library_name, diagnosis)
colnames(counts$counts) <- info$library_name
count_info <- counts$counts

# filter to remove amino acid seqs that are 0 across an entire group
info <- info[order(match(info$library_name, colnames(count_info))), ]
# check that order matches
all.equal(info$library_name, colnames(count_info))
tmp <- filterByExpr(y = as.matrix(count_info), 
             group = paste0(info$study_accession, "_", info$diagnosis),
             min.count = 1, min.total.count = 2)
tmp2 <- count_info[tmp, ]
count_info <- tmp2
dim(count_info)

# make counts integers
count_info <- apply(count_info, 1:2, round) 

# transpose counts
w <- count_info %>% 
  as.data.frame %>% 
  t %>% 
  as_tibble

# re-add sample names
w %<>%
  add_column("sample" = colnames(count_info), .before=TRUE)

## Join up the data -- observations are rows and variables are columns
# join takes too long, use a different method
# df <- full_join(info, w, by = c("library_name" = "sample")) 
info <- info[order(match(info$library_name, w$sample)), ]
# check that order matches
all.equal(info$library_name, w$sample)
df <- as.data.frame(cbind(info, w[ , -1])) 
df_new <- df
## create a new column, mm, that sums up over the kmer abundances
# df$mm <- df %>% 
#   select(-study_accession, -library_name, -diagnosis) %>% 
#   rowSums()
## rearrange df so mm column is at the beginning
# df_new <- df %>% 
#   select(1:3, mm, everything())


## Look at how many gene counts are not observed across groups
# tmp <- pivot_longer(df_new, cols = `1`:`692555`, names_to = "aa_seq", values_to = "count")
# tmp <- tmp %>% 
#   group_by(study_accession, diagnosis, aa_seq) %>%
#   summarise(total = sum(count))

## analyse first amino acid sequence
tmp <-df_new %>%
  select(study_accession, mm, `6391`, diagnosis) %>%
  rename(ww = `6391`) %>%
  corncob::bbdml(formula = cbind(ww, mm - ww) ~ study_accession + diagnosis,
                 formula_null = cbind(ww, mm - ww) ~ study_accession,
                 phi.formula=~1,
                 phi.formula_null = ~ 1,
                 data = .,
                 test = "LRT", boot = FALSE)
## Plot!
my_gene_counts <- count_info["228419",]
my_gene_counts <- data.frame(x228419 = my_gene_counts, sample = names(my_gene_counts))
my_gene_counts <- left_join(my_gene_counts, info, by = c("sample" = "library_name"))
ggplot(my_gene_counts, aes(x = diagnosis, y = x228419)) +
  geom_violin() +
  theme_minimal()

## do the same thing with all of the aa seqs
fit_corncob <- function(col_num) {
  # record protein being queried
  aa_name <- df_new %>%
    select(col_num) %>%
    colnames()
  # capture the warning message output by corncob when a group contains zero counts
  tc <- textConnection("messages","w")
  sink(tc, type="message")
  # run corncob with LRT
  corncob_out <- df_new %>%
    select(study_accession, mm, diagnosis, col_num) %>%
    rename(ww = 4) %>%
    corncob::bbdml(formula = cbind(ww, mm - ww) ~ study_accession + diagnosis,
                   formula_null = cbind(ww, mm - ww) ~ study_accession,
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
                              corncob_out$coefficients, 
                              aa_seq = aa_name,
                              separation_in_abund_model = ifelse(length(messages) == 0, "none", "separation"))
}
# system.time(sapply(6:50, fit_corncob)) # 4 seconds
# system.time({all_ccs <- sapply(5:ncol(df_new), fit_corncob, simplify=F)}) # 46 s
system.time({all_ccs <- sapply(5:ncol(df_new), fit_corncob, simplify=F)})

all_ccs <- do.call(rbind, all_ccs)

sig <- all_ccs %>%
  filter(mu %in% c("mu.diagnosisnonIBD", "mu.diagnosisUC")) %>%
  group_by(mu) %>%
  mutate(bonferroni = p.adjust(Pr...t.., method = "bonferroni")) %>%
  filter(bonferroni < .05)
  

# could parallelise with `parallel::parSapply
# library(parallel)
# cl <- makeCluster(getOption("cl.cores", 2))
# system.time({all_ccs <- parallel::parSapply(cl, 5:20, fit_corncob)})


                          

# TEST --------------------------------------------------------------------

col_num <- 5
aa_name <- df_new %>%
  select(col_num) %>%
  colnames()


corncob_out <- df_new %>%
  select(study_accession, mm, diagnosis, col_num) %>%
  rename(ww = 4) %>%
  corncob::bbdml(formula = cbind(ww, mm - ww) ~ study_accession + diagnosis,
                 formula_null = cbind(ww, mm - ww) ~ study_accession,
                 phi.formula = ~ 1,
                 phi.formula_null = ~ 1,
                 data = .,
                 test = "LRT", 
                 boot = FALSE) %>%
  summary()

corncob_coeff <- data.frame(mu = rownames(corncob_out$coefficients), corncob_out$coefficients)
corncob_coeff$aa_seq <- aa_name