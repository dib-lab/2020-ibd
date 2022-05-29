library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

#  scp -P 2022 -i ~/.ssh/key_name.pem tereiter@farm.cse.ucdavis.edu:/home/tereiter/github/2020-ibd/all-query-results.csv .
query_res <- read_csv("~/Downloads/all-query-results.csv") %>%
  mutate(filename = gsub("/home/tereiter/github/2020-ibd/outputs/gather_matches_loso_prokka/", "", filename)) %>%
  mutate(filename = gsub(".ffn", "", filename)) %>%
  mutate(sample = gsub("_k31_r1", "", catlas_base)) %>%
  separate(col = record_name, into = c("orf_name", "annotation"), sep = " ", remove = F, extra = "merge")

# join with GTDBTK annotation information ------------------
lca <- read_csv("inputs/at_least_5_studies_vita_vars_gather_all_lca.csv")
query_res <- left_join(query_res, lca, by = c("filename" = "name_no_gz"))

# join with hash importance information --------------------
read_varimp <- function(path_optimal_rf){
  study <- gsub("outputs\\/optimal_rf\\/", "", path_optimal_rf)
  study <- gsub("_optimal_rf\\.RDS", "", study)
  optimal_rf <- readRDS(path_optimal_rf)
  varimp <- data.frame(hash = names(optimal_rf$variable.importance), 
                       importance = optimal_rf$variable.importance,
                       study = study)
  rownames(varimp) <- seq(1:nrow(varimp))
  # add a column where varimp is normalized by the total var imp
  # e.g., divide by the sum of all variable importances
  # this will make all variable importance measures sum to 1
  varimp$norm <- varimp$importance / sum(varimp$importance)
  return(varimp)
}
ihmp_varimp <- read_varimp("outputs/optimal_rf/iHMP_optimal_rf.RDS")
prjeb2054_varimp <- read_varimp("outputs/optimal_rf/PRJEB2054_optimal_rf.RDS")
prjna237362_varimp <- read_varimp("outputs/optimal_rf/PRJNA237362_optimal_rf.RDS")
prjna385949_varimp <- read_varimp("outputs/optimal_rf/PRJNA385949_optimal_rf.RDS")
prjna400072_varimp <- read_varimp("outputs/optimal_rf/PRJNA400072_optimal_rf.RDS")
srp057027_varimp <- read_varimp("outputs/optimal_rf/SRP057027_optimal_rf.RDS")
varimp <- rbind(ihmp_varimp, prjeb2054_varimp, prjna237362_varimp, 
                prjna385949_varimp, prjna400072_varimp, srp057027_varimp)
varimp_cum <- varimp %>%
  group_by(hash) %>%
  summarise(total_imp = sum(norm)) %>%
  arrange(desc(total_imp))
varimp_cum$hash <- as.numeric(as.character(varimp_cum$hash))

query_res <- left_join(query_res, varimp_cum, by = c("hashval" = "hash"))



# summarize! --------------------------------------------------------------

tmp <- query_res %>%
  group_by(hashval, filename, annotation, GTDB, total_imp) %>%
  tally()
View(tmp)


tmp2 <- query_res %>%
  select(hashval, annotation) %>%
  distinct() %>%
  group_by(annotation) %>%
  tally()
View(tmp2)  
# Questions we can ask ---------

# 1. Do the 16s/30s/gyrA etc. species matches match the checkm results?
#    (preliminarily they mostly look like they do)
# 2. Will linking the checkm results (e.g. whole microbe went up vs. whole microbe went down)
#    make these results more interpretable? E.g., remove orgs that checkm told us
#    disappear in IBD, and then only look at genes that are from other orgs and what those
#    genes represent? this should maybe get at strain variation?
# 3. When we see multiple 16s or 30s genes from a single query, is one more abundant
#    in IBD while one is more abundant in nonIBD? Can we tell if they come from different
#    16s and 30s sequences?
# 4. Why does Tyrosine recombinase XerC pop up like it's having a party? But
#    actually 42 DIFFERENT hashes annotate to this gene...why? And other than
#    16s/30s, why does this happen for other genes, too? (see tmp2 above)
# 5. Do we detect any genes from the biosynthetic gene cluster for R. gnavus?
#    Answer: No we do not. Corncob tells us these are differentially abundant,
#            meaning our method is still lossy (e.g. maybe the cluster doesn't
#            have a representative hash), or the cluster isn't important for 5/6
#            studies, or it's not important at all.
# 6. Do we see enriched pathways in RF genes? (probably need to corncob # reads
#    from each gene nbhd to determine if it's more abundant in CD/UC/nonIBD).
#    I already see oxidation genes, which fits with the idea that the gut becomes
#    leaky and oxidized in IBD, and I saw some flagellin and AMR/ABR genes.
# 7. Other than 16s/30s/gyrA etc., what genes do we see that mulitple hashes encode?
# 8. If a hash is in the nbhd of a gene from one genome, is it more likely to
#    to also be in the nbhd of another specific genome? E.g., can we detect
#    potential horizontal gene transfer partners
# 9. Should we do something with diversity of nbhd of the gene? Like compare the
#    number of k-mers (normalized for nubmer of reads?) between IBD and nonIBD...?


# explore -----------------------------------------------------------------

rgn <- query_res %>%
  filter(filename == "GCF_900036035.1_RGNV35913_genomic.fna")
View(rgn)

tmp3 <- rgn %>%
  select(hashval, total_imp, record_name, orf_name, annotation) %>%
  distinct
View(tmp3)
rgn_ps <- read_csv("~/Downloads/HPBUSUK011N-Alignment-HitTable.csv", col_names = F) %>%
  filter(X3 > 75)

table(rgn_ps$X2 %in% rgn$orf_name)
View(rgn_ps)


# ribo --------------------------------------------------------------------

ribo16s <- query_res %>%
  filter(stringr::str_detect(string = annotation, pattern = "16S ribosomal"))
length(unique(ribo16s$filename))
View(ribo16s)

ribo23s <- query_res %>%
  filter(stringr::str_detect(string = annotation, pattern = "23S ribosomal"))
length(unique(ribo23s$filename))
table(unique(ribo16s$filename) %in% unique(ribo23s$filename))

ribo50s <- query_res %>%
  filter(stringr::str_detect(string = annotation, pattern = "50S ribosomal"))
length(unique(ribo50s$filename))

ribo <- rbind(ribo16s, ribo23s, ribo50s)
View(ribo)
length(unique(ribo$filename))
