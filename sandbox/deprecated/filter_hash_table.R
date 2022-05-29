library(feather)
library(edgeR)


# read in and format data -------------------------------------------------

hashes <- read_feather("all_unnormalized_abund_hashes_wide.feather")
samples <- hashes$sample                          # same sample information
samples <- gsub("_filt_named\\.sig", "", samples) # remove file suffix
#head(samples)
rownames(hashes) <- samples                       # set rownames to sample names
hashes <- hashes[ , -ncol(hashes)]                # remove the samples column
hashes <- t(hashes)                               # transpose hashes
#head(rownames(hashes))                            # check rownames are still hash
#head(colnames(hashes))                           # check colnames are still samples
colnames(hashes) <- samples                       # set colnames of hashes to samples


# filter by expression ----------------------------------------------------

#tmp <- filterByExpr(hashes, min.count = 5, min.total.count = 15)
tmp <- hashes[rowSums(hashes) > 15, ]
write_feather(as.data.frame(tmp), "all_unnomarlized_hashes_filtered.feather")
dim(tmp)
# [1] 2768882    2292
tmp2 <- filterByExpr(tmp, min.count = 5, min.total.count = 15)
table(tmp2)
# tmp2
# FALSE    TRUE
# 2768348     534
hashes_filtbyexpr <- tmp[tmp2, ]
write.csv(hashes_filtbyexpr, "all_unnormalized_abund_hashes_filtbyexpr.csv", quote = F)


# don't filter bc corncob -------------------------------------------------

hashes$hashes <- rownames(hashes)
write_feather("")
