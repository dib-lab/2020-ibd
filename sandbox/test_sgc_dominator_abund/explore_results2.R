library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
setwd("~/github/2020-ibd")



# tmp <- read_csv(domset_abund_files[1]) %>%
#   filter(level == 1)
# hist(tmp$size)
# 
# ggplot(tmp, aes(x = size)) +
#   geom_histogram() +
#   theme_minimal() +
#   scale_y_sqrt() +
#   labs(y = "square root count")
# 
# View(head(tmp))
domset_abund_files <- Sys.glob("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10_abund/*dom_abund.csv")
dom_info <- read_csv(domset_abund_files[1]) %>%
  select(-abund) %>%
  mutate(dom_id = as.character(dom_id))
df <- do.call(cbind, lapply(domset_abund_files,
                            function(x) read_csv(x)[ , "abund"]))
df <- cbind(dom_info, df)

colnames(df) <- c("dom_id", "level", "size", domset_abund_files)
colnames(df) <- gsub("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10_abund/", "", colnames(df))
colnames(df) <- gsub("_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz.dom_abund.csv", "", colnames(df))

df <- df %>%
  filter(level == 1)

row_sums_total <- rowSums(df[ , -c(1:3)])
row_sums_total <-  as.data.frame(row_sums_total) 
ggplot(row_sums_total, aes(x = row_sums_total)) +
  geom_histogram() +
  theme_minimal() +
  scale_y_sqrt()

tmp_pa <- df
rownames(tmp_pa) <- df$dom_id
tmp_pa <- tmp_pa[ , -c(1:3)]
tmp_pa[tmp_pa>0] <-1
row_sums_pa <- rowSums(tmp_pa)
row_sums_pa <- as.data.frame(row_sums_pa)
rownames(row_sums_pa) <- rownames(tmp_pa)
# filter to those that occur in 100 samples or more -- 25,910 pieces
# set filter = 100 to get pieces used for dda
filter = 0
row_sums_pa <- row_sums_pa %>%
  rownames_to_column("dom_id") %>%
  filter(row_sums_pa >= filter) 
row_sums_pa <- left_join(row_sums_pa, dom_info, by = "dom_id")

# convert wide to long to get plots in one pdf
row_sums_pa_long <- pivot_longer(row_sums_pa, cols = -dom_id, names_to = "attribute", values_to = "value") 
ggplot(row_sums_pa_long %>%
         filter(attribute %in% c("size", "row_sums_pa")),
       aes(x = value)) +
  geom_histogram() +
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  facet_wrap(~ attribute, scales = "free", ncol = 1,
             strip.position = "bottom",
             labeller = as_labeller(c(row_sums_pa = "number of query neighborhoods", size = "number of k-mers"))) +
  scale_y_sqrt(labels = scales::comma) +
  labs(x = " ", y = "number of pieces")

tmp <- row_sums_pa_long %>%
  filter(attribute == "row_sums_pa") %>%
  group_by(attribute, value) %>%
  tally()
View(tmp)

mean(df$size)
sd(df$size)
