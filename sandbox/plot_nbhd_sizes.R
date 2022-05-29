setwd("~/github/ibd")
library(readr)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(broom)
library(car)

## format file paths
files_libs <- list.files("outputs/sgc_genome_queries", pattern = ".csv$", 
                    full.names = T, recursive = T)
files_hmp <- list.files("outputs/sgc_genome_queries_hmp", pattern = ".csv$", 
                                 full.names = T, recursive = T)
files <- c(files_libs, files_hmp)

## parse sample names from file paths
names <- gsub("outputs\\/sgc_genome_queries(_hmp)?\\/", "", files)
names <- gsub("_k31_r1_search_oh0/results.csv", "", names)

## read in files
res_all <- list()
for(i in 1:length(files)){
  print(i)
  res <- read_csv(files[i])
  sample <- names[i]
  res$sample <- sample
  res_all[[i]] <- res
}

res_all <- do.call(rbind, res_all)
res_all$query <- gsub("\\/home\\/tereiter\\/github\\/ibd\\/outputs/gather_genomes\\/", "", res_all$query)
res_all$query <- gsub(".fna", "", res_all$query)


res_filt <- res_all %>%
  filter(best_containment > 0.95)


similarity = 1 / res_filt$similarity
group = 'all'
fraction = 100 / length(similarity)
cnt = 1
df2 <- data.frame(group, similarity, fraction, cnt)
head(df2)
df2$similarity <- ifelse(df2$similarity == Inf, 0, df2$similarity)
df2$bin <- cut(df2$similarity, breaks = c(0, 5.36, 18.41, 160.59, 75.118, 9859), 
               labels = c(5, 18, 160, 75, 9859))
gb <- df2 %>%
  group_by(group, similarity) %>%
  tally() 
head(gb)


# new code ----------------------------------------------------------------

lib_info <- read_tsv("inputs/working_metadata.tsv") %>% 
  select(study_accession, library_name, diagnosis) %>%
  #mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  distinct()

hmp_info <- read_tsv("inputs/hmp2_mgx_metadata.tsv") %>%
  mutate(study_accession = "hmp") %>%
  select(study_accession, External.ID, diagnosis)

colnames(hmp_info) <- colnames(lib_info)
info <- rbind(lib_info, hmp_info)

res_all <- left_join(res_all, info, by = c("sample" = "library_name"))

res_all$bp_div_kmer_count = res_all$bp / res_all$num_query_kmers
res_filt <- filter(res_all, containment > .90)
ggplot(res_filt, aes(x = diagnosis, y = bp_div_kmer_count)) +
  geom_violin()

a <- aov(formula = bp_div_kmer_count ~ diagnosis, data = res_filt)
summary(a)
posthoc <- TukeyHSD(x=a, 'diagnosis', conf.level=0.95)
posthoc


# with titus' sim code ----------------------------------------------------

res_con <- filter(res_all, best_containment > .90)
res_con$similarity <- 1/res_con$similarity
ggplot(res_con, aes(x = diagnosis, y = similarity)) +
  geom_violin()

ggplot(res_con, aes(x = diagnosis, y = query, size = similarity)) +
  geom_point()

a <- aov(formula = similarity ~ diagnosis + study_accession, data = res_con)
summary(a)
posthoc <- TukeyHSD(x=a, 'diagnosis', conf.level=0.95)
posthoc


# split by microbe and disease --------------------------------------------

res_con$diagnosis <- factor(res_con$diagnosis, levels = c("nonIBD", "CD", "UC"))

res_all$similarity2 <- 1/res_all$similarity
ggplot(res_all, aes(x = diagnosis, y = query, size = similarity2)) +
  geom_point()


# single query density --------------------------------------------------

## R. gnavus
rg <- res_all %>%
  filter(query == "GCA_900036035.1_RGNV35913_genomic") %>%
  mutate(similarity2 = 1/similarity)


ggplot(rg, aes(x = diagnosis, y = bp)) +
  geom_violin()

a <- aov(formula = contigs ~ diagnosis + study_accession, data = rg)
summary(a)
posthoc <- TukeyHSD(x=a, 'diagnosis', conf.level=0.95)
posthoc

ggplot(rg, aes(x = bp, y = contigs, color = diagnosis)) +
  geom_point(aes(alpha = .1))

# P ex...

pe <- res_all %>%
  filter(query == "ERS473343_4") %>%
  mutate(similarity2 = 1/similarity)


ggplot(pe, aes(x = diagnosis, y = bp)) +
  geom_violin()

a <- aov(formula = contigs ~ diagnosis + study_accession, data = pe)
summary(a)
posthoc <- TukeyHSD(x=a, 'diagnosis', conf.level=0.95)
posthoc

ggplot(pe, aes(x = bp, y = contigs, color = diagnosis)) +
  geom_point(aes(alpha = .1))

# heatmap -----------------------------------------------------------------

#mat_row <- data.frame(diagnosis = info$diagnosis, study = info$study_accession)
mat_row <- data.frame(diagnosis = info$diagnosis)

rownames(mat_row) <- info$library_name

tmp_pal <- brewer.pal(12, "Paired")
mat_colors <- list(diagnosis = tmp_pal[c(11, 12, 2)])
names(mat_colors$diagnosis) <- unique(mat_row$diagnosis)

mat <- res_all %>%
  select(sample, query, bp) %>%
  pivot_wider(id_cols = sample, names_from = query, values_from = bp)
rownames(mat) <- mat$sample
mat <- select(mat, -sample)

pheatmap(mat                  = log(mat), 
         show_rownames        = F,
         show_colnames        = F, 
         cluster_rows         = T, 
         cluster_cols         = T,
         color                = viridis(18),
         annotation_row       = mat_row,
         annotation_colors    = mat_colors,
         fontsize             = 14,
         annotation_names_row = F)




# scatter plot ------------------------------------------------------------
# read in tax info
gtdb <- read_tsv("sandbox/gather_vita_hashes/gtdbtk.bac120.summary.tsv") %>%
  mutate(accession = user_genome) %>%
  select(accession, classification) %>%
  mutate(classification = gsub("[dkpcofgs]__", "", classification)) %>%
  separate(data = ., col = classification, sep = ";",
           into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
  mutate(strain = NA) %>%
  mutate(source = "gtdbtk") 

res_all <- left_join(res_all, gtdb, by = c("query" = "accession"))

pdf("tmp.pdf", height = 100, width = 10)
ggplot(res_all, aes(x = query, y = log(bp), fill = diagnosis)) +
  geom_violin() +
  facet_wrap(~order, ncol = 1, drop = T, scales = "free_x") +
  theme_minimal()
dev.off()



# do anovas between num BP for each query nbhd ----------------------------


anova_wrapper <- function(data, model_expression_as_string, grouping_variable,...) {
  f_wrap <- paste0('function(.) {',model_expression_as_string,'}') %>%    parse(text=.) %>% eval
  data %>% group_by_(grouping_variable) %>% 
    do(f_wrap(.) %>% Anova(...=...) %>% tidy) %>% return
}

data=data.frame()
for (i in 1:10) {data=rbind(data,cbind(replicate=i,iris))}

aov_model_expression_as_string = 'aov(bp ~ diagnosis, data = .)'
#lm_model_expression_as_string = 'lm(Sepal.Length ~ Sepal.Width + Petal.Length , data = .)'
grouping_variable = 'query'


test <- res_all %>% 
  anova_wrapper(model_expression_as_string = aov_model_expression_as_string,
                grouping_variable = grouping_variable,type="III")
head(test)


# anova in a different way... ---------------------------------------------

res_aov <- res_all %>%
  group_by(query) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = res_all %>% map(.x = ., .f = ~ summary(aov(bp ~ diagnosis, data = .))))

res_aov <- lapply(split(res_all, res_all$query), 
                  function(d){summary(aov(bp~diagnosis, data=d))})
res_tukey <- lapply(split(res_all, res_all$query), 
                    function(d){summary(aov(bp~diagnosis, data=d))})
res_tukey <-
  res_all %>% 
  group_by(query) %>% 
  do(multitst = TukeyHSD(aov(bp ~ diagnosis, data = .))) %>%
  tidy(multitst)

res_tukey$bonferroni <- p.adjust(res_tukey$adj.p.value, method = "bonferroni", n = length(res_tukey$adj.p.value))

res_tukey_filt <- filter(res_tukey, bonferroni < .05)
