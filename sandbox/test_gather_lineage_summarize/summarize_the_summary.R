library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(tidyr)
library(ComplexUpset)
library(UpSetR)
setwd("~/github/2020-ibd")

lineages <- list.files("sandbox/test_gather_lineage_summarize", full.names = T, pattern = ".tsv$") %>%
  set_names() %>%
  map_dfr(read_delim, delim = ";", skip=3, col_names = c("tmp", "phylum", "class", "order", "family", "genus", "species"), .id = "source")  %>%
  separate(col = tmp, into = c("level", "fraction", "superkingdom"), sep = " ") %>%
  mutate(source = gsub("sandbox/test_gather_lineage_summarize/", "", source)) %>%
  mutate(source = gsub("_lineage_summarized.tsv", "", source)) %>%
  mutate(source = gsub("_vita_vars_genbank", "", source )) %>%
  separate(col = source, into = c("study", "seed"), sep = "_") %>%
  mutate(lineage = paste(superkingdom, phylum, class, order, family, genus, species, sep = ";"))

species <- lineages %>%
  filter(level == "species") 

shared_lineages <- species %>%
  group_by(lineage) %>%
  tally() %>%
  filter(n == 36)


# upset plot to visualize intersections of lineages in all results --------

# generate a fromlist command
for(row in 1:nrow(tmp)){
  cat("'", tmp$study[row], "_", tmp$seed[row], "'= ", 
      "lineages %>% filter(study == '", tmp$study[row], "') %>% filter(seed == '", tmp$seed[row], "') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),\n",
      sep = "")
}
upset_input <- fromList(list('iHMP_seed1'= lineages %>% filter(study == 'iHMP') %>% filter(seed == 'seed1') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'iHMP_seed2'= lineages %>% filter(study == 'iHMP') %>% filter(seed == 'seed2') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'iHMP_seed3'= lineages %>% filter(study == 'iHMP') %>% filter(seed == 'seed3') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'iHMP_seed4'= lineages %>% filter(study == 'iHMP') %>% filter(seed == 'seed4') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'iHMP_seed5'= lineages %>% filter(study == 'iHMP') %>% filter(seed == 'seed5') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'iHMP_seed6'= lineages %>% filter(study == 'iHMP') %>% filter(seed == 'seed6') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJEB2054_seed1'= lineages %>% filter(study == 'PRJEB2054') %>% filter(seed == 'seed1') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJEB2054_seed2'= lineages %>% filter(study == 'PRJEB2054') %>% filter(seed == 'seed2') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJEB2054_seed3'= lineages %>% filter(study == 'PRJEB2054') %>% filter(seed == 'seed3') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJEB2054_seed4'= lineages %>% filter(study == 'PRJEB2054') %>% filter(seed == 'seed4') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJEB2054_seed5'= lineages %>% filter(study == 'PRJEB2054') %>% filter(seed == 'seed5') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJEB2054_seed6'= lineages %>% filter(study == 'PRJEB2054') %>% filter(seed == 'seed6') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA237362_seed1'= lineages %>% filter(study == 'PRJNA237362') %>% filter(seed == 'seed1') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA237362_seed2'= lineages %>% filter(study == 'PRJNA237362') %>% filter(seed == 'seed2') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA237362_seed3'= lineages %>% filter(study == 'PRJNA237362') %>% filter(seed == 'seed3') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA237362_seed4'= lineages %>% filter(study == 'PRJNA237362') %>% filter(seed == 'seed4') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA237362_seed5'= lineages %>% filter(study == 'PRJNA237362') %>% filter(seed == 'seed5') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA237362_seed6'= lineages %>% filter(study == 'PRJNA237362') %>% filter(seed == 'seed6') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA385949_seed1'= lineages %>% filter(study == 'PRJNA385949') %>% filter(seed == 'seed1') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA385949_seed2'= lineages %>% filter(study == 'PRJNA385949') %>% filter(seed == 'seed2') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA385949_seed3'= lineages %>% filter(study == 'PRJNA385949') %>% filter(seed == 'seed3') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA385949_seed4'= lineages %>% filter(study == 'PRJNA385949') %>% filter(seed == 'seed4') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA385949_seed5'= lineages %>% filter(study == 'PRJNA385949') %>% filter(seed == 'seed5') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA385949_seed6'= lineages %>% filter(study == 'PRJNA385949') %>% filter(seed == 'seed6') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA400072_seed1'= lineages %>% filter(study == 'PRJNA400072') %>% filter(seed == 'seed1') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA400072_seed2'= lineages %>% filter(study == 'PRJNA400072') %>% filter(seed == 'seed2') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA400072_seed3'= lineages %>% filter(study == 'PRJNA400072') %>% filter(seed == 'seed3') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA400072_seed4'= lineages %>% filter(study == 'PRJNA400072') %>% filter(seed == 'seed4') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA400072_seed5'= lineages %>% filter(study == 'PRJNA400072') %>% filter(seed == 'seed5') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'PRJNA400072_seed6'= lineages %>% filter(study == 'PRJNA400072') %>% filter(seed == 'seed6') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'SRP057027_seed1'= lineages %>% filter(study == 'SRP057027') %>% filter(seed == 'seed1') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'SRP057027_seed2'= lineages %>% filter(study == 'SRP057027') %>% filter(seed == 'seed2') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'SRP057027_seed3'= lineages %>% filter(study == 'SRP057027') %>% filter(seed == 'seed3') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'SRP057027_seed4'= lineages %>% filter(study == 'SRP057027') %>% filter(seed == 'seed4') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'SRP057027_seed5'= lineages %>% filter(study == 'SRP057027') %>% filter(seed == 'seed5') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage),
                 'SRP057027_seed6'= lineages %>% filter(study == 'SRP057027') %>% filter(seed == 'seed6') %>% filter(level =='species')%>%  select(lineage) %>% pull(lineage)))
#pdf("tmp_upset.pdf", height = 20, width = 150)
UpSetR::upset(data = upset_input, nsets = 50, nintersects = 10000, order.by = "degree")
#dev.off()

#pdf("tmp_upset_small.pdf", height = 10, width = 5)
UpSetR::upset(data = upset_input, nsets = 50, nintersects = 10, order.by = "degree")
#dev.off()
