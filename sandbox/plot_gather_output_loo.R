setwd("~/github/2020-ibd")
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(UpSetR)
library(wesanderson)

files_all <- list.files("outputs/gather", "all.csv$", full.names = T)
files_refseq <- list.files("outputs/gather", "refseq.csv", full.names = T)
files_genbank <- list.files("outputs/gather", "genbank.csv", full.names = T)

refseq <- files_refseq %>%
  set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/gather\\/", "", source)) %>%
  mutate(source = gsub("_vita_vars_refseq\\.csv", "", source))

refseq_totals <- refseq %>%
  group_by(source) %>%
  summarize(total = sum(f_unique_to_query)) %>%
  mutate(db = "refseq")


all <- files_all %>%
  set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/gather\\/", "", source)) %>%
  mutate(source = gsub("_vita_vars_all\\.csv", "", source))

all_totals <- all %>%
  group_by(source) %>%
  summarize(total = sum(f_unique_to_query)) %>%
  mutate(db = "all")

genbank <- files_genbank %>%
  set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/gather\\/", "", source)) %>%
  mutate(source = gsub("_vita_vars_genbank\\.csv", "", source))

genbank_totals <- genbank %>%
  group_by(source) %>%
  summarize(total = sum(f_unique_to_query)) %>%
  mutate(db = "genbank")

totals <- rbind(refseq_totals, genbank_totals, all_totals)

pdf(file = "new_figure_scripts/percent_var_imp_in_db.pdf", height = 4, width = 8.5)
ggplot(totals, aes(x = reorder(source, -total), y = total * 100, fill = db, 
                   label = round(total * 100, digits = 1))) +
  geom_col(position = "dodge") +
  theme_minimal()+
  geom_text(position = position_dodge(width = .9), vjust = 3, color = "white") +
  scale_fill_manual(values = c(refseq = "#DD8D29", genbank = "#E2D200", all = "#46ACC8"),
                    labels = c(refseq = "RefSeq", genbank = "GenBank", all = "All DBs"),
                    name = "database") +
  labs(x = "study", y = "percent identifiable")
dev.off()
# wesanderson::wes_palette(name = "FantasticFox1")[3]

# upset -------------------------------------------------------------------

iHMP <- all %>% filter(source == "iHMP")
iHMP <- as.character(iHMP$name)

PRJEB2054 <- all %>% filter(source == "PRJEB2054")
PRJEB2054 <- as.character(PRJEB2054$name)

PRJNA237362 <- all %>% filter(source == "PRJNA237362")
PRJNA237362 <- as.character(PRJNA237362$name)

PRJNA237362 <- all %>% filter(source == "PRJNA237362")
PRJNA237362 <- as.character(PRJNA237362$name)

PRJNA385949 <- all %>% filter(source == "PRJNA385949")
PRJNA385949 <- as.character(PRJNA385949$name)

PRJNA400072  <- all %>% filter(source == "PRJNA400072")
PRJNA400072  <- as.character(PRJNA400072$name)

SRP057027  <- all %>% filter(source == "SRP057027")
SRP057027  <- as.character(SRP057027$name)

strain_list <- list(iHMP = iHMP,
                    PRJEB2054 = PRJEB2054,
                    PRJNA237362 = PRJNA237362,
                    PRJNA385949 = PRJNA385949,
                    PRJNA400072 = PRJNA400072,
                    SRP057027 = SRP057027)
UpSetR::upset(fromList(strain_list), nsets = 6, nintersects = 100, order.by = "degree")


intersect(intersect(intersect(intersect(intersect(iHMP, PRJEB2054), PRJNA237362), PRJNA385949), PRJNA400072), SRP057027)
