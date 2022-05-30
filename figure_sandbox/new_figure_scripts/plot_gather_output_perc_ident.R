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
totals$db <- factor(totals$db, levels = c("refseq", "genbank", "all"))
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

pdf(file = "new_figure_scripts/upset_var_imp_genomes.pdf", width =13, height = 5)
upset(fromList(strain_list), nsets = 6, nintersects = 100, order.by = "degree",
      number.angles = 30, text.scale = c(2, 1.8, 1.8, 1.8, 1.5, 1.2))
dev.off()

pdf(file = "new_figure_scripts/upset_var_imp_genomes_small.pdf", width = 8, height = 5)
upset(fromList(strain_list), nsets = 6, nintersects = 6, order.by = "degree",
      number.angles = 30, text.scale = c(2, 1.8, 1.8, 1.35, 1.5, 1.6))
dev.off()
intersect(intersect(intersect(intersect(intersect(iHMP, PRJEB2054), PRJNA237362), PRJNA385949), PRJNA400072), SRP057027)


# count intersection genomes  ---------------------------------------------

intersect_genomes_all <- Reduce(intersect, list(PRJEB2054, PRJNA237362,
                                               iHMP, SRP057027, 
                                               PRJNA400072, PRJNA385949))

intersect_genomes_no_PRJEB2054 <- Reduce(intersect, list(PRJNA237362,
                                                        iHMP, SRP057027, 
                                                        PRJNA400072, PRJNA385949))

intersect_genomes_no_PRJNA237362 <- Reduce(intersect, list(PRJEB2054,
                                                          iHMP, SRP057027, 
                                                          PRJNA400072, PRJNA385949))

intersect_genomes_no_iHMP <- Reduce(intersect, list(PRJEB2054, PRJNA237362,
                                                   SRP057027, 
                                                   PRJNA400072, PRJNA385949))

intersect_genomes_no_SRP057027 <- Reduce(intersect, list(PRJEB2054, PRJNA237362,
                                                        iHMP,
                                                        PRJNA400072, PRJNA385949))

intersect_genomes_no_PRJNA400072 <- Reduce(intersect, list(PRJEB2054, PRJNA237362,
                                                          iHMP, SRP057027, 
                                                          PRJNA385949))

intersect_genomes_no_PRJNA385949 <- Reduce(intersect, list(PRJEB2054, PRJNA237362,
                                                          iHMP, SRP057027, 
                                                          PRJNA400072))

length(unique(c(intersect_genomes_all, intersect_genomes_no_iHMP, intersect_genomes_no_PRJEB2054,
                intersect_genomes_no_PRJNA237362, intersect_genomes_no_PRJNA385949, 
                intersect_genomes_no_PRJNA400072, intersect_genomes_no_SRP057027)))

at_least_5_genomes <- unique(c(intersect_genomes_all, intersect_genomes_no_iHMP, intersect_genomes_no_PRJEB2054,
                              intersect_genomes_no_PRJNA237362, intersect_genomes_no_PRJNA385949, 
                              intersect_genomes_no_PRJNA400072, intersect_genomes_no_SRP057027))



# add pangenomes ----------------------------------------------------------

# add in % of signature that matches to JUST the pangenomes, aka nbhd queries, 
# of the 41 genomes

files <- list.files("sandbox/41pangenomes/", ".csv$", full.names = T)

nbhd_queries <- files %>%
  set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("sandbox\\/41pangenomes\\/", "", source)) %>%
  mutate(source = gsub("_vita_vars_vs_merged_sgc_sig.csv", "", source))

nbhd_queries_totals <- nbhd_queries %>%
  group_by(source) %>%
  summarize(total = sum(f_unique_to_query)) %>%
  mutate(db = "nbhd_queries")


totals2 <- rbind(refseq_totals, genbank_totals, all_totals, nbhd_queries_totals)
totals2$db <- factor(totals2$db, levels = c("refseq", "genbank", "all", "nbhd_queries"))
pdf(file = "new_figure_scripts/percent_var_imp_in_db_with_nbhd_queries.pdf", height = 4, width = 8.5)
ggplot(totals2, aes(x = reorder(source, -total), y = total * 100, fill = db, 
                   label = round(total * 100, digits = 1))) +
  geom_col(position = "dodge") +
  theme_minimal()+
  geom_text(position = position_dodge(width = .9), vjust = 3, color = "white") +
  scale_fill_manual(values = c(refseq = "#DD8D29", genbank = "#E2D200", all = "#46ACC8", nbhd_queries = "#B40F20"),
                    labels = c(refseq = "RefSeq", genbank = "GenBank", all = "All DBs", nbhd_queries = "Nbhd Queries"),
                    name = "database") +
  labs(x = "study", y = "percent identifiable")
dev.off()


# add pangenomes + all dbs ----------------------------------------------------

# add in % of signature that matches to JUST the pangenomes, aka nbhd queries, 
# of the 41 genomes

files <- list.files("sandbox/41pangenomes/gather_results", "all_dbs.csv$", full.names = T)
nbhd_queries <- files %>%
  set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("sandbox\\/41pangenomes\\/gather_results\\/", "", source)) %>%
  mutate(source = gsub("_vita_vars_vs_merged_sgc_sig_and_all_dbs.csv", "", source))

nbhd_queries_totals <- nbhd_queries %>%
  group_by(source) %>%
  summarize(total = sum(f_unique_to_query)) %>%
  mutate(db = "nbhd_queries")


totals2 <- rbind(refseq_totals, genbank_totals, all_totals, nbhd_queries_totals)
totals2$db <- factor(totals2$db, levels = c("refseq", "genbank", "all", "nbhd_queries"))

# subtract previous db from each total to calculate the percent added by each database
totals_added <- totals2 %>%
  group_by(source, db) %>%
  arrange(match(db, c("refseq", "genbank", "all", "nbhd_queries")), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(added = total - lag(total, default = dplyr::first(total)))
totals_added$added <- ifelse(totals_added$db == "refseq", totals_added$total, totals_added$added)

pdf(file = "new_figure_scripts/percent_ident_in_db_stacked_with_nbhd_queries.pdf", height = 4, width = 8)
ggplot(totals_added, aes(x = reorder(source, -total), y = added * 100, fill = db, 
                label = round(added * 100, digits = 1))) +
  geom_col(position = position_stack(reverse = TRUE)) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 12)) +
  geom_text(color = "white", position = position_stack(vjust = 0.5, reverse = T)) + 
  ylim(c(0, 100)) +
  #scale_fill_manual(values = c(refseq = "#B997C7", genbank = "#824D99", all = "#4E78C4", nbhd_queries = "#57A2AC"),
  scale_fill_manual(values = c(refseq = "#1B7837", genbank = "#5AAE61", all = "#ACD39E", nbhd_queries = "#762A83"),
                    labels = c(refseq = "RefSeq", genbank = "GenBank", all = "All DBs", nbhd_queries = "Nbhd Queries"),
                    name = "database") +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "study", y = "percent identifiable")
dev.off()



pdf(file = "new_figure_scripts/percent_ident_in_db_stacked.pdf", height = 4, width = 8)
ggplot(totals_added %>% filter(db != "nbhd_queries"), aes(x = reorder(source, -total), y = added * 100, fill = db, 
                label = round(added * 100, digits = 1))) +
  geom_col(position = position_stack(reverse = TRUE)) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 12)) +
  geom_text(color = "white", position = position_stack(vjust = 0.5, reverse = T)) + 
  ylim(c(0, 100)) +
  scale_fill_manual(values = c(refseq = "#1B7837", genbank = "#5AAE61", all = "#ACD39E", nbhd_queries = "#762A83"),
                    labels = c(refseq = "RefSeq", genbank = "GenBank", all = "All DBs"),
                    name = "database") +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "study", y = "percent identifiable")
dev.off()



# add a plot for number of organisms identified in each study -------------

# check that there are still 1161
all %>% 
  select(name) %>%
  distinct() %>%
  nrow()

# make a plot for number of orgs in each study
all_count_orgs <- all %>% 
  select(source, name) %>%
  distinct() %>%
  group_by(source) %>%
  tally() %>%
  mutate(database = "All DBs")

all_count_orgs$source <- factor(all_count_orgs$source,
                                levels =  c("SRP057027", "PRJNA400072", "PRJNA237362", "PRJNA385949", "iHMP", "PRJEB2054"))
pdf(file = "new_figure_scripts/num_orgs_in_db.pdf", height = 4, width = 8)
ggplot(all_count_orgs, aes(x = source, y = n, fill = database, label = n)) +
  geom_col() +
  geom_text(color = "white", position = position_stack(vjust = 0.5)) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  scale_fill_manual(values = "grey") +
  labs(x = "study", y = "organisms identified")
dev.off()






