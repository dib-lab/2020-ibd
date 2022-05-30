setwd("~/github/2020-ibd")
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
#library(wesanderson)
library(ggupset)


# example -----------------------------------------------------------------

tidy_movies %>%
  distinct(title, year, length, .keep_all=TRUE) %>%
  ggplot(aes(x=Genres)) +
  geom_bar() +
  scale_x_upset(n_intersections = 20)

# try with my data --------------------------------------------------------

# first need to create a dataframe that tells us which genomes are present in
# which studies.

files_all <- list.files("outputs/gather", "all.csv$", full.names = T)

all <- files_all %>%
  set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/gather\\/", "", source)) %>%
  mutate(source = gsub("_vita_vars_all\\.csv", "", source)) %>%
  select(source, name)

# try basic upset with just this information
all %>%
  dplyr::distinct(name, .keep_all=TRUE) %>%
  ggplot(aes(x=source)) +
  geom_bar() +
  scale_x_upset()

# Then need to add information to dataframe that says how much variable 
# importance is attributable to that genome for that study.

all_totals <- all %>%
  group_by(source) %>%
  summarize(total = sum(f_unique_to_query)) %>%
  mutate(db = "all")

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