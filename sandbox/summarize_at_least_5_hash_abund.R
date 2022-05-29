library(dplyr)
library(tidyr)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, diagnosis) %>%
  mutate(library_name = gsub("-", "\\.", library_name)) %>%
  distinct()
table(info$library_name %in% all$X1)

all <- left_join(info, all, by = c("library_name" = "X1"))

head(colnames(all))

ggplot(all, aes(x = diagnosis, y = `4121156555200726`*100)) +
  geom_violin() +
  theme_minimal()

all_long_summarized <- all %>%
  pivot_longer(cols = `1119237125513`:`9141524607182576`, names_to = "hash",
               values_to = "norm_abund") %>%
  group_by(diagnosis, hash) %>%
  mutate(norm_abund = norm_abund * 100) %>%
  summarize(mean = mean(norm_abund),
            sd = sd(norm_abund)) 
View(all_long_summarized)

all_wide_summarized <- all_long_summarized %>%
  select(-sd) %>%
  pivot_wider(names_from = diagnosis, values_from = mean)

cd_up <- all_wide_summarized %>%
  filter(CD > nonIBD)

cd_down <- all_wide_summarized %>%
  filter(CD < nonIBD) %>%
  select(-UC) %>%
  mutate(diff = nonIBD - CD)


