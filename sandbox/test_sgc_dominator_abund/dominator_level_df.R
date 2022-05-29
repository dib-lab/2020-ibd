library(dplyr)
library(readr)
library(tidyr)
# assign each level 1 dominator to their higher level dominators. 

# catlas <- read_csv("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/catlas.csv",
#                    col_names = c("node_id", "cdbg_vertex_id", "level", "node_id_children"))
catlas <- read_csv("~/github/spacegraphcats/dory_k21_r1/catlas.csv",
                   col_names = c("node_id", "cdbg_vertex_id", "level", "node_id_children"))
catlas <- catlas %>%
  select(-cdbg_vertex_id) %>%
  separate_rows(node_id_children, convert = T) %>%
  mutate(level = paste0("level", level))

catlas2 <-  catlas %>%
  filter(!is.na(node_id_children)) %>%
  pivot_wider(id_cols = node_id_children, names_from = level, values_from = node_id)

catlas3 <- catlas2 %>% 
  group_by(node_id_children) %>%
  summarize(across(everything(), ~ first(na.omit(.))))
