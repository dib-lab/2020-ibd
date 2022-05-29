library(ggraph)
library(tidygraph)
library(igraph)
library(readr)
library(dplyr)
setwd("~/github/2020-ibd")

# there are two coordinate reference systems in spacegraphcats -- the cDBG nodes and the catlas nodes. 
# cDBG nodes are labeled by the original cDBG vertex ID. 
# Each node only has one vertex ID.
# Conversely, a node may have multiple catlas IDs; 
# if a node in the cDBG is a dominator at multiple levels in the catlas, the node will have a unique id at each level of the catlas.
# Most data inputs are based on the catlas ID, but the viz graphs are based on the cDBG vertex ID. 
# Therefore, the catlas is needed to act as a map between the two labelling systems.

dom_info_level1 <- read_tsv("sandbox/test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_dom_info.tsv") %>%
  filter(level == 4)

catlas_level1 <- read_csv("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/catlas.csv",
                          col_names = c("catlas_id", "cdbg_id", "level"),
                          col_types = "ddd") %>%
  filter(level == 4) %>%
  left_join(dom_info_level1, by = c("catlas_id" = "dom_id", "level"))


sig_ccs_level1 <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_no_cdbg.tsv") %>%
  filter(mu == "mu.diagnosisCD") %>%
  mutate(catlas_id = dom_id) %>%
  filter(level == 4) %>%
  left_join(catlas_level1, by = c("catlas_id", "level")) %>%
  #mutate(catlas_id = as.character(catlas_id)) %>%
  select(catlas_id, cdbg_id, level, mu, estimate, bonferroni) %>%
  mutate(direction = ifelse(estimate > 0, "increased", "decreased")) %>%
  distinct()



level1 <- read.graph("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/level_4_atlas.gml", format = "gml")
level1 <- as_tbl_graph(level1)
# tmp <- level1 %>% 
#   as_tibble() 
# table(tmp$label %in% catlas_level1$catlas_id)
# table(tmp$label %in% catlas_level1$cdbg_id)

level1 <- igraph::as_data_frame(level1, what = "both")
level1$vertices <- level1$vertices %>%
  mutate(label = as.numeric(label)) %>%
  left_join(catlas_level1, by = c("label" = "cdbg_id")) %>%
  left_join(sig_ccs_level1, by = c("label" = "cdbg_id", "catlas_id", "level")) %>%
  mutate(id = id + 1)
level1 <- graph_from_data_frame(d = level1$edges, vertices = level1$vertices)
level1 <- as_tbl_graph(level1)

ggraph(level1) +
  geom_edge_link() +
  geom_node_point(aes(fill = direction, size = size), alpha = .5, shape=21)

ggraph(simplify(level1), layout = 'graphopt') +
  geom_edge_link0(edge_colour = "grey66", edge_alpha=0.2) +
  #geom_node_point(aes(fill = direction), size = 12, alpha = .5, shape=21) +
  geom_node_point(aes(color = direction, size = size)) +
  theme_classic(base_size=14) +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
