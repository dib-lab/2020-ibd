#library(tidygraph)
library(ggraph)
library(igraph)
library(readr)
library(dplyr)
setwd("~/github/2020-ibd")

sig_ccs <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_no_cdbg.tsv") %>%
  filter(mu == "mu.diagnosisCD") %>%
  mutate(label = as.character(dom_id)) %>%
  select(label, mu, estimate, bonferroni) %>%
  mutate(direction = ifelse(estimate > 0, "increased", "decreased")) %>%
  distinct()

level1 <- read.graph("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/level_1_atlas.gml", format = "gml")
level1 <- igraph::as_data_frame(level1, what = "both")
level1$vertices <- level1$vertices %>%
  left_join(sig_ccs, by = "label") %>%
  mutate(id = id + 1)
level1 <- graph_from_data_frame(d = level1$edges, vertices = level1$vertices)

pdf("tmp.pdf", height = 20, width = 20)
ggraph(level1) +
  geom_edge_link() +
  geom_node_point(size = 1, colour = direction)
dev.off()

pdf("tmp.pdf", height = 200, width = 200)
ggraph(level1, 'graphopt') + 
  geom_edge_link() +
  geom_node_point(size = 10, colour = 'steelblue')
dev.off()


# miloR graphics ----------------------------------------------------------
pdf("tmp.pdf", height = 500, width = 500)
ggraph(simplify(level1), layout = 'graphopt') +
  geom_edge_link0(edge_colour = "grey66", edge_alpha=0.2) +
  geom_node_point(aes(fill = direction), size = 12, alpha = .5, shape=21) +
  #scale_size(range =size_range, name="Nhood size") +
  #scale_edge_width(range = c(0.2,3), name="overlap size") +
  theme_classic(base_size=14) +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
dev.off()

pdf("tmp_biofabric.pdf", height =200, width = 200)
ggraph(simplify(level1), 'fabric', sort.by = node_rank_fabric()) + 
  geom_node_range(aes(color = direction)) + 
  geom_edge_span(end_shape = 'square') + 
  coord_fixed()
dev.off()
