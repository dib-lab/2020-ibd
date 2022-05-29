library(tidygraph)
library(ggraph)
library(igraph)
setwd("~/github/2020-ibd")

level7 <- read.graph("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/level_7_atlas.gml", format = "gml")
level7 <- as_tbl_graph(level7)
ggraph(level7) +
  geom_edge_link() +
  geom_node_point(size = 1, colour = 'steelblue')

pdf("tmp.pdf", height = 200, width = 200)
ggraph(level7, 'graphopt') + 
  geom_edge_link() +
  geom_node_point(size = 10, colour = 'steelblue')
dev.off()
