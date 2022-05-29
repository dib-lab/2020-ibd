library(tidygraph)
library(ggraph)
library(igraph)
# graph <- read_csv("~/Downloads/dory_k21_r1_domset_df.csv") %>%
#   select(from = dominator, to = set_node)
# graph <- as_tbl_graph(graph)
# graph
# ggraph(graph) +
#   geom_edge_link() + 
#   geom_node_point(size = 2, colour = 'steelblue')

graph <- read_csv("~/Downloads/cdbg_to_layer1.csv") %>%
  select(-X1)
ggraph(as_tbl_graph(graph)) +
  geom_edge_link() +
  geom_node_point(size = 1, colour = 'steelblue')

graph <- read_csv("~/Downloads/layer1_to_cdbg.csv") %>%
  select(-X1)
ggraph(as_tbl_graph(graph)) +
  geom_edge_link() +
  geom_node_point(size = 1, colour = 'steelblue')

as_tbl_graph(graph) %>%
  activate(nodes) %>%
  mutate(num_nodes = graph_order())

# -------------------------------------------------------------------------

graph <- read_csv("~/Downloads/_cdbg_to_catlas.csv") %>%
  select(-X1)
ggraph(as_tbl_graph(graph)) +
  geom_edge_link() +
  geom_node_point(size = 1, colour = 'steelblue')


# cdbg gxt ---------------------------------------------------------------------


graph <- read_delim("~/Downloads/cdbg.gxt", skip =1, col_names = c("from", "to"),
                    delim = " ")
ggraph(as_tbl_graph(graph)) +
  geom_edge_link() +
  geom_node_point(size = 1, colour = 'steelblue')


# yo's level 2 ------------------------------------------------------------

graph <- read.graph("~/Downloads/tmp_dory_level_2.gml", format = "gml")
ggraph(as_tbl_graph(graph)) +
  geom_edge_link() +
  geom_node_point(size = 1, colour = 'steelblue')
