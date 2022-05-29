# library(ggraph)
# library(tidygraph)
library(igraph)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
setwd("~/github/2020-ibd")

# read in and join graph annotations --------------------------------------

# read in dom size info
dom_info_level1 <- read_tsv("sandbox/test_corncob_dda/rgnv_nbhd_catlas_diginorm_hardtrim_piece_size_k31_r10_dom_info.tsv") %>%
  filter(level == 1)

# read in cdbg annotations
cdbg_annotations <- read_csv("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10_multifasta/multifasta.cdbg_annot.csv") %>%
  select(cdbg_id, record_name) %>%
  separate(record_name, into = "query_name", sep = " ", remove = F)

# transfer cdbg annotations to level1 dominators
cdbg_to_pieces <- read_csv("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/cdbg_to_pieces.csv")
cdbg_to_pieces <- left_join(cdbg_to_pieces, cdbg_annotations, by = c("cdbg_node" = "cdbg_id"))
dom_annotations <- cdbg_to_pieces %>%
  select(dominator, record_name, query_name) %>%
  distinct()

# add dom diff abund
sig_ccs_level1 <- read_tsv("sandbox/test_corncob_dda/corncob_results_sig_ccs_annotated_no_cdbg.tsv") %>%
  filter(mu == "mu.diagnosisCD") %>%
  mutate(catlas_id = dom_id) %>%
  filter(level == 1) %>%
  select(catlas_id, level, mu, estimate, bonferroni) %>%
  mutate(direction = ifelse(estimate > 0, "increased", "decreased")) %>%
  distinct()

# eggnog annotations
eggnog <- read_tsv("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/GCA_900036035.1_eggnog/GCA_900036035.1.emapper.annotations",
                   skip = 3) %>%
  mutate(query_name = `#query_name`) %>%
  select(-`#query_name`)

# map from catlas id (dom id) to cdbg vertex id; join with other information
catlas_level1 <- read_csv("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/catlas.csv",
                          col_names = c("catlas_id", "cdbg_id", "level"),
                          col_types = "ddd") %>%
  filter(level == 1) %>%
  left_join(dom_info_level1, by = c("catlas_id" = "dom_id", "level")) %>%
  left_join(dom_annotations, by = c("catlas_id" = "dominator")) %>%
  left_join(sig_ccs_level1, by = c("catlas_id", "level")) %>%
  left_join(eggnog, by = "query_name")



# read in graph and join with annotations ---------------------------------

# currently there are nodes missing in the gml that are present in the catlas;
# continue writing code based on current level 1 graph, and fix later when gml
# contains all nodes from all levels in the catlas

# level	nodes_in_catlas	nodes_in_gml
# 0	  NA	    1813984
# 1	  285342	28548
# 2	  11147	  11125
# 3	  4556	  4550
# 4	  1918	  1918
# 5	  841	    841
# 6	  374	    374
# 7	  150	    150
# 8	  48	    48
# 9	  30	    30
# 10	30	    NA
# 11	1	      NA

# level1 <- read.graph("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/level_1_atlas.gml", format = "gml")
# # add metadata to graph
# level1 <- igraph::as_data_frame(level1, what = "both")
# level1$vertices <- level1$vertices %>%
#   mutate(label = as.numeric(label)) %>%
#   left_join(catlas_level1, by = c("label" = "cdbg_id")) %>%
#   mutate(id = id + 1)
# # convert to tibble graph; supports repeated node entries for nodes with 
# # multiple annotations
# level1 <- as_tbl_graph(level1)

# What is the average shorthest path length between nodes with the same annotation? ----

# Step 0: Find names of nodes for pairs of annotations

duplicate_anno_names <- catlas_level1 %>%
  select(record_name, direction) %>%
  group_by(record_name, direction) %>%
  filter(!is.na(direction)) %>%
  tally() %>%
  group_by(record_name) %>%
  tally() %>%
  filter(n == 2)

duplicate_anno_cdbg_nodes <- catlas_level1 %>%
  filter(!is.na(record_name)) %>%
  filter(record_name %in% duplicate_anno_names$record_name)

# Strategy 1: identify pairs of shortest nodes (e.g. in a dataframe).
#             loop over pairs, calculate shortest path, save information

# identify nodes with the same annotation in the level 1 graph
# level1_tbl <- level1 %>% as_tibble()
# tmp <- level1_tbl %>% 
#   filter(!is.na(direction)) %>%
#   filter(!is.na(record_name)) 
# View(tmp)
# ids <- level1 %>%
#   activate(nodes) %>%
#   pull(id)
# from <- which(ids == 4708)
# to <-  which(ids == 7559)
# 
# shortest <- level1 %>%
#   morph(to_shortest_path, from, to)
# 
# shortest <- shortest %>%
#   mutate(selected_node = TRUE) %>%
#   activate(edges) %>%
#   mutate(selected_edge = TRUE) %>%
#   unmorph()
# 
# shortest <- shortest %>%
#   activate(nodes) %>%
#   mutate(selected_node = ifelse(is.na(selected_node), 1, 2)) %>%
#   activate(edges) %>%
#   mutate(selected_edge = ifelse(is.na(selected_edge), 1, 2)) %>%
#   arrange(selected_edge)
# 
# shortest %>%
#   activate(edges) %>%
#   filter(selected_edge == 2) %>%
#   as_tibble() %>%
#   summarise(total_nodes = n())

# Strategy 2: identify all shortest paths between all nodes in a graph. 
#             post process to subset to "interesting" pairs of nodes

level1 <- read.graph("sandbox/test_sgc_dominator_abund/try_preprocess_diginorm_trim_low_abund/rgnv_nbhd_diginorm_hardtrim_piece_size_k31_r10/level_1_atlas_baby.gml", format = "gml")
# check -- which diff abund nodes aren't in the graph?
table(catlas_level1$cdbg_id %in% V(level1)$label)
sig_ccs_level1_tmp <- left_join(sig_ccs_level1, catlas_level1, by = "catlas_id") %>%
  select(catlas_id, cdbg_id) %>%
  distinct()
# all but 8 of the nodes appear; maybe stick with the small graph
table(sig_ccs_level1_tmp$cdbg_id %in% V(level1)$label)
# calculate length of all shortest paths
# produces a symmetrical dataframe where only the columns are labelled by cdbg node name
distances_level1 <- distances(level1)
colnames(distances_level1) <- V(level1)$label
distances_level1 <- as.data.frame(distances_level1)

duplicate_anno_keep <- colnames(distances_level1) %in% as.character(duplicate_anno_cdbg_nodes$cdbg_id)
duplicate_anno_distances <- distances_level1[duplicate_anno_keep, duplicate_anno_keep]

duplicate_anno_distances <- duplicate_anno_distances %>%
  mutate(to_cdbg_id = colnames(.)) %>%
  pivot_longer(cols = -to_cdbg_id, names_to = "from_cdbg_id", values_to  = "shortest_path_length")

# add in to/from record_names and direction, then filter to rows where record names are the same
duplicate_anno_distances <- duplicate_anno_distances %>%
  filter(shortest_path_length != 0) %>%
  mutate(to_cdbg_id = as.numeric(to_cdbg_id)) %>%
  mutate(from_cdbg_id = as.numeric(from_cdbg_id)) 

duplicate_anno_distances <-  duplicate_anno_distances %>%
  left_join(catlas_level1, by = c("to_cdbg_id" = "cdbg_id")) %>%
  select(to_cdbg_id, to_record_name = record_name, to_direction = direction, from_cdbg_id, shortest_path_length) %>%
  left_join(catlas_level1, by = c("from_cdbg_id" = "cdbg_id")) %>%
  select(to_cdbg_id, to_record_name, to_direction, 
         from_cdbg_id, from_record_name = record_name, from_direction = direction, 
         shortest_path_length) %>%
  filter(!is.na(to_direction)) %>%
  filter(!is.na(from_direction))

duplicate_anno_distances <- duplicate_anno_distances %>%
  filter(to_record_name == from_record_name) %>%
  filter(to_direction != from_direction)

# what would contextualize these results? 
# + longest contiguous path of annotated nodes?  e.g. there's a path through 15
#   nodes all annotated as RNA polymerase?

# what is the average shortest path between an increased abundant node and the next increased abundant node? ----

# subset distances_level1 to only increased abundance nodes
increased_anno <- catlas_level1 %>%
  filter(direction == "increased")
increased_anno_keep <- colnames(distances_level1) %in% as.character(increased_anno$cdbg_id)
increased_anno_distances <- distances_level1[increased_anno_keep, increased_anno_keep]

increased_anno_distances <- increased_anno_distances %>%
  mutate(to_cdbg_id = colnames(.)) %>%
  pivot_longer(cols = -to_cdbg_id, names_to = "from_cdbg_id", values_to  = "shortest_path_length")

increased_anno_distances <- increased_anno_distances %>%
  mutate(to_cdbg_id = as.numeric(to_cdbg_id)) %>%
  mutate(from_cdbg_id = as.numeric(from_cdbg_id)) %>%
  filter(shortest_path_length != 0) %>%
  group_by(to_cdbg_id) %>%
  filter(shortest_path_length == min(shortest_path_length)) %>%
  left_join(catlas_level1, by = c("to_cdbg_id" = "cdbg_id")) %>%
  select(to_cdbg_id, to_record_name = record_name, to_direction = direction, from_cdbg_id, shortest_path_length) %>%
  left_join(catlas_level1, by = c("from_cdbg_id" = "cdbg_id")) %>%
  select(to_cdbg_id, to_record_name, to_direction, 
         from_cdbg_id, from_record_name = record_name, from_direction = direction, 
         shortest_path_length)

mean(increased_anno_distances$shortest_path_length)
sd(increased_anno_distances$shortest_path_length)
min(increased_anno_distances$shortest_path_length)
max(increased_anno_distances$shortest_path_length)

# subset distances_level1 to decreased abundance nodes
decreased_anno <- catlas_level1 %>%
  filter(direction == "decreased")
decreased_anno_keep <- colnames(distances_level1) %in% as.character(decreased_anno$cdbg_id)
decreased_anno_distances <- distances_level1[decreased_anno_keep, decreased_anno_keep]

decreased_anno_distances <- decreased_anno_distances %>%
  mutate(to_cdbg_id = colnames(.)) %>%
  pivot_longer(cols = -to_cdbg_id, names_to = "from_cdbg_id", values_to  = "shortest_path_length")

decreased_anno_distances <- decreased_anno_distances %>%
  mutate(to_cdbg_id = as.numeric(to_cdbg_id)) %>%
  mutate(from_cdbg_id = as.numeric(from_cdbg_id)) %>%
  filter(shortest_path_length != 0) %>%
  group_by(to_cdbg_id) %>%
  filter(shortest_path_length == min(shortest_path_length)) %>%
  left_join(catlas_level1, by = c("to_cdbg_id" = "cdbg_id")) %>%
  select(to_cdbg_id, to_record_name = record_name, to_direction = direction, from_cdbg_id, shortest_path_length) %>%
  left_join(catlas_level1, by = c("from_cdbg_id" = "cdbg_id")) %>%
  select(to_cdbg_id, to_record_name, to_direction, 
         from_cdbg_id, from_record_name = record_name, from_direction = direction, 
         shortest_path_length)

mean(decreased_anno_distances$shortest_path_length)
sd(decreased_anno_distances$shortest_path_length)
min(decreased_anno_distances$shortest_path_length)
max(decreased_anno_distances$shortest_path_length)

# what is the average shortest path between a differentially abundant node and a marker gene?

# what is the longest contiguous path of increased nodes in the graph? Decreased?
