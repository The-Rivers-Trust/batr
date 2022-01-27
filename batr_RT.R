# Take points representing instream barriers (e.g. dams, weirs and culverts)
# and lines representing a river network (must have consistent topology 
# e.g. Ecrins, HydroSHEDS or some other modelled drainage network). 
# Use those points to cut the lines and group the segments between the points.
# Attribute points with info like distance to source and distance to mouth.

# Josh Jones
# 8/03/2021

# Packages
# Specify the packages of interest
packages <- c("sf", "tidygraph", "tidyverse", "igraph",
              "sfnetworks", "tictoc")

# Load, or install & load, all packages
package_check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

################
# Input required
################

tic()

# data in
# river as lines
ln_path <- "./derived_data/tamar_river_clean_v3_nolocal_primacy1_buildpolylines.gpkg"

# ln_path <- "./derived_data/OSWaterNetwork_BasinID_lt5km_v2.gpkg"

# ln_path <- "./derived_data/tamar_oprvrs_polylines.gpkg"

# barriers as points
pt_path <- "./derived_data/tamar_river_obstacles_v3.gpkg"

# data out
root_source_path <- "./derived_data/river_nodes_batr.gpkg"
out_barr_path <- "./derived_data/river_batr.gpkg"
out_river_path <- "./derived_data/barriers_batr.gpkg"

# root_source_path <- "./derived_data/openrivers_river_nodes_batr.gpkg"
# out_barr_path <- "./derived_data/openrivers_river_batr.gpkg"
# out_river_path <- "./derived_data/openrivers_barriers_batr.gpkg"

# Your CRS
epsg <- 27700

# snapping distance
snap_dist <- 1000
################

# river as lines
ln <- st_read(ln_path)  %>%
  st_cast("LINESTRING") 

# barriers as points
pt <- st_read(pt_path) %>%
  mutate(is_point = TRUE) %>%
  filter(origin == "man_made")

################
# Input required
# trimming original columns
# these are the attributes to be recovered at the end
################
pt_attr <- pt
coords <- st_coordinates(pt_attr)
pt_attr$easting <- as.data.frame(coords)$X
pt_attr$northing <- as.data.frame(coords)$Y
pt_attr$geometry <- NULL
################

# function to group segments into BFL units/fragments
group_custom <- function(mode = "out") {
  # First we get the node indices of the nodes we want to route from.
  # These are:
  # --> All points that were blended into the network.
  # --> The root node at the start of the network tree.
  # Including the root will group all edges
  # that don't have a blended point downstreams.
  origins <- which(.N()$is_point | with_graph(.G(), node_is_root()))
  # Calculate the cost matrix from the origins, to all nodes in the network.
  costs <- st_network_cost(.G(), from = origins, Inf_as_NaN = TRUE, mode = mode)
  # For each node in the network:
  # --> Define which of the origins is the first
  # to be reached when travelling downstreams.
  # Remember that in the cost matrix:
  # --> The origins (the blended points + the root node) are the rows.
  # --> The destinations (all nodes in the network) are the columns.
  # Hence, we loop over the columns and keep only the minimum cost value per column.
  # We should first remove the zeros, which are the cost values from and to the same node.
  keep_minimum_cost <- function(i) {
    i[i == 0] <- NaN
    if (any(!is.na(i))) i[i != min(i, na.rm = TRUE)] <- NaN
    i
  }
  costs <- apply(costs, 2, keep_minimum_cost)
  # For each origin we know now which nodes are in its group.
  # However, we want to know which edges are in the group.
  # The cost matrix does not provide that information.
  # Finding the paths from the origins to the nodes in its group will give us this.
  # Hence, for each origin:
  # --> We compute the paths to all nodes in its group.
  # --> We extract the edge indices that are part of these paths.
  get_edge_indices <- function(i) {
    orig <- origins[i]
    dest <- which(!is.na(costs[i, ]))
    if (length(dest) > 0) {
      paths <- st_network_paths(.G(), from = orig, to = dest, mode = mode)
      edge_idxs <- do.call(c, paths$edge_paths)
      unique(edge_idxs)
    }
  }
  groups <- lapply(seq_len(nrow(costs)), get_edge_indices)
  # In tidygraph the largest group always gets group index number 1.
  # We can achieve the same by ordering the groups by number of edges.
  groups <- groups[order(lengths(groups), decreasing = TRUE)]
  # Now we can assign a group index to each edge.
  edge_idxs <- do.call(c, groups)
  group_idxs <- rep(seq_along(groups), lengths(groups))
  # The last thing left to do is to return the group indices in the correct order.
  # That is: the order of the edges in the edge table of the network.
  group_idxs[order(edge_idxs)]
}

# Snap barriers to river
st_snap_points = function(x, y, max_dist = 1000) {
  
  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)
  
  out = do.call(c,
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) 
                    return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  return(out)
}

# convert river to sfnetwork
# ln_sfnetwork <- as_sfnetwork(ln) 

# assign new basin IDs outlets to rivers 
# i.e. basins with only 1 outlet to sea, not multiple
net <- as_sfnetwork(ln) %>%
  mutate(is_point = FALSE) 

basin_sfnetwork <- net %>%
  activate("edges") %>%
  mutate(outlet = (group_custom(mode = "in")))

# extract basin rivers
edges <- st_as_sf(basin_sfnetwork, "edges") %>%
  mutate(seglen = st_length(.))

edgelen <- edges %>%
  st_set_geometry(NULL) %>%
  group_by(outlet) %>%
  summarise(outletlen = sum(seglen)) %>%
  as.data.frame()

edges_joined <- left_join(edges, edgelen, by = "outlet") %>% 
  filter(as.numeric(outletlen) > 5000)

basin_sfnetwork <- as_sfnetwork(edges_joined)

# Snap barriers
##############
## snapping is quick in QGIS so skip this step and use QGIS
## if your dataset is large
##############

# set point to TRUE so they can be distinguished from nodes later
# only keep points landing on network
# message("Snapping barriers snapped to the river...")
# pt_snap <- pt %>%
#   st_snap_points(ln, snap_dist) %>%  
#   st_as_sf() %>%
#   st_transform(epsg) %>%
#   dplyr::mutate(is_point = TRUE) %>% 
#   st_intersection(st_geometry(ln))

# combine barriers and rivers
message("Combining barriers and rivers...")
# snapping distance
tol <- units::set_units(50, "m")

# combine barriers and rivers
message("Combining barriers and rivers...")
ln_sfnetwork <- st_network_blend(basin_sfnetwork, 
                                 pt, 
                                 tol) %>%
  morph(to_subgraph, is.na(is_point)) %>%
  mutate(is_point = if_else(is_point, is_point, FALSE)) %>% # keeping barriers
  unmorph()

# plot(ln_sfnetwork)

# NOTE! Network is directed towards the root of the tree!
# i.e. flowing downstream
# Therefore we route over inbound edges instead of outbound edges to get correct results.
ln_sfnetwork <- ln_sfnetwork %>%
  activate("edges") %>%
  mutate(group = group_custom(mode = "in")) %>%
  mutate(weight = edge_length())

# # Plot the fragmented river network
# G <- ln_sfnetwork
# 
# nodes = st_as_sf(G, "nodes")
# edges = st_as_sf(G, "edges")
# edges$group = as.factor(edges$group) # Such that sf uses categorical colors.
# 
# plot(st_geometry(edges))
# plot(edges["group"], lwd = 4, key.pos = NULL, add = TRUE)
# plot(nodes[nodes$is_point, ], pch = 8, add = TRUE)
# plot(nodes[!nodes$is_point, ], pch = 20, add = TRUE)

# distance to mouth and source is only calculated for each
# barrier but it can be done for each node if necessary

# distance to source
# get a list of all the source nodes
source_vec <- with_graph(ln_sfnetwork, node_is_source())
# get their index
source_idx <- as.list(which(source_vec %in% TRUE))

# get a list of all the nodes
origins <- with_graph(ln_sfnetwork, .N()$is_point)
# get their index
origins_idx <- as.list(which(origins %in% TRUE))

# find the distance from all the nodes to all of the sources
d2s_all <- st_network_cost(ln_sfnetwork,
                           from = origins_idx,
                           to = source_idx,
                           Inf_as_NaN = TRUE,
                           mode = "in"
)

# only keep the furthest source for each node
d2s <- apply(d2s_all, 1, function(x) max(na.omit(x[x != 0])))
d2s_df <- data.frame(unlist(origins_idx), d2s) %>%
  rename(node_id = unlist.origins_idx.,
         d2s = d2s)

# distance to mouth
# find the root (river mouth)
root_vec <- with_graph(ln_sfnetwork, node_is_root())
# get its index
root_idx <- match(TRUE, root_vec)

# find the distance from each node to the root (river mouth)
d2m <- as.vector(st_network_cost(ln_sfnetwork,
                                 to = origins_idx,
                                 from = root_idx,
                                 mode = "in"
))

d2m_df <- data.frame(unlist(origins_idx), d2m) %>%
  rename(node_id = unlist.origins_idx.)

# extract fragmented rivers
edges <- st_as_sf(ln_sfnetwork, "edges") %>%
  mutate(seglen = st_length(.))
edges$group <- as.factor(edges$group) # Such that sf uses categorical colors.

edgelen <- edges %>%
  st_set_geometry(NULL) %>%
  group_by(group) %>%
  summarise(fraglen = sum(seglen)) %>%
  as.data.frame()

edges <- left_join(edges, edgelen, by = "group")

# distance upstream and downstream to the next barrier(s) (USHAB in O'Hanley models)
ds_edges <- edges %>%
  select(-to, -weight) %>%
  rename(
    DS_seglen = seglen,
    DS_fraglen = fraglen,
    DS_group = group
  ) %>%
  select(from, DS_seglen, DS_fraglen, DS_group) %>%
  as.data.frame()
ds_edges$geom <- NULL

us_edges <- edges %>%
  rename(
    US_seglen = seglen,
    US_fraglen = fraglen,
    US_group = group
  ) %>%
  select(to, US_seglen, US_fraglen, US_group) %>%
  as.data.frame()
us_edges$geom <- NULL

# add attributes to the barriers
# distance to mouth, distance to source,
# upstream length and downstream length
nodes <- st_as_sf(ln_sfnetwork, "nodes") %>%
  tibble::rowid_to_column("node_id") %>%
  left_join(d2m_df, by = "node_id") %>%
  left_join(d2s_df, by = "node_id") %>%
  left_join(ds_edges, by = c("node_id" = "from")) %>%
  left_join(us_edges, by = c("node_id" = "to"))

# Downstream ID = DSID
# DSID is the id with the minimum d2m
ds_id <- nodes %>%
  group_by(DS_group) %>%
  slice(which.min(d2m)) %>%
  select(node_id, DS_group) %>%
  rename(DSID = node_id) %>%
  as.data.frame()
ds_id$geom <- NULL

nodes <- left_join(nodes, ds_id, by = "DS_group")

# root (river mouth) and sources
root_sources <- ln_sfnetwork %>%
  activate("nodes") %>%
  mutate(root = node_is_root(), source = node_is_source()) %>%
  st_as_sf("nodes") %>% 
  select(root, source, is_point)

# Number of downstream barriers to the mouth
(barr_idx <- which(root_sources$is_point))

(root_idx <- which(root_sources$root))

paths <- st_network_paths(as_sfnetwork(ln_sfnetwork, directed = FALSE), 
                           from = root_idx,
                           to = barr_idx)["node_paths"] %>%
  as.list()

# filter out non-barrier nodes
f <- function(i){
  i[i %in% barr_idx]
}

fancyFilter <- function(f, x){
  if ( is.list( x[[1]] ) ) #only testing the first element... bad practice
    lapply( x, fancyFilter, f = f ) #recursion FTW!!
  else
    return(lapply(x, Filter, f = f ) )
}

n_ds_barr <- fancyFilter(f, list$node_paths) %>%
  tibble(barr = .) %>%
  mutate(nbarr_ds = lengths(barr)) %>% select(nbarr_ds)


cbind(out_barr, n_ds_barr)

# outputs

out_barr <- nodes %>%
  filter(is_point) %>%
  select(-is_point) %>%
  mutate(d2s = if_else(is.infinite(d2s), 0, d2s),
         US_fraglen = as.numeric(US_fraglen),
         US_fraglen = if_else(is.na(US_fraglen), 0, US_fraglen)) %>%
  cbind(n_ds_barr)

out_barr

st_write(out_barr, out_barr_path, append = FALSE)

# fragmented river
out_river <- edges %>%
  select(-weight)

out_river

st_write(out_river, out_river_path, append = FALSE)


root_sources

st_write(root_sources, root_source_path)

toc()
# OSWN
# total length 478 km
# 99 barriers
# 48.33 sec elapsed

# OSOR
# total length 730 km
# 97 barriers
# 13.72 sec elapsed

x = out_barr 
x$geom <- NULL
x %>% 
  summarise(LNEXT = mean((US_fraglen)),
            LSOURCE = mean((d2s)),
            SLENGTH = mean(d2m)) %>%
  mutate(PrimaryPathTotal = LSOURCE + SLENGTH,
         percentLNEXT = LNEXT/PrimaryPathTotal*100)

###


