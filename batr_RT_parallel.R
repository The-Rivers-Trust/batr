# Take points representing instream barriers (e.g. dams, weirs and culverts)
# and lines representing a river network (must have consistent topology 
# e.g. Ecrins, HydroSHEDS or some other modelled drainage network). 
# Use those points to cut the lines and group the segments between the points.
# Attribute points with info like distance to source and distance to mouth.

# in parallel

# Josh Jones
# 7/12/2021

# Packages
# Specify the packages of interest
packages <- c("sf", "tidygraph", "tidyverse", "igraph",
              "sfnetworks", "tictoc", "parallel")

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

# data in
# river as lines
ln_path <- "./derived_data/OSWaterNetwork_BasinID_gt5km_v2_multi2singlepart.gpkg"

# barriers as points
pt_path <- "./derived_data/RiverObstacles_Eng_Wal_50m_BasinID.gpkg"

# Your CRS
epsg <- 27700

# snapping distance
snap_dist <- 50
################

# river as lines
ln <- st_read(ln_path)  %>%
  st_cast("LINESTRING") %>%
  st_zm(drop = TRUE, what = "ZM")

# barriers as points
pt <- st_read(pt_path) %>%
  mutate(is_point = TRUE) %>%
  filter(origin == "man_made") %>%
  filter(!is.na(BasinID))

################
# Input required
# trimming original columns
# these are the attributes to be recovered at the end
################
# pt_attr <- pt
# coords <- st_coordinates(pt_attr)
# pt_attr$easting <- as.data.frame(coords)$X
# pt_attr$northing <- as.data.frame(coords)$Y
# pt_attr$geometry <- NULL
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

frag <- function(x) {
  message("Filtering barriers and rivers for BasinID ", x)
  # only work on one basin at a time
  ln <- filter(ln, BasinID == x) 
  pt <- pt %>%
    filter(BasinID == x)

  # check if there are barriers
  n_pt <- nrow(pt)
  if (n_pt < 1) {
    out_river_path <- paste0("./derived_data/batr_output/OSWaterNetwork_batr_BasinID_", x, ".gpkg")
    message("Writing outputs to ", out_river_path, " ...")
    st_write(ln, out_river_path, append = FALSE)
    return(pt) # no fragmentation, skip to next basin
  }

  # convert river to sfnetwork
  ln_sfnetwork <- as_sfnetwork(ln) 
  
  # combine barriers and rivers
  message("Combining barriers and rivers...")
  # snapping distance
  tol <- units::set_units(snap_dist, "m")

  ln_sfnetwork <- st_network_blend(ln_sfnetwork, 
                                   pt, 
                                   tol) %>%
    morph(to_subgraph, is.na(is_point)) %>%
    mutate(is_point = FALSE) %>%
    unmorph()

  # NOTE! Network is directed towards the root of the tree!
  # i.e. flowing downstream
  # Therefore we route over inbound edges instead of outbound edges to get correct results.
  message("Grouping fragments...")
  # tryCatch({
  #   ln_sfnetwork <- ln_sfnetwork %>%
  #     activate("edges") %>%
  #     mutate(group = group_custom(mode = "in")) %>%
  #     mutate(weight = edge_length())
  # },
  # error = function(e) {
  #   message("Error... Trying to smooth network...")
  #   # convert river to sfnetwork
  #   ln_sfnetwork <- as_sfnetwork(ln) %>% 
  #     convert(to_spatial_smooth)
  #   return(ln_sfnetwork)
  #   message("Success!")
  # })
  
  ln_sfnetwork <- ln_sfnetwork %>%
    activate("edges") %>%
    mutate(group = group_custom(mode = "in")) %>%
    mutate(weight = edge_length())
  
  # distance to mouth and source is only calculated for each
  # barrier but it can be done for each node if necessary
  
  message("Finding distance to source and mouth...")
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
                             mode = "in")
  
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
                                   mode = "in"))
  
  d2m_df <- data.frame(unlist(origins_idx), d2m) %>%
    rename(node_id = unlist.origins_idx.)
  
  message("Extracting and joining outputs...")
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
    mutate(root = node_is_root(), 
           source = node_is_source()) %>%
    st_as_sf("nodes") %>% 
    select(root, source, is_point)
  
  # # number of upstream barriers towards all sources
  # barr_idx <- which(root_sources$is_point)
  # 
  # source_idx <- which(root_sources$source)
  # 
  # us_paths <- st_network_paths(as_sfnetwork(ln_sfnetwork, directed = FALSE), 
  #                           from = barr_idx,
  #                           to = source_idx)["node_paths"] %>%
  #   as.list()
  # 
  # # filter out non-barrier nodes
  # f <- function(i){
  #   i[i %in% barr_idx]
  # }
  # 
  # fancyFilter <- function(f, x){
  #   if ( is.list( x[[1]] ) ) # only testing the first element... bad practice
  #     lapply( x, fancyFilter, f = f ) # recursion FTW!!
  #   else
  #     return(lapply(x, Filter, f = f ) )
  # }
  # 
  # n_us_barr <- fancyFilter(f, us_paths$node_paths) %>%
  #   tibble(barr = .) 
  # 
  # message("n_us_barr")
  # message(head(n_us_barr))
  
  # %>%
  #   mutate(nbarr_ds = lengths(barr)) %>% 
  #   select(nbarr_ds)
  
  
  # Number of downstream barriers towards the mouth (root)
  barr_idx <- which(root_sources$is_point)
  
  root_idx <- which(root_sources$root)
  
  paths <- st_network_paths(as_sfnetwork(ln_sfnetwork, directed = FALSE), 
                            from = root_idx,
                            to = barr_idx)["node_paths"] %>%
    as.list()
  
  
  n_ds_barr <- fancyFilter(f, paths$node_paths) %>%
    tibble(barr = .) %>%
    mutate(nbarr_ds = lengths(barr)) %>% 
    select(nbarr_ds)
  
  # outputs
  # barriers with original attributes back to the points
  out_barr <- nodes %>%
    filter(is_point) %>%
    select(-is_point) %>%
    mutate(d2s = if_else(is.infinite(d2s), 0, d2s),
           US_fraglen = as.numeric(US_fraglen),
           US_fraglen = if_else(is.na(US_fraglen), 0, US_fraglen))
  
  # sometimes there is a mismatch in row numbers
  if (nrow(n_ds_barr) == nrow(out_barr)) {
    out_barr <- cbind(out_barr, n_ds_barr)
  } else {
    message("Mismatch in row number for N downstream barriers. Using st_join() instead of cbind()...")
    nbarr_ds_geom <- cbind(filter(root_sources, is_point), n_ds_barr) %>%
      select(nbarr_ds, geom)
    out_barr <- st_join(out_barr, nbarr_ds_geom)
  }

  # fragmented river
  out_river <- edges %>%
    select(-weight)
  
  # data out
  out_river_path <- paste0("./derived_data/batr_output/OSWaterNetwork_batr_BasinID_", x, ".gpkg")
  out_barr_path <- paste0("./derived_data/batr_output/RiverObstacles_Eng_Wal_batr_BasinID_", x, ".gpkg")
  
  message("Writing outputs to ", out_river_path, " and ", out_barr_path, " ...")
  st_write(out_barr, out_barr_path, append = FALSE)
  st_write(out_river, out_river_path, append = FALSE)
  # st_write(root_sources, root_source_path)
}

basinIDs <- unique(ln$BasinID)

tic()
frag("7016_1")
toc()

# remove BasinIDs listed below because they're not a "true" outlet basin
# they're endorheic

# not sure why this doesn't work
# I think it's because the point is snapping to either source or root node
# 7021_4 7028_6 7162_11 7070_13

# also singular straight lines...

# when processing is all done extract the BasinIDs from this list and figure
# it out

badchannels <- c("69_10", "7156_10", "7100_6", "7167_16", 
                 "7167_19", "7167_10", 
                 "7100_3",
                 "7156_8", "7132_137", "7100_4", 
                 "7051_3", "7051_5", "7156_6", "7013_3", 
                 "7013_4", "7019_3", "7152_7", "7028_7", "7021_4", "7028_6",
                 "7162_11", "7123_7", "7152_9", "7070_13", "7040_3", "7042_9",
                 "7045_2","7095_14", "7127_5", "7133_6", "7134_4", "7140_2",
                 "7143_12", "7152_48", "7167_17") 


anti_list <- read_csv("./derived_data/anti_list_8.csv",
                      show_col_types = FALSE)  %>%
  filter(!value %in% badchannels) %>%
  as_vector()

# done <- c("7040_1", "7041_1", "7040_5", "7041_2", "7042_11",
#           "7045_5", "7045_4", "7045_3", "7092_7", "7093_3", "7093_2", "7093_1",
#           "7092_13", "7092_8", "7094_1", "7162_8", "7153_4",
#           "7153_2", "7153_1", "7095_6", "7141_5", "7095_8",
#           "7095_9", "7095_5", "7095_8", "7134_1", "7141_1",
#           "7141_3", "7141_2", "7144_2", "7144_1", "7153_3")

exclude <- append(badchannels, anti_list)
# exclude <- append(exclude, done)

# exclude <- badchannels

basinIDs <- as_tibble(basinIDs) %>%
  filter(!value %in% exclude) %>%
  filter(!is.na(value)) %>%
  as_vector()

basinIDs <- anti_list

basinIDs <- c("7013_1",
              "7104_1",
              "7140_1",
              "7152_1")

# bmoved <- c('7028_6', '7045_2', '7040_3', '7042_9')
# 
# 
# out <- lapply(bmoved,
#               frag)

tic()

# initiate cluster and give each node the things
cl <- makeCluster(2, outfile = "bugfix_batr_200122.txt") # no_cores
clusterExport(cl, c("group_custom", "ln", "pt", "snap_dist"))
clusterEvalQ(cl, c(library(sf),
                   library(tidygraph),
                   library(igraph),
                   library(dplyr),
                   library(readr),
                   library(sfnetworks)))

# Run
out <- parLapply(cl,
                 basinIDs,
                 frag)

# Finish
stopCluster(cl)

toc()


