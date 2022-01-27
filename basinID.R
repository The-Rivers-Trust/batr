# creating proper basin ID on OS Water Network
# exclude basins with total length less than 5 km
# this is the input for batrRT.R

# Josh Jones 
# 01/12/2021

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

tic()

# Calculate the number of cores
no_cores <- detectCores() - 1

################
# Input required
################

# data in
# river as lines
ln_path <- "./derived_data/OSWaterNetwork_primacy1_notlocal.gpkg"

# data out
out_river_path <- "./derived_data/OSWaterNetworkBasinID.gpkg"

# Your CRS
epsg <- 27700

################

# river as lines
ln <- st_read(ln_path)  %>%
  st_cast("LINESTRING") 

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

assignBasin <- function(i) {
  message("Extracting catchments from ln with ID ", i)
  ln <- ln %>% 
    filter(catchmentid == i)
  
  # convert river to sfnetwork
  # assign new basin IDs outlets to rivers 
  # i.e. basins with only 1 outlet to sea, not multiple
  message("Creating sfnetwork...")
  net <- as_sfnetwork(ln) %>%
    mutate(is_point = FALSE) 
  
  message("Gathering edges into basins...")
  basin_sfnetwork <- net %>%
    activate("edges") %>%
    mutate(BasinID = paste0(catchmentid, "_", group_custom(mode = "in")))
  
  message("unique BasinIDs created: ", unique(basin_sfnetwork$BasinID))
  
  # extract basin rivers
  edges <- st_as_sf(basin_sfnetwork, "edges") %>%
    mutate(seglen = st_length(.))
  
  message("Do some stats...")
  edgelen <- edges %>%
    st_set_geometry(NULL) %>%
    group_by(BasinID) %>%
    summarise(basinlen = sum(seglen)) %>%
    as.data.frame()
  
  smallbasins <- edgelen  %>% 
    filter(as.numeric(basinlen) < 5000) %>% 
    distinct(BasinID) %>% 
    nrow()
  
  bigbasins <- edgelen  %>% 
    filter(as.numeric(basinlen) > 5000) %>% 
    distinct(BasinID) %>% 
    nrow()
  
  bigbasinslength <- edgelen  %>% 
    filter(as.numeric(basinlen) > 5000) %>% 
    distinct(BasinID, .keep_all = TRUE) %>% 
    summarise(sumbasinlen = sum(basinlen))
  
  biglen <- round(bigbasinslength$sumbasinlen/1000, 1)
  
  message("Excluding ", smallbasins, " basins with <5 km of total length.")
  message("Leaving ", bigbasins, " basins with >5 km of total length.")
  message("With a total of ", biglen, " km of stream length.")
  message("")
  
  edges_big <- left_join(edges, edgelen, by = "BasinID") %>% 
    filter(as.numeric(basinlen) > 5000)
  
  # keep small seperate
  edges_small <- left_join(edges, edgelen, by = "BasinID") %>% 
    filter(as.numeric(basinlen) < 5000)
  
  out_path_big <- paste0("./derived_data/basin_network/OSWaterNetwork_BasinID_", i, ".gpkg")
  out_path_small <- paste0("./derived_data/basin_network/OSWaterNetwork_BasinID_", i, "lt5km.gpkg")

  message("Writing networks with basins to OSWaterNetwork_BasinID_", i, ".gpkg")
  if (bigbasins != 0) {
    st_write(edges_big, out_path_big, append = FALSE)
  }
  
  if (smallbasins != 0) {
    st_write(edges_small, out_path_small, append = FALSE)
  }
  
}

basin_list <- na.omit(unique(ln$catchmentid))

anti_list <- c(66, 67, 68, 69, 7005, 7006, 7007, 7009, 7010, 7011, 7012, 7013, 7014, 7015, 7016, 7017, 7018, 7019, 7020, 7021, 7022, 7023, 7024, 7025, 7026, 7027, 7028, 7029, 7030, 7031, 7032, 7033, 7034, 7035, 7036, 7037, 7038, 7039, 7040, 7041, 7042, 7043, 7044, 7045, 7046, 7047, 7048, 7049, 7050, 7051, 7052, 7053, 7054, 7055, 7056, 7057, 7058, 7059, 7060, 7061, 7062, 7063, 7064, 7065, 7066, 7067, 7068, 7069, 7070, 7071, 7072, 7073, 7074, 7075, 7076, 7077, 7078, 7079, 7080, 7081, 7082, 7083, 7084, 7085, 7086, 7087, 7088, 7089, 7090, 7091, 7092, 7093, 7094, 7095, 7096, 7097, 7098, 7099, 70, 7100, 7101, 7102, 7103, 7104, 7105, 7106, 7107, 7108, 7109, 7110, 7111, 7112, 7113, 7114, 7115, 7116, 7117, 7118, 7119, 7120, 7121, 7122, 7123, 7124, 7126, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7134, 7135, 7136, 7137, 7138, 7139, 7140, 7141, 7142, 7143, 7144, 7145, 7146, 7147, 7148, 7149, 7150, 7151, 7152, 7153, 7154, 7155, 7156, 7157, 7158, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 7168, 7169, 7170, 7171, 7172, 71, 72, 73, 74)

basin_list <- as_tibble(basin_list) %>%
  filter(!value %in% anti_list) %>%
  filter(!is.na(value)) %>%
  as_vector()

basin_list <- c(7051, 7114, 7167, 7152, 7099, 7010)

tic()

# assignBasin(7152)

# out <- lapply(basin_list, assignBasin)

# initiate cluster and give each node the things
cl <- makeCluster(6, outfile = "bugfix_BasinID_071221.txt") # no_cores
clusterExport(cl, c("group_custom", "ln"))
clusterEvalQ(cl, c(library(sf),
                   library(tidygraph),
                   library(igraph),
                   library(dplyr),
                   library(readr),
                   library(sfnetworks)))

# Run
out <- parLapply(cl,
          basin_list,
          assignBasin)

# Finish
stopCluster(cl)

# bigbasinsonly <- st_as_sf(do.call(rbind, out))

# out <- lapply(basin_list, assignBasin)
# 
# bigbasinsonly <- st_as_sf(do.call(rbind, out))

# st_write(bigbasinsonly, out_river_path, append = FALSE)

toc()
