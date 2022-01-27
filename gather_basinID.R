# gather small

file_list <- list.files("./derived_data/basin_network/", pattern = "*lt5km.gpkg", full.names = TRUE)

shapefile_list <- lapply(file_list, read_sf)

lt5km <- sf::st_as_sf(data.table::rbindlist(shapefile_list))

# gather large

file_list <- grep(list.files("./derived_data/basin_network/", full.names = TRUE), pattern = "*lt5km.gpkg",invert=TRUE, value=TRUE)

shapefile_list <- lapply(file_list, read_sf)

shapefile_list_not_empty <- shapefile_list[sapply(shapefile_list, function(x) dim(x)[1]) > 0]

gt5km <- sf::st_as_sf(data.table::rbindlist(shapefile_list_not_empty))

# write

out_path_big <- paste0("./derived_data/OSWaterNetwork_BasinID_gt5km_v2.gpkg")

out_path_small <- "./derived_data/OSWaterNetwork_BasinID_lt5km_v2.gpkg"

st_write(gt5km, out_path_big, append = FALSE)

st_write(lt5km, out_path_small, append = FALSE)

length(unique(gt5km$BasinID))
length(unique(lt5km$BasinID))
