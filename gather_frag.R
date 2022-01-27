# gather barr

file_list <- list.files("./derived_data/batr_output/", pattern = "RiverObstacles*", full.names = TRUE) 

file_list <- file_list[file_list != "./derived_data/batr_output/RiverObstacles_Eng_Wal_batr.gpkg" ]
file_list <- file_list[file_list != "./derived_data/batr_output/RiverObstacles_Eng_Wal_batr_v3.gpkg" ]

shapefile_list <- lapply(file_list, read_sf)

barr <- sf::st_as_sf(data.table::rbindlist(shapefile_list, fill = TRUE))

# gather riv

file_list <- list.files("./derived_data/batr_output/", pattern = "OSWaterNetwork*", full.names = TRUE)

file_list <- file_list[file_list != "./derived_data/batr_output/OSWaterNetwork_batr.gpkg"]
file_list <- file_list[file_list != "./derived_data/batr_output/OSWaterNetwork_batr_v3.gpkg"]

shapefile_list <- lapply(file_list, read_sf)

shapefile_list_not_empty <- shapefile_list[sapply(shapefile_list, function(x) dim(x)[1]) > 0]

riv <- sf::st_as_sf(data.table::rbindlist(shapefile_list_not_empty, fill = TRUE))

# write

out_river_path <- "./derived_data/batr_output/OSWaterNetwork_batr_v3.gpkg"
out_barr_path <- "./derived_data/batr_output/RiverObstacles_Eng_Wal_batr_v3.gpkg"

st_write(riv, out_river_path, append = FALSE)

st_write(barr, out_barr_path, append = FALSE)

