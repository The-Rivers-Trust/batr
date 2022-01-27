# analysis batr output vs inputs

# Josh Jones 
# 21/12/2021

# packages
library("tidyverse", "sf")

# data in
# all river obs
rivobs <- read_csv("raw_data/River_Obstacles.csv")

# just england and wales
rivobs_eng_wal <- sf::st_read("./derived_data/RiverObstacles_Eng_Wal.gpkg")
rivobs_eng_wal$geom <- NULL

# only batr outputs
batrout <- sf::st_read("derived_data/batr_output/RiverObstacles_Eng_Wal_batr_v3.gpkg") %>%
  mutate(d2m = if_else(!is.finite(DS_fraglen), 0, DS_fraglen))
batrout_nogeom <- batrout
batrout_nogeom$geom <- NULL

# river for stats
oswn <- sf::st_read("./derived_data/batr_output/OSWaterNetwork_batr_v3.gpkg")
oswn_nogeom <- oswn
oswn_nogeom$geom <- NULL

oswn_catchlen <- oswn_nogeom %>% 
  group_by(catchmentid) %>%
  summarise(SYSlength = sum(length))

colnames(rivobs)
colnames(batrout)

###
# rename, join, manipulate etc
###
oswn_join <- oswn %>% 
  select(inspireid, catchmentid, geom) %>%
  left_join(oswn_catchlen, by = "catchmentid")
batrout_buffer <- sf::st_buffer(batrout %>% select(CoreoID, geom), 1)
batr_oswn <- st_join(batrout_buffer, oswn_join)
batr_oswn$geom <- NULL
batr_oswn <- batr_oswn %>% distinct(CoreoID, .keep_all = TRUE)

batrout_stats <- batrout %>%
  select(CoreoID, BasinID ,d2m, d2s, US_fraglen, nbarr_ds) %>%
  left_join(batr_oswn, by = "CoreoID") %>%
  rename(UQIDO = inspireid,
         SYSID = catchmentid,
         LNEXT = US_fraglen,
         LSOURCE = d2s,
         SLENGTH = d2m,
         STOTAL = nbarr_ds) %>%
  mutate(PrimaryPathTotal = LSOURCE + SLENGTH,
         `%toNext` = round(LNEXT/PrimaryPathTotal * 100, 2),
         `%DS` = round(SLENGTH/PrimaryPathTotal * 100, 2),
         `%PPT_SystemLength` = round(PrimaryPathTotal/SYSlength * 100, 2)) %>%
  select(CoreoID, 
         UQIDO,
         SYSID,
         SYSlength,
         LNEXT,
         LSOURCE,
         SLENGTH,
         STOTAL,
         PrimaryPathTotal,
         `%toNext`,
         `%DS`,
         `%PPT_SystemLength`, geom)

# batrout_stats$geom <- NULL

rivobs_v2 <- left_join(rivobs, batrout_stats, by = "CoreoID")

write_csv(rivobs_v2, "./derived_data/River_Obstacles_frag_stats.csv")

st_write(rivobs_v2 %>% 
           select(-FID), "./derived_data/River_Obstacles_frag_stats.gpkg")

### summaries ###
rivobs_v2 %>%
  summarise(total = n())

rivobs_v2 %>% 
  filter(origin == "man_made") %>%
  summarise(withstats = n())

rivobs_v2 %>% 
  filter(!is.na(UQIDO)) %>%
  summarise(withstats = n())

mean <- rivobs_v2 %>% 
  filter(!is.na(UQIDO)) %>%
  select(LNEXT, LSOURCE, SLENGTH, 
         STOTAL, PrimaryPathTotal, 
         `%toNext`, `%DS`, `%PPT_SystemLength`) %>% 
  summarise_all(funs(mean(., na.rm = TRUE))) %>%
  mutate(stat = "mean") %>%
  pivot_longer(!stat, names_to = "metric", values_to = "val")

median <- rivobs_v2 %>% 
  filter(!is.na(UQIDO)) %>%
  select(LNEXT, LSOURCE, SLENGTH, 
         STOTAL, PrimaryPathTotal, 
         `%toNext`, `%DS`, `%PPT_SystemLength`) %>% 
  summarise_all(funs(median(., na.rm = TRUE))) %>%
  mutate(stat = "median") %>%
  pivot_longer(!stat, names_to = "metric", values_to = "val")

sd <- rivobs_v2 %>% 
  filter(!is.na(UQIDO)) %>%
  select(LNEXT, LSOURCE, SLENGTH, 
         STOTAL, PrimaryPathTotal, 
         `%toNext`, `%DS`, `%PPT_SystemLength`) %>% 
  summarise_all(funs(sd(., na.rm = TRUE))) %>%
  mutate(stat = "sd") %>%
  pivot_longer(!stat, names_to = "metric", values_to = "val")

fragstats <- rbind(mean, median, sd) %>%
  pivot_wider(names_from = stat, values_from = val)

write_csv(fragstats, "./derived_data/fragstats_summary.csv")

rivobs_eng_wal %>% summarise(n = n())

rivobs_eng_wal %>% 
  filter(origin == "man_made") %>% 
  summarise(n = n())


##################################

# check batr column integrity
# first check NA
# batrout %>% filter(is.na(d2m)) %>% distinct(BasinID)
# batrout %>% filter(is.na(d2s)) %>% distinct(BasinID)
batrout %>% filter(is.na(DS_fraglen)) %>% 
  group_by(BasinID) %>% 
  summarise(n = n()) # NAs are barriers at the mouth
# batrout %>% filter(is.na(US_fraglen)) %>% 
#   group_by(BasinID) %>% 
#   summarise(n = n())
batrout %>% filter(is.na(nbarr_ds)) %>% 
  group_by(BasinID) %>% 
  summarise(n = n()) # NAs are on canals

# then check Inf
batrout %>% filter(!is.finite(d2m)) %>% 
  group_by(BasinID) %>% 
  summarise(n = n()) # on stranded channel
# batrout %>% filter(!is.finite(d2s)) %>% 
#   distinct(BasinID)
batrout %>% filter(!is.finite(DS_fraglen)) %>% 
  group_by(BasinID) %>% 
  summarise(n = n()) # Infs are barriers at the mouth ### fix these
# batrout %>% filter(!is.finite(US_fraglen)) %>% 
#   group_by(BasinID) %>% 
#   summarise(n = n())
batrout %>% filter(!is.finite(nbarr_ds)) %>% 
  group_by(BasinID) %>% 
  summarise(n = n()) # stranded canals