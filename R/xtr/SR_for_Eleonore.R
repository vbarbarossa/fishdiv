
#' *****************************************************************************************************************************
#' # Species richness data sampling
#' 
#' 
#' **Q: should we exclude species occurring in endoreic basins?**

#' load packages needed
valerioUtils::libinv(c('dplyr','sf'))

#' intersect the catchments with the HB12 points to assign the gsim id to each catchment

# load hydrobasins 12 shpefile points layer
hb12_p <- read_sf('~/surfdrive/Documents/projects/fishsuit/data/hybas12_points_nolakes.gpkg') %>%
  dplyr::select(HYBAS_ID,MAIN_BAS)

# load the previously produced catchments layer
catch <- read_sf('~/surfdrive/Documents/projects/Hydrographies/ws/basins_5min_pcrglobwb_adjusted.gpkg')

# check crs
st_crs(hb12_p) == st_crs(catch)

# get the sparse list of intersections
lst <- st_intersects(catch,hb12_p)

# number of catchments without hydrobasins:
length(lst) - (lst %>% sapply(.,function(x) length(x) != 0) %>% sum)

# ids of hydrobasins
hybas_id <- hb12_p$HYBAS_ID

# compile table with HYBAS_ID and corresponding ID
names(lst) <- catch$ID
catch_hb <- lst[which(unname(lst %>% sapply(.,function(x) length(x) != 0)))] %>% 
  lapply(seq_along(.), function(n,v,i){data.frame(HYBAS_ID = hybas_id[v[[i]]],ID = n[i])},n=names(.),v = .) %>%
  do.call('rbind',.) %>%
  as_tibble()

# load Tedesco data for filtering out non-native species
# sample Tedesco basins
bas <- read_sf('~/surfdrive/Documents/projects/fishsuit/data/Tedesco/Basin042017_3119.shp') %>% mutate(id = 1:nrow(.)) 
bas_ras <-  fasterize::fasterize(sf = bas,raster = raster::raster(res = 1/12), field = 'id')
hb12_p$ws <- raster::extract(bas_ras,hb12_p)
# check most recurring ted basin id per main bas, and assign that value to all the sub-basins belonging to the main bas
hb12_main_bas <- hb12_p %>%
  as_tibble() %>% select(-geom) %>%
  filter(!is.na(ws)) %>%
  group_by(MAIN_BAS) %>%
  summarise(
    bas = table(ws) %>% sort(decreasing = T) %>% .[1] %>% names
  )

# merge back with hb12_p
hb12 <- hb12_p %>% as_tibble() %>% select(-geom) %>% left_join(hb12_main_bas) %>% select(-ws)

ted_occ <- read.csv('~/surfdrive/Documents/projects/fishsuit/data/Tedesco/Occurrence_Table.csv',sep = ';') %>%
  as_tibble() %>%
  rename(BasinName = X1.Basin.Name) %>%
  left_join(bas %>% as_tibble() %>% select(-geometry,BasinName,id)) %>%
  select(id,BasinName,binomial = X6.Fishbase.Valid.Species.Name, Status = X3.Native.Exotic.Status)
# exotic native 
# 8128 102441 
# complement the table with Su et al., 2021 updated info on exotic species and extinctions
su_occ <- read.csv('data/Su_et_al_Science_table.csv') %>%
  as_tibble() %>%
  rename(BasinName = Basin.Name) %>%
  left_join(bas %>% as_tibble() %>% select(-geometry) %>% select(BasinName,id)) %>%
  select(id,BasinName,binomial = Species.Latin.Name, Status)
su_occ$Status[su_occ$Status == 'extinction/extirpation'] <- 'extinct'
su_occ$Status[su_occ$Status == 'introduction'] <- 'exotic'
# exotic extinct 
# 7219     292 
#' check fish endemism at the basin scale: if a fish species occurs in only one main basin then it is considered endemic

# all species names excluding exclusively lentic and other filtering steps as in fishsuit
# correspondence table with fishbase names (same used in fish_sp_name)
# fish_sp_names <- foreign::read.dbf('../fishsuit/proc/species_ranges_raw.dbf') %>%
#   as_tibble() %>% select(binomial,fishbase_1) %>% 
#   filter(fishbase_1 %in% (read.csv('../fishsuit/proc/thresholds_average_filtered.csv')%>%
#                             pull(binomial) %>% as.character()))

fish_sp_names <- read_sf('~/surfdrive/Documents/projects/fishsuit/proc/species_ranges_merged.gpkg') %>%
  as_tibble() %>% select(-geom) %>% 
  filter(binomial %in% (read.csv('~/surfdrive/Documents/projects/fishsuit/proc/thresholds_average_filtered.csv')%>%
                          pull(binomial) %>% as.character()))


# get the fish data and filter out exclusively lentic
fish <- readRDS('~/surfdrive/Documents/projects/fishsuit/proc/species_ranges_merged_on_hybas12.rds') %>%
  inner_join(fish_sp_names) %>%
  as_tibble() %>%
  select(HYBAS_ID, binomial) %>%
  distinct()
length(unique(fish$binomial)) #11,425

# filter out non-native and extinct based on ted_occ and su_occ
# prepare a list of exotic fish per ted basin id
exotic <- ted_occ %>%
  bind_rows(su_occ) %>%
  distinct() %>%
  filter(Status == 'exotic') %>%
  select(id,binomial) %>%
  mutate(binomial  = gsub('\\.',' ',binomial)) %>%
  arrange(id)

extinct <- su_occ %>%
  distinct() %>%
  filter(Status == 'extinct') %>%
  select(id,binomial) %>%
  mutate(binomial  = gsub('\\.',' ',binomial)) %>%
  arrange(id)

fish <- fish %>% left_join(hb12)
fish$bas[is.na(fish$bas)] <- 0

fish_split <- fish %>%
  split(fish$bas)

fish_1 <- fish_split[names(fish_split) %in% as.character(unique(c(exotic$id,extinct$id)))]
fish_2 <- fish_split[!names(fish_split) %in% names(fish_1)]

fish_filtered <- lapply(names(fish_1),function(x){
  t <- fish_1[x][[1]]
  e <- exotic[as.character(exotic$id) == x,] %>% pull(binomial)
  e <- unique(c(e,extinct[as.character(extinct$id) == x,] %>% pull(binomial)))
  t <- t[!t$binomial %in% e,]
  return(t)
}) %>% do.call('rbind',.) %>%
  bind_rows(fish_2 %>% do.call('rbind',.))

fish_exotic <- lapply(names(fish_1),function(x){
  t <- fish_1[x][[1]]
  e <- exotic[as.character(exotic$id) == x,] %>% pull(binomial)
  t <- t[t$binomial %in% e,]
  return(t)
}) %>% do.call('rbind',.)


# difference in number of hbunits covered
# > nrow(fish) - nrow(fish_filtered)
# [1] 339476 #~1% of occurrence records
# based on updated exotic and extinct
# > nrow(fish) - nrow(fish_filtered)
# [1] 355872
# == 16,396 extra records removed due to extinctions

# merge with hb12_p data to have info on MAIN_BAS
# fish_end <- fish %>%
#   full_join(.,hb12_p %>% as_tibble() %>% select(-geom)) %>%
#   group_by(binomial) %>%
#   summarize(endemism = as.integer(length(unique(MAIN_BAS)) == 1))

#' join table with the fish data
fish <- fish_filtered %>% 
  # full_join(.,fish_end) %>%
  filter(HYBAS_ID %in% catch_hb$HYBAS_ID) %>%
  inner_join(.,catch_hb)
# total species from 11422 to 8507

fish_exo <- fish_exotic %>% 
  # full_join(.,fish_end) %>%
  filter(HYBAS_ID %in% catch_hb$HYBAS_ID) %>%
  inner_join(.,catch_hb)



#' calculate species richness per catchment
fish_sr <- fish %>%
  group_by(ID) %>%
  summarise(
    SR_tot = length(unique(binomial))
  )

fish_sr_exo <- fish_exo %>%
  group_by(ID) %>%
  summarise(
    SR_exo = length(unique(binomial))
  )

fish_sr <- left_join(fish_sr,fish_sr_exo)
fish_sr$SR_exo[is.na(fish_sr$SR_exo)] <- 0
# range of SR values
summary(fish_sr)

# save table
write.csv(fish_sr,'proc/pcrglobwb_5min_basins_SR.csv',row.names = F)

# bind back with catch and plot
fish_sr_sp <- inner_join(catch %>% mutate(ID = as.character(ID)),fish_sr)
fish_sr_sp$log_sr <- log10(fish_sr_sp$SR_tot)
plot(fish_sr_sp["log_sr"])

