# set wd
# knitr::opts_knit$set(root.dir = '~/projects/fishdiv')

#' load assigned variables from MASTER.R
source('R/MASTER.R')

# what has been loaded
as.list(.GlobalEnv)

#' load packages needed
valerioUtils::libinv(c('dplyr','sf'))


#' *****************************************************************************************************************************
#' # Compile monitoring stations data
#' 
#' 
#' ## Flow data
#' 
#' 
#' Exclude stations with less than `r MIN_MONITORING_YEARS` years of monitoring.
#' 
#' Reliable yearly values are selected following Gudmundsson et al., 2018:
#' *Index values at a yearly time step are reliable if at least 350 daily observations are declared reliable.*

yr_th <- 10
# use the homogeneity table to get a list of stations with MIN_MONITORING_YEARS or more good years of record
stations_selection <- read.csv(paste0(GSIM_dir,'GSIM_indices/HOMOGENEITY/yearly_homogeneity.csv')) %>%
  as_tibble() %>%
  filter(number.good.time.steps >= yr_th) %>%
  select(gsim.no,ends_with('step'),ends_with('steps'))

# scan through these stations and compute long-term mean of indices provided
start_time <- Sys.time() #monitor time elapsed
stations_indices <- lapply( # use lapply function that can be easily parallelized on linux
  pull(stations_selection,gsim.no), # get the vector of names from the tibble
  function(x){
    tab <- read.csv(paste0(GSIM_dir,'GSIM_indices/TIMESERIES/yearly/',x,'.year'),skip = 21) %>% #skip first 21 rows which are not tabulated
      dplyr::filter(n.available >= 350) %>%
      mutate(year = lapply(date,function(x) strsplit(x,'-')[[1]][1]) %>% as.numeric) %>% 
      # filter(year >= 1950) %>%
      select(-date,-n.missing,-n.available,-year) %>%
      #some files have been encoded as both csv and tab-delim
      #therefore I need to clean the strings from extra '\t' separators and read them again as numeric
      mutate_all(function(y) gsub('[\r\n\t]', '',y)) %>%
      #some records have n.available = 365 but 0 values everywhere
      #need to remove those. They have NA in CV (0/0) --> use that
      filter(CV != 'NA') %>% # no need to use is.na() as at this stage they are still all characters
      # and then convert to numeric to avoid warning messages
      mutate_all(as.numeric)
    
    # if(nrow(tab) >= yr_th){
    return(
      bind_cols(data.frame(gsim.id = x,no_years = nrow(tab)),as.data.frame(t(colMeans(tab))) )
    ) 
    # }
  }
) %>% do.call('rbind',.) %>%
  right_join(stations_selection,., by = c('gsim.no' = 'gsim.id'))
Sys.time() - start_time #~6 min

# # check total no of years
# sum(stations_indices$no_years[stations_indices$gsim.no %in% tab$ID]) # tab is final tab

stations_indices <- stations_indices %>% electe(-no_years)

#' ## Physical attributes
#' 

stations_attributes <- read.csv(paste0(GSIM_dir,'GSIM_catalog/GSIM_catchment_characteristics.csv')) %>%
  as_tibble() %>%
  filter(gsim.no %in% stations_indices$gsim.no) %>%
  filter(!is.na(long.org)) %>% #make sure there are no missing values from the originla coordinates (but this is the case)
  # create new coordinate columns without missing values
  # paste the new lon, lat and where these are missing vals paste the original coordinates
  # this is needed to later plot the data points
  # since GSIM provides an estimated catchment boundary also for those stations 
  # that have not beed realocated by the automated routine
  mutate(LON = ifelse(is.na(long.new),long.org,long.new),
         LAT = ifelse(is.na(lat.new),lat.org,lat.new))

#' visualize the stations distribution

plot(st_geometry(st_as_sf(stations_attributes,coords = c('LON','LAT'),crs=4326)),
     pch = 21, cex = 0.5, main = paste0('n = ',nrow(stations_attributes)))

#' exclude stations with an upstream catchment area < `r MIN_UP_AREA` km2,
#' and check the plot again

stations_attributes <- stations_attributes %>% filter(area.est > MIN_UP_AREA)

plot(st_geometry(st_as_sf(stations_attributes,coords = c('LON','LAT'),crs=4326)),
     pch = 21, cex = 0.5, main = paste0('n = ',nrow(stations_attributes %>% filter(area.est > MIN_UP_AREA))))

#' exclude stations whose catchment boundary overlaps for more than `r MAX_OVERLAP_PERC`%

crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") #%>% st_transform(crs_custom)

# assign a main basin value to each station
start_time <- Sys.time() #monitor time elapsed
basins <- lapply(
  c('af','ar','as','au','eu','gr','na','sa','si'), 
  function(x){
    readRDS(paste0(BAS_dir,x,'.rds')) %>%
      st_buffer(0) %>%
      st_crop(bb) %>%
      st_transform(crs_custom) %>%
      mutate(MAIN_BAS_AREA = as.numeric(st_area(.)/10**6)) %>% #include area in km2
      # make the shapefile lighter, anyway there are no stations catchments < than the threshold
      # do minus 100 km2 since there might be differences in the area calculation depending on the projection used
      filter(MAIN_BAS_AREA >= (MIN_UP_AREA-100)) 
  }
) %>%  do.call('rbind',.)
Sys.time() - start_time #~1.5 min

# get the MAIN_BAS id for each station by intersecting the polygons with the points
start_time <- Sys.time()
stations_filter1 <- st_intersection(
  st_as_sf(stations_attributes,coords = c('LON','LAT'),crs=4326) %>% st_transform(crs_custom),
  basins %>% st_buffer(0) # buffer to avoid TopologyException error
) %>% as_tibble() %>% select(-geometry)
Sys.time() - start_time # ~6 min


# core function that compile the catchments for the stations based on the maximum area overlap criteria
start_time <- Sys.time()
stations_catchments <- lapply(
  # from the main basin with more stations
  table(stations_filter1$MAIN_BAS %>% as.character) %>% sort(.,decreasing = T) %>% names,
  function(ws){
    
    # sub-table with only the needed stations belonging to the main basin
    t <- stations_filter1 %>% filter(MAIN_BAS == ws) %>% select(gsim.no,area.est,MAIN_BAS)
    
    if(nrow(t) == 1){
      # simply return the catchment shpefile
      catch <- read_sf(paste0(GSIM_dir,'GSIM_catchments/',tolower(as.character(t$gsim.no)),'.shp'))%>% 
        mutate(gsim.no = t$gsim.no) %>% select(gsim.no) %>%
        full_join(.,t) %>%
        st_buffer(0) %>% # buffer for any further geometry operation needed
        select(gsim.no,area.est,MAIN_BAS)
      
      return(catch)
      
    }else{
      
      # compile a shapefile with all the catchments
      catch <- lapply(
        as.character(t$gsim.no), function(x) read_sf(paste0(GSIM_dir,'GSIM_catchments/',tolower(x),'.shp'))
      ) %>% do.call('rbind',.) %>% mutate(gsim.no = t$gsim.no) %>% select(gsim.no) %>%
        full_join(.,t) %>%
        # order from most to least downstream (using catchment area)
        arrange(desc(area.est)) %>%
        st_buffer(0) %>% 
        select(gsim.no,area.est,MAIN_BAS)
      
      # set index for iterative routine
      i = 1
      
      # iterate from most downstream to most upstream station
      # and exclude stations based on catchment overlap criteria
      while(!is.na(catch$gsim.no[i])){ # this stops at the end of the rows of the catch shapefile 
        
        # calculate difference polygon of most downstream catchment with upstream catchments
        d <- st_difference(catch[i,],catch[i+1:nrow(catch),]) %>% 
          # transform to projected coordinates to calculate area
          st_transform(crs_custom) %>% 
          # calculate perc_overlap with intersecting points as (1 - area_diff/area)*100 = (area_up/area)*100
          mutate(perc_overlap = (1 - as.numeric(st_area(.)/10**6)/area.est)*100) %>% 
          # filter out stations with overlap > MAX_OVERLAP_PERC
          filter(perc_overlap <= MAX_OVERLAP_PERC)
        
        # update the catchment shapefile accordingly
        catch <- catch %>%
          filter(gsim.no %in% c(as.character(gsim.no[1:i]),as.character(d$gsim.no.1))) %>%
          arrange(desc(area.est))
        
        # update index
        i = i + 1
      }
      
      return(catch)
    }
    
  }
) %>% do.call('rbind',.)
Sys.time() - start_time #~ 26 min

#' exclude stations falling within areas with less than X% species coverage
#' <span style="color:red"> **TBD** </span>

#' save important layers

# tables
write.csv(stations_indices %>% filter(gsim.no %in% stations_catchments$gsim.no),'tabs/stations_indices.csv',row.names = F)
write.csv(stations_attributes %>% filter(gsim.no %in% stations_catchments$gsim.no),'tabs/stations_attributes.csv',row.names = F)

# shapefiles
st_write(stations_attributes %>% filter(gsim.no %in% stations_catchments$gsim.no) %>% 
           st_as_sf(.,coords = c('LON','LAT'),crs=4326) %>% select(gsim.no,area.est,quality),
         'spatial/stations_points.gpkg'
)
st_write(stations_catchments,
         'spatial/stations_catchments.gpkg'
)
st_write(basins %>% filter(MAIN_BAS %in% stations_catchments$MAIN_BAS) %>% st_transform(4326),
         'spatial/main_basins.gpkg'
)

# clear the console
rm(list = ls())

#' *****************************************************************************************************************************
#' # Number of dams from GOOD2 dataset
#' 

# # read dams data
# dams <- read_sf('data/GOOD2/GOOD2_dams.shp')
# 
# # read GSIM catchments
# catch <- read_sf('spatial/stations_catchments.gpkg')
# 
# # check crs
# st_crs(dams) == st_crs(catch)
# 
# # get the sparse list of intersections
# lst <- st_intersects(catch,dams)
# 
# # number of catchments without dams:
# length(lst) - lst %>% sapply(.,function(x) length(x) != 0) %>% sum
# 
# # ids of dams
# dams_id <- dams$DAM_ID
# 
# # compile table with HYBAS_ID and corresponding gsim.no
# names(lst) <- catch$gsim.no
# catch_dams <- lst[which(unname(lst %>% sapply(.,function(x) length(x) != 0)))] %>% 
#   lapply(seq_along(.), function(n,v,i){data.frame(DAMS_ID = dams_id[v[[i]]],gsim.no = n[i])},n=names(.),v = .) %>%
#   do.call('rbind',.) %>%
#   as_tibble()
# 
# # count number of dams in each catchment
# dams_no <- catch_dams %>%
#   group_by(gsim.no) %>%
#   summarise(dams_good2_no = n())
# 
# # join to the catchments table
# catch_new <- left_join(catch,dams_no)
# # set to 0 the NAs
# catch_new$dams_good2_no[is.na(catch_new$dams_good2_no)] <- 0
# 
# # and update the previous catchment gpkg file
# st_write(catch_new,'spatial/stations_catchments2.gpkg',delete_dsn = TRUE)

#' ********************************************************************************
#' # Sample climate data
#' 
source('R/MASTER.R')
library(raster)
# read GSIM catchments
catch <- read_sf('spatial/stations_catchments.gpkg')

# check crs
st_crs(catch)

# sample worldclim layers
lyr = c('spatial/temp_hist.tif','spatial/temp_cur.tif','spatial/prec_hist.tif','spatial/prec_cur.tif')
name = c('t_hist','t_cur','p_hist','p_cur')
df <- data.frame(gsim.no = catch$gsim.no)
for(i in 1:length(lyr)){
  r <- raster(lyr[i]) %>% projectRaster(raster(crs = 4326))
  df[,name[i]] <- extract(r,catch,fun=mean,na.rm=T)
}

# save sampled data
write.csv(df,'tabs/stations_climate.csv',row.names = F)

# clear the console
rm(list = ls())

#' *****************************************************************************************************************************
#' # Species richness data sampling
#' 
#' 
#' **Q: should we exclude species occurring in endoreic basins?**

source('R/MASTER.R')

# what has been loaded
as.list(.GlobalEnv)

#' load packages needed
valerioUtils::libinv(c('dplyr','sf'))

#' intersect the catchments with the HB12 points to assign the gsim id to each catchment

# load hydrobasins 12 shpefile points layer
hb12_p <- read_sf('~/surfdrive/Documents/projects/fishsuit/data/hybas12_points_nolakes.gpkg') %>%
  dplyr::select(HYBAS_ID,MAIN_BAS)

# load the previously produced catchments layer
catch <- read_sf('spatial/stations_catchments.gpkg')

# check crs
st_crs(hb12_p) == st_crs(catch)

# get the sparse list of intersections
lst <- st_intersects(catch,hb12_p)

# number of catchments without hydrobasins:
length(lst) - (lst %>% sapply(.,function(x) length(x) != 0) %>% sum)

# ids of hydrobasins
hybas_id <- hb12_p$HYBAS_ID

# compile table with HYBAS_ID and corresponding gsim.no
names(lst) <- catch$gsim.no
catch_hb <- lst[which(unname(lst %>% sapply(.,function(x) length(x) != 0)))] %>% 
  lapply(seq_along(.), function(n,v,i){data.frame(HYBAS_ID = hybas_id[v[[i]]],gsim.no = n[i])},n=names(.),v = .) %>%
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
  group_by(gsim.no) %>%
  summarize(
    SR_tot = length(unique(binomial))
    # ,SR_end = length(unique(binomial[endemism == 1]))
  )

fish_sr_exo <- fish_exo %>%
  group_by(gsim.no) %>%
  summarize(
    SR_exo = length(unique(binomial))
  )

fish_sr <- left_join(fish_sr,fish_sr_exo)
fish_sr$SR_exo[is.na(fish_sr$SR_exo)] <- 0
# range of SR values
summary(fish_sr)

# save table
write.csv(fish_sr,'tabs/stations_SR.csv',row.names = F)

# clear the console
rm(list = ls())


#' ***************************************************************************
#' # Paleoconnectivity
#' 

# get matrix of values
M <- read.csv('data/palaeoconnection_Dias et al 2014_adj.csv',sep=',',header = T) %>%
  apply(1,function(x){
    s <- strsplit(x,',')[[1]]
    s[s == "\"1\""] <- 1
    s[s == "NA"] <- 0
    return(s[-c(1,2)] %>% as.numeric)
  })

# and vector of names
names <- read.csv('data/palaeoconnection_Dias et al 2014_adj.csv',sep=',',header = T) %>%
  apply(1,function(x){
    strsplit(x,',')[[1]][1]
  })

# compute the area of the basins based on ted dataset
bas <- read_sf('~/surfdrive/Documents/projects/fishsuit/data/Tedesco/Basin042017_3119.shp') %>% mutate(id = 1:nrow(.)) %>%
  select(BasinName,Area=Surf_area) %>% as_tibble() %>% select(-geometry) %>%
  filter(BasinName %in% names) %>%
  as.data.frame()
row.names(bas) <- bas$BasinName

# filter out extra rows and col
to_exclude <- which(!names %in% bas$BasinName)
M <- M[-to_exclude,-to_exclude]

# area vector
v <- bas[names[-to_exclude],'Area']

# multiply the matrix by the vector to get the total sum of each row
tot <- data.frame(basin = names[-to_exclude],
                  paleoarea_extra = M %*% v)

# merge back with hydrobasins/gsim data ---
# load hydrobasins 12 shpefile points layer
library(raster)
hb12_p <- read_sf('~/surfdrive/Documents/projects/fishsuit/data/hybas12_points_nolakes.gpkg') %>%
  dplyr::select(HYBAS_ID,MAIN_BAS)

# load the previously produced catchments layer
catch <- read_sf('spatial/stations_catchments.gpkg')

# check crs
st_crs(hb12_p) == st_crs(catch)

# get the sparse list of intersections
lst <- st_intersects(catch,hb12_p)

# number of catchments without hydrobasins:
length(lst) - (lst %>% sapply(.,function(x) length(x) != 0) %>% sum)

# ids of hydrobasins
hybas_id <- hb12_p$HYBAS_ID

# compile table with HYBAS_ID and corresponding gsim.no
names(lst) <- catch$gsim.no
catch_hb <- lst[which(unname(lst %>% sapply(.,function(x) length(x) != 0)))] %>% 
  lapply(seq_along(.), function(n,v,i){data.frame(HYBAS_ID = hybas_id[v[[i]]],gsim.no = n[i])},n=names(.),v = .) %>%
  do.call('rbind',.) %>%
  as_tibble()

bas <- read_sf('~/surfdrive/Documents/projects/fishsuit/data/Tedesco/Basin042017_3119.shp') %>% mutate(id = 1:nrow(.)) 
bas_ras <-  fasterize::fasterize(sf = bas,raster = raster::raster(res = 1/12), field = 'id')
hb12_p$ws <- raster::extract(bas_ras,hb12_p)
# check most recurring ted basin id per main bas, and assign that value to all the sub-basins belonging to the main bas
hb12_main_bas <- hb12_p %>%
  as_tibble() %>% dplyr::select(-geom) %>%
  filter(!is.na(ws)) %>%
  group_by(MAIN_BAS) %>%
  summarise(
    bas = table(ws) %>% sort(decreasing = T) %>% .[1] %>% names
  )

# merge back with hb12_p
hb12 <- hb12_p %>% as_tibble() %>% dplyr::select(-geom) %>% left_join(hb12_main_bas) %>% dplyr::select(-ws)

cb <- left_join(catch_hb,hb12) %>%
  dplyr::select(gsim.no,bas) %>%
  distinct() %>%
  left_join(bas %>% mutate(bas = as.character(id)) %>% dplyr::select(bas,basin = BasinName)) %>%
  dplyr::select(-geometry)

res <- inner_join(cb,tot)

write.csv(res,'tabs/stations_paleoarea.csv',row.names = F)

