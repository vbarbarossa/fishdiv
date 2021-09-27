library(sf); library(dplyr)

#' ***************************************************************************


source('R/MASTER.R')
library(raster)
# read GSIM catchments
catch <- read_sf('spatial/stations_catchments.gpkg')

# check crs
st_crs(catch)

# sample worldclim layers
lyr = c('data/Koeppen-Geiger-Classification-Reclassfied_5min_moderesampling.tif')
name = c('kg_class')
df <- data.frame(gsim.no = catch$gsim.no)
i=1
r <- raster(lyr[i]) %>% projectRaster(raster(crs = 4326))
r[] <- floor(r[]) # there are some weird values (very small percentage)

# Create the function.
getmode <- function(v) {
  uniqv <- unique(v[!is.na(v)])
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

df[,name[i]] <- extract(r,catch) %>% lapply(getmode) %>% unlist

# save sampled data
write.csv(df,'tabs/stations_KG.csv',row.names = F)

#' ***************************************************************************
# sample biomes
name = c('realm_class')

# convert ecoregions to raster file
er <- st_read('data/GIS_hs_snapped/feow_hydrosheds.shp') %>%
  left_join(read.csv('data/GIS_hs_snapped/legend.csv') %>% rename(FEOW_ID = ID))
er$Realm_ID[is.na(er$Realm_ID)] <- 9
st_crs(er) <- 4326

r <- fasterize::fasterize(er,raster(res=1/12),field = 'Realm_ID')

df <- data.frame(gsim.no = catch$gsim.no)
i=1

# Create the function.
getmode <- function(v) {
  uniqv <- unique(v[!is.na(v)])
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

df[,name[i]] <- extract(r,catch) %>% lapply(getmode) %>% unlist

# save sampled data
write.csv(df,'tabs/stations_realm.csv',row.names = F)




