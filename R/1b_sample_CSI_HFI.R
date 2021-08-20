library(sf); library(dplyr)

#' ***************************************************************************
#' # Sample CSI

bas <- read_sf('spatial/stations_catchments.gpkg')

csi <- read_sf('data/ffr_data.gpkg')

#' intersect the catchments with CSI
# check crs
st_crs(bas) == st_crs(csi)

# get the sparse list of intersections
lst <- st_intersects(bas,csi)
names(lst) <- bas$gsim.no

# number of catchments without csi
length(lst) - (lst %>% sapply(.,function(x) length(x) != 0) %>% sum)
# exclude them
to_exclude <- which(unname(lst %>% sapply(.,function(x) length(x) == 0)))

# table for quicker reading
csi_tab <- csi %>% as_tibble() %>% select(-geom)

# compile CSI
bas <- bas[-to_exclude,]
bas$csi <- lst[-to_exclude] %>% 
  lapply(., function(x){
    
    t <- csi_tab[x,]
    
    val = sum(t$VOLUME_TCM*t$CSI)/sum(t$VOLUME_TCM)
    
    return(val)
    
    }) %>%
  do.call('c',.)

# save sampled data
write.csv(bas %>% as_tibble() %>% select(-geom) %>% select(gsim.no,csi),'tabs/stations_CSI.csv',row.names = F)

#' ***************************************************************************
#' # Sample HFI

library(raster); library(sf); library(dplyr)

# read GSIM catchments
catch <- read_sf('spatial/stations_catchments.gpkg')

# check crs
st_crs(catch)

# sample worldclim layers
lyr = c('data/HFI/Maps/HFP2009.tif','data/HFI/Maps/HFP1993.tif')
name = c('HFP2009_mean','HFP1993_mean')
df <- data.frame(gsim.no = catch$gsim.no)
for(i in 1:length(lyr)){
  r <- raster(lyr[i])
  poly <- catch %>% st_transform(st_crs(r))
  df[,name[i]] <- exactextractr::exact_extract(r,poly,fun='mean')
}

# save sampled data
write.csv(df,'tabs/stations_HFI.csv',row.names = F)


