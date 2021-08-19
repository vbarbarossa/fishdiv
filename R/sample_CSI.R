library(sf); library(dplyr)

bas <- read_sf('spatial/stations_catchments2.gpkg')

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

st_write(bas,'spatial/stations_catchments.gpkg')






