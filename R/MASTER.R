
# folders and files ----------------------------------------------------
# GSIM main directory
GSIM_dir <- 'data/GSIM/'

# HYDROSHEDS BASINS DIR
BAS_dir <- '~/surfdrive/data/watersheds_hybas12_rds/'

# thresholds -----------------------------------------------------------
# minimum number of years the station has been monitored
MIN_MONITORING_YEARS <- 10 #-

# hydrobasins level 8 (mostly used to draw the fish ranges)
# lapply(list.files('F:/hydroBASINS/dbf_lev8/',full.names = T),foreign::read.dbf) %>% do.call('rbind',.) %>% summary()
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0    223.4    475.7    708.0    952.0 115569.9 
MIN_UP_AREA <- 500 #km2

# maximum allowed overlap between catchment boundaries of stations 
# (therefore belonging to the same main basin)
MAX_OVERLAP_PERC <- 50