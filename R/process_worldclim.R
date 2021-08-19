
source('R/MASTER.R')

valerioUtils::libinv(c('raster','dplyr','foreach'))
# 3 models
# cc, mr, me
# lgm
# bi1, pr

# precipitation annual mean of historical in mm/month
p_hist <- foreach(mod = c('cc', 'mr', 'me')) %do% {
  lapply(list.files('data/worldclim_v1/',pattern = paste0(mod,'lgmpr'),full.names = T),raster) %>% brick() %>% mean()
} %>% brick() %>% mean()

# precipitation annual mean of present in  mm/month
p_cur <- lapply(paste0('data/worldclim_v1/prec',1:12,'.bil'),raster) %>% brick() %>% mean()

# temperature historical in K
t_hist <- foreach(mod = c('cc', 'mr', 'me')) %do% {
  raster(paste0('data/worldclim_v1/',mod,'lgmbi1.tif'))/10 + 273.15 #rescale from worldclim and convert to Kelvin
} %>% brick() %>% mean()

# temperature present in K
t_cur <- raster('data/worldclim_v1/bio1.bil')/10 + 273.15

# save
writeRaster(p_hist,'spatial/prec_hist.tif')
writeRaster(p_cur,'spatial/prec_cur.tif')
writeRaster(t_hist,'spatial/temp_hist.tif')
writeRaster(t_cur,'spatial/temp_cur.tif')
