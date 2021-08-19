
source('R/MASTER.R')
valerioUtils::libinv(c('dplyr','tidyr','ggplot2','sf'))

#' load table with all atributes
tab <- read.csv('tabs/stations_all_attributes.csv') %>%
  as_tibble()

#' check bivariate plots
p <- tab %>%
  select(colnames(tab)[92:ncol(tab)]) %>%
  gather(-SR_tot, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = SR_tot)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~ var, scales = "free",ncol = 4) +
  theme_bw()
p

p_log <- tab %>%
  select(colnames(tab)[92:ncol(tab)]) %>%
  gather(-SR_tot, key = "var", value = "value") %>% 
  ggplot(aes(x = log10(value), y = log10(SR_tot))) +
  geom_point(alpha = 0.2) +
  facet_wrap(~ var, scales = "free",ncol = 4) +
  theme_bw()
p_log

ggsave('figs/check_bivariate.pdf',p,width = 200,height = 300,units = 'mm')
ggsave('figs/check_bivariate_all_log.pdf',p_log,width = 200,height = 300,units = 'mm')


# map all the variables used spatially

# load the table
# t <- read.csv('tabs/input_tab_transformed.csv')
# t2 <- read.csv('tabs/input_tab.csv')
t <- read.csv('tabs/input_tab_transformed_divAREA.csv')

# read spatial data and bind
s <- read_sf('spatial/stations_catchments.gpkg')

s <- s %>% select(ID = gsim.no) %>% right_join(t)

st_write(s,'spatial/input_tab.gpkg')


# library(raster)
# tc <- raster('spatial/temp_cur.tif')
# th <- raster('spatial/temp_hist.tif')
# plot(tc-th)
# plot(raster('spatial/prec_cur.tif') - raster('spatial/prec_hist.tif'))


# main centroid map-----------------------------------------------------------------------------------------
library(sf)
s <- sf::read_sf('spatial/stations_catchments.gpkg')

# calculate centroids of s
sc <- s %>% select(ID = gsim.no) %>% st_centroid()

# bind centroids to variables
# t2 <- read.csv('tabs/input_tab.csv')
t2 <- t
sc <- sc %>% right_join(t2)
s <- s %>% select(ID = gsim.no) %>% right_join(t2 %>% select(ID))


# base layers
crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)

# and draw
# p <- ggplot() +
#   geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
#   geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
#   geom_sf(data = world, fill = "grey90", lwd = NA) +
#   geom_sf(data = s, fill = "grey60", lwd = NA) +
#   geom_sf(data = sc, aes(color = log10(SR_tot))) +
#   # scale_fill_viridis_c(breaks = seq(0,1,0.1),
#   #                      labels = seq(0,1,0.1),
#   #                      limits = c(0,1),
#   #                      option = 'C',na.value = "transparent") +
#   # facet_grid(warming~scenario) +
#   theme_minimal() +
#   theme(text = element_text(size = 13),
#         panel.grid.major = element_line(color=NA),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         legend.position = 'bottom',
#         legend.key.width = unit(6,'line'),
#         strip.background = element_rect('white'),
#         strip.background.x = element_blank(),
#         strip.background.y = element_blank(),
#         strip.text = element_text(angle = 0, vjust = -1, size = 13),
#         legend.title = element_blank()
#   )
# p
# ggsave(paste0('figs/maps_RC.jpg'),p,
#        width = 183,height = 200,dpi = 600,units = 'mm')


#----------------------------------------------------------------------------------
library(dplyr); library(sf);

Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')
Q_name <- c('Mean flow','Minimum flow','Maximum flow')

random_term <- 'BAS'
response_selection = 'SR_tot'
dfu <- read.csv('tabs/input_tab_divAREA.csv')

res_filtered <- lapply(Q_magnitude,function(x) read.csv(paste0('tabs/dredge_coefficients_',response_selection,'_',x,'_FILTERED.csv')))

df_ <- list()
df_line <- list()
for(i in 1:3){
  
  df <- read.csv('tabs/input_tab_transformed_divAREA.csv')
  
  Qvar = Q_magnitude[i]
  Qn <- Q_name[i]
  
  
  mod <- paste0('SR_tot ~',
                res_filtered[[i]] %>% filter(BIC == min(BIC)) %>% dplyr::select(model) %>% pull %>% as.character,
                ' + ',
                paste0("(1|",random_term,")"))
  mod <- gsub(Q_magnitude[i],'Q',mod)
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  fit <- lme4::lmer(mod, data = df)
  df$resid <- resid(fit)
  
  # untransformed residuals
  # trans_resp <- function(x) predict(readRDS('proc/covariates_BN.rds')$SR_tot,x,inverse = T) %>% round(0) %>% log10
  # df$residu <- (dfu$SR_tot %>% log10) - trans_resp(predict(fit))
  
  library(effects)
  est<-Effect("HFP2009", partial.residuals=T, fit, quantiles = seq(0,1,by=0.02))
  est2 <- Effect("FSI", partial.residuals=T, xlevels = list(x1=c(-1,0,1)), fit)
  # est <- predictorEffect("HFP2009", partial.residuals=T, fit)
  # 
  # v_est2 <- est2$x[,1]
  # v_est2[1] <- -1.13230750292444449 #adjust decimals of first number, otherwise produces NAs --> need to check transformation
  
  df_line[[i]] <- 
    rbind(
      data.frame(
        y = est$fit,
        x = predict(readRDS('proc/covariates_BN.rds')$HFP2009,c(est$x[,1]),inverse = T),
        lo = est$lower, up = est$upper, Q = Qn, var = "HFI"
      ),
      data.frame(
        y = est2$fit,
        x = predict(readRDS('proc/covariates_BN.rds')$FSI,est2$x[,1],inverse = T),
        lo = est2$lower, up = est2$upper, Q = Qn, var = "FSI"
      )
    )
  
  df_[[i]] <- rbind(
    data.frame(resid = est$residuals + abs(min(est$residuals)),var = "HFI", Q = Qn, x = dfu$HFP2009),
    data.frame(resid = est2$residuals + abs(min(est2$residuals)),var = "FSI", Q = Qn, x = dfu$FSI)
  )
  
  
}

df_line <- do.call('rbind',df_line)
df_ <- do.call('rbind',df_)

library(ggplot2)
# p <- ggplot() +
#   geom_point(data = df_,aes(x = HFP2009, y = resid), alpha = 0.3, color = 'gray70',stroke = 0, size = 2) +
#   geom_line(data = df_line,aes(x=x,y=y)) +
#   geom_ribbon(data = df_line,aes(x = x, ymin=lo, ymax=up), linetype=2, alpha=0.1, color = 'gray20') +
#   # scale_color_manual(values = c('blue','green','red')) +
#   scale_y_continuous(breaks = c(0,1,2),labels = c(1,10,100)) +
#   xlab('Human footprint index (-)') +
#   ylab('Species richness (-)') +
#   facet_wrap('var') +
#   theme_bw() +
#   theme(strip.background = element_blank())

p <- ggplot() +
  geom_point(data = df_ %>% filter(Q == "Mean flow"),aes(x = x, y = resid), alpha = 0.3, color = 'gray70',stroke = 0, size = 2) +
  geom_line(data = df_line %>% filter(Q == "Mean flow"),aes(x=x,y=y)) +
  geom_ribbon(data = df_line %>% filter(Q == "Mean flow"),aes(x = x, ymin=lo, ymax=up), linetype=2, alpha=0.1, color = 'gray20') +
  # scale_color_manual(values = c('blue','green','red')) +
  scale_y_continuous(breaks = c(0,1,2),labels = c(1,10,100)) +
  xlab('') +
  ylab('Species richness') +
  facet_wrap('var',scales = 'free') +
  theme_bw() +
  theme(strip.background = element_blank())
p

ggsave('figs/partial_dependence.jpg',p,
       width = 80,height = 45,dpi = 600,units = 'mm',scale = 2)

#---------------------------------------------------------------------------------------------------
# RESIDUALS MAP

library(dplyr); library(sf);

Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')
random_term <- 'BAS'
response_selection = 'SR_tot'

# calculate residuals
for(i in 1:3){
  
  df <- read.csv('tabs/input_tab_transformed_divAREA.csv')
  dfu <- read.csv('tabs/input_tab_divAREA.csv')
  
  res_filtered <- lapply(Q_magnitude,function(x) read.csv(paste0('tabs/dredge_coefficients_',response_selection,'_',x,'_FILTERED.csv')))
  
  Qvar <- Q_magnitude[i]
  mod <- paste0('SR_tot ~',
                res_filtered[[i]] %>% filter(BIC == min(BIC)) %>% dplyr::select(model) %>% pull %>% as.character,
                ' + ',
                paste0("(1|",random_term,")"))
  mod <- gsub(Q_magnitude[i],'Q',mod)
  
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  fit <- lme4::lmer(mod, data = df)
  
  # predicted vs observed figure
  library(ggplot2)
  m1 <- res_filtered[[i]] %>% filter(BIC == min(BIC))
  pg<-ggplot() +
    geom_ribbon(aes(x = c(0,3.5), ymax = c(0.5,4),ymin = c(-0.5,3)),fill = 'gray90') +
    geom_point(data = df,aes(x = SR_tot, y = predict(fit)), color = 'gray20', alpha = 0.5) +
    geom_abline(intercept = 0,slope = 1, color = 'blue') +
    geom_text(aes(x = c(3.4,3.4), y = c(0.2,0.4),label = c(paste0('R^2~marginal == ',round(m1$R2_marginal,2)),
                                                           paste0('R^2~conditional == ',round(m1$R2_conditional,2)))),
              hjust = 1,parse=T) +
    # annotate(geom = 'text', x = 2, y = 0.2, label = paste('R^2~marginal == ', round(m1$R2_marginal,2)), hjust = 0, vjust = 1, parse = TRUE) +
    xlab('Observed Species Richness') +
    ylab('Predicted Species Richness') +
    scale_x_continuous(breaks = c(0,1,2,3),labels = c(1,10,100,1000)) +
    scale_y_continuous(breaks = c(0,1,2,3),labels = c(1,10,100,1000)) +
    coord_cartesian(xlim = c(-0.05,3.5),ylim = c(-0.05,3.5), expand = F) +
    theme_bw() +
    theme(panel.grid = element_blank())
  # pg
  
  ggsave(paste0('figs/goodness_of_fit_',Qvar,'.jpg'),pg,
         width = 100,height = 100,dpi = 600,units = 'mm',scale = 1)

  }






df$resid <- resid(fit)

# link to spatial data for lat long
s <- read_sf('spatial/stations_catchments2.gpkg')
s <- s %>% select(ID = gsim.no) %>% right_join(df[,c('ID','resid')])

# calculate centroids of s
sc <- st_centroid(s %>% st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# base layers
crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
crs_custom <- st_crs(4326)
crop_main <- c(-180, 180, -60, 90)
cropping_poly <- st_bbox(raster::extent(crop_main), crs = st_crs(4326)) %>% st_as_sfc(.)

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)

# and draw
library(ggplot2)
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = world, fill = "grey40", lwd = NA) +
  # geom_sf(data = s, fill = "grey60", lwd = NA) +
  geom_sf(data = sc, aes(color = resid),alpha = 0.7, shape = 19, size = 1.5, stroke = 0) +
  scale_color_gradient2(mid = 'grey90',
                        low = scales::muted('blue'),
                        high = scales::muted('red'),
                        midpoint = 0,
                        breaks = seq(-1,1,0.5),
                        labels = seq(-1,1,0.5),
                        na.value = 'transparent') +
  coord_sf(expand = F) +
  theme_minimal() +
  theme(text = element_text(size = 13),
        panel.grid.major = element_line(color=NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 13),
        legend.title = element_blank()
  )
# p

# map residuals with latitude
d <- data.frame(resid = sc$resid, lat = st_coordinates(sc %>% st_transform(4326))[,'Y'])
pl <- ggplot(d) +
  geom_point(aes(y = resid, x = lat, color = resid), alpha = 0.7, shape = 19, size = 2, stroke = 0) +
  scale_color_gradient2(mid = 'grey90',
                        low = scales::muted('blue'),
                        high = scales::muted('red'),
                        midpoint = 0,
                        breaks = seq(-1,1,0.5),
                        # labels = paste0(c('<5',seq(-5,0,1),paste0('+',seq(1,5,1)),'>5'),'%'),
                        na.value = 'transparent') +
  geom_smooth(aes(y = resid, x = lat), method = 'loess', color = 'grey40') +
  # geom_smooth(aes(y = abs(resid), x = lat), method = 'loess', color = 'grey40', lty = 2) +
  # coord_cartesian(xlim = c(-60,90)) +
  scale_x_continuous(breaks = seq(-90,90,30), labels = paste0(seq(-90,90,30),'°'),position = 'top') +
  xlab('') +
  ylab('') +
  coord_flip(xlim = c(-60,90), ylim = c(-1.2,1.2),expand = F) +
  theme_bw() +
  theme(
    text = element_text(color = 'black'),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = 'black',size=0.1),
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    legend.position = 'none'
  )
# pl
library(ggpubr)
pm <- ggarrange(
  ggarrange(p + theme(legend.position = 'none', plot.margin = unit(c(-1,0.1,-1,0.1),'cm')),
            pl + theme(plot.margin = unit(c(0.46,-0.3,-0.40,0.2),'cm')),
            nrow=1,ncol = 2,
            widths = c(1,0.3)), 
  get_legend(p),
  nrow = 2, heights = c(1,0.2)
)
ggsave('figs/map_resid.jpg', pm,
       width = 200,height = 90,dpi = 300,units = 'mm')
ggsave('figs/map_resid.tiff', pm,
       width = 200,height = 90,dpi = 600,units = 'mm')

# PLOT SR
sc <- left_join(sc,dfu %>% select(ID,SR_tot))

# and draw
library(ggplot2)
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  # geom_sf(data = s, fill = "grey60", lwd = NA) +
  geom_sf(data = sc, aes(color = log10(SR_tot)),alpha = 0.7, shape = 19, size = 1.5, stroke = 0) +
  scale_color_viridis_c(breaks = c(0,1,2,3),labels = c(1,10,100,1000)) +
  coord_sf(expand = F) +
  theme_minimal() +
  theme(text = element_text(size = 13),
        panel.grid.major = element_line(color=NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 13),
        legend.title = element_blank()
  )
p

# map residuals with latitude
d <- data.frame(resid = log10(sc$SR_tot), lat = st_coordinates(sc %>% st_transform(4326))[,'Y'])
pl <- ggplot(d) +
  geom_point(aes(y = resid, x = lat, color = resid), alpha = 0.7, shape = 19, size = 2, stroke = 0) +
  scale_color_viridis_c(breaks = c(0,1,2,3)) +
  geom_smooth(aes(y = resid, x = lat), method = 'loess', color = 'grey40') +
  # geom_smooth(aes(y = abs(resid), x = lat), method = 'loess', color = 'grey40', lty = 2) +
  # coord_cartesian(xlim = c(-60,90)) +
  scale_x_continuous(breaks = seq(-90,90,30), labels = paste0(seq(-90,90,30),'°'),position = 'top') +
  scale_y_continuous(breaks = c(0,1,2,3),labels = c(1,10,100,1000)) +
  xlab('') +
  ylab('') +
  coord_flip(xlim = c(-60,90), ylim=c(-0.1,3.4),expand = F) +
  theme_bw() +
  theme(
    text = element_text(color = 'black'),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = 'black',size=0.1),
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    legend.position = 'none'
  )
pl
library(ggpubr)
pm <- ggarrange(
  ggarrange(p + theme(legend.position = 'none', plot.margin = unit(c(-1,0.1,-1,0.1),'cm')),
            pl + theme(plot.margin = unit(c(0.46,-0.3,-0.40,0.2),'cm')),
            nrow=1,ncol = 2,
            widths = c(1,0.3)), 
  get_legend(p),
  nrow = 2, heights = c(1,0.2)
)
ggsave('figs/map_SR.jpg', pm,
       width = 200,height = 90,dpi = 300,units = 'mm')
# ggsave('figs/map_SR.tiff', pm,
#        width = 200,height = 90,dpi = 600,units = 'mm')

