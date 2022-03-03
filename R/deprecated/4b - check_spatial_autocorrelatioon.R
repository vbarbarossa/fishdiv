
library(dplyr); library(sf);

# preprocess data to get residuals at centroid of catchments -------------------

# refit best models
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')
random_term <- 'BAS'
response_selection = 'SR_tot'
df <- read.csv('tabs/input_tab_transformed.csv')

res_filtered <- lapply(Q_magnitude,function(x) read.csv(paste0('tabs/dredge_coefficients_',response_selection,'_',x,'_FILTERED.csv')))

# calculate residuals
i = 1
Qvar <- Q_magnitude[i]
mod <- paste0('SR_tot ~',
              res_filtered[[i]] %>% filter(BIC == min(BIC)) %>% dplyr::select(model) %>% pull %>% as.character,
              ' + ',
              paste0("(1|",random_term,")"))
mod <- gsub(Q_magnitude[i],'Q',mod)
colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
fit <- lme4::lmer(mod, data = df)
df$resid <- resid(fit)

# link to spatial data for lat long
s <- read_sf('spatial/stations_catchments.gpkg')
s <- s %>% select(ID = gsim.no) %>% right_join(df[,c('ID','resid')])

# calculate centroids of s to calculate distances
sc <- st_centroid(s %>% st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# inverse of distance among centroids for Moran's I
dist <- st_distance(sc,sc)
dist.inv <- 1/dist
diag(dist.inv) <- 0
attr(dist.inv,"units") <- NULL
class(dist.inv) <- setdiff(class(dist.inv),"units")

# check Moran's I --------------------------------------------------------------
mi <- ape::Moran.I(sc$resid, dist.inv)
mi

# check correlogram ------------------------------------------------------------

d <- sc

t0=Sys.time()
ncf.cor <- ncf::correlog(st_coordinates(d)[,1] %>% as.numeric(), st_coordinates(d)[,2] %>% as.numeric(), d$resid,
                         increment=100, resamp=1)
Sys.time()-t0

library(ggplot2)
dp <- as.data.frame(sp.cor)
p <- ggplot(data.frame(distance = dp$dist.class,
                       correlation = dp$coef,
                       pval = ifelse(dp$p.value >= 0.05, 'not-significant', 'significant'))) +
  geom_point(aes(x = distance, y = correlation, color = pval)) +
  theme_minimal()
p  

p <- ggplot(data.frame(distance = unname(ncf.cor$mean.of.class),
                       correlation = unname(ncf.cor$correlation))) +
  geom_point(aes(x = distance, y = correlation),alpha = 0.1, fill = 'grey70',shape = 21, size = 2, stroke = 0) +
  theme_minimal()
# p  
ggsave('figs/residuals_correlogram_Qmean.jpg',p,width = 10, height = 5)

# try alternative method (no resampling needed)
sp.cor <- pgirmess::correlog(coords = st_coordinates(sc),z = sc$resid)
a <- plot(sp.cor,title='') # saved manually

# check semivariogram ----------------------------------------------------------
library(geoR)
dists <- dist(st_coordinates(sc)) 
summary(dists)
br = seq(0,max(dists), l = 100)

variog(coords = st_coordinates(sc),data = sc$resid, breaks = br) %>% plot
