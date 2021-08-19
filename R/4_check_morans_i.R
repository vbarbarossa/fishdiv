
library(dplyr); library(sf);

Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')
random_term <- 'BAS'
response_selection = 'SR_tot'
df <- read.csv('tabs/input_tab_transformed_divAREA.csv')


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
s <- read_sf('spatial/stations_catchments2.gpkg')
s <- s %>% select(ID = gsim.no) %>% right_join(df[,c('ID','resid')])

# calculate centroids of s
sc <- st_centroid(s %>% st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# inverse of distance among centroids
dist <- st_distance(sc,sc)
dist.inv <- 1/dist
diag(dist.inv) <- 0
attr(dist.inv,"units") <- NULL
class(dist.inv) <- setdiff(class(dist.inv),"units")

# dist.inv[1:5, 1:5]

# moran's i
ape::Moran.I(sc$resid, dist.inv)

# check the variogram
# install.packages("geoR")
library(geoR)
dists <- dist(st_coordinates(sc)) 
summary(dists)
br = seq(0,max(dists), l = 500)

variog(coords = st_coordinates(sc),data = sc$resid, breaks = br) %>% plot


# # untransform data
# # can use yeojohnson from bestnormalize
# bn <- readRDS('proc/covariates_BN.rds')
# 
# # try
# xf = predict(bn$SR_tot, 
#              newdata = df$SR_tot, #bn$SR_tot$BN$chosen_transform$x.t, 
#              inverse = T)
# 
# plot(xf,bn$SR_tot$x)
# data.frame(xf,bn$SR_tot$x)
