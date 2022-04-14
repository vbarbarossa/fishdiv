#' load assigned variables from MASTER.R
source('R/MASTER.R')

# what has been loaded
as.list(.GlobalEnv)

#' load packages needed
valerioUtils::libinv(c(
  'dplyr','tidyr','ggplot2','foreach',
  'caret',
  'lme4',# for mixed effect modelling
  'Hmisc',
  'glmulti'# for dredging models
))

# ------------------------------------------------------------------------------
# CUSTOM FUNCTIONS & PARAMS ----------------------------------------------------

# function to plot and print table of coeffs and other diagnostics
plot_and_print <- function(fit_list,name, coefs_names = NULL){
  
  if(!dir.exists(name)) dir.create(name,recursive = T)
  
  library(jtools); library(ggplot2)
  
  # data frames for plotting with correct names
  t <- tab %>% select(-ID,-BAS,-KG,-REALM,-SR_exo)
  tt <- tab.t %>% select(-ID,-BAS,-KG,-REALM,-SR_exo)
  
  if(!is.null(coefs_names)){
    t <- t %>% rename(c(
      'Species Richness' = 'SR_tot',
      'Discharge mean' = 'Q_MEAN',
      'Discharge min.' = 'Q_MIN',
      'Discharge max.' = 'Q_MAX',
      coefCorr[-which(coefCorr %in% c('I(FSI * Q)','I(HFP2009 * Q)','Q'))]))
    tt <- tt %>% rename(c(
      'Species Richness' = 'SR_tot',
      'Discharge mean' = 'Q_MEAN',
      'Discharge min.' = 'Q_MIN',
      'Discharge max.' = 'Q_MAX',
      coefCorr[-which(coefCorr %in% c('I(FSI * Q)','I(HFP2009 * Q)','Q'))]))
  }
  
  # histograms
  ph1 <- t %>%
    gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_minimal()
  ph2 <- tt %>%
    gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_minimal()
  
  ggsave(paste0(name,'/hist.jpg'),ph1,width = 200,height = 200,units = 'mm')
  ggsave(paste0(name,'/hist_trans.jpg'),ph2,width = 200,height = 200,units = 'mm')
  
  
  # bivariate plots
  pb1 <- t %>%
    rename(SR_tot = 'Species Richness') %>%
    gather(-SR_tot, key = "var", value = "value") %>% 
    ggplot(aes(x = value, y = SR_tot)) +
    geom_point(alpha = 0.2) +
    ylab('Species Richness') +
    facet_wrap(~ var, scales = "free",ncol = 4) +
    theme_bw()
  pb2 <- tt %>%
    rename(SR_tot = 'Species Richness') %>%
    gather(-SR_tot, key = "var", value = "value") %>% 
    ggplot(aes(x = value, y = SR_tot)) +
    geom_point(alpha = 0.2) +
    ylab('Species Richness') +
    facet_wrap(~ var, scales = "free",ncol = 4) +
    theme_bw()
  ggsave(paste0(name,'/bivariate.jpg'),pb1,width = 200,height = 200,units = 'mm')
  ggsave(paste0(name,'/bivariate_trans.jpg'),pb2,width = 200,height = 200,units = 'mm')
  
  # correlation matrix
  jpeg(paste0(name,'/corrplot.jpg'),width = 200, height = 200, res = 600, units = 'mm')
  if(!is.null(coefs_names)){
    corrplot::corrplot(cor(tt %>% select(-'Species Richness'), 
                           method = 'spearman', use = "pairwise.complete.obs"), 
                       method = 'number', type = 'lower',number.cex = 0.8)
    
  }else{
    corrplot::corrplot(cor(tt %>% select(-'SR_tot'), 
                           method = 'spearman', use = "pairwise.complete.obs"), 
                       method = 'number', type = 'lower',number.cex = 0.8)
  }
  dev.off()
  
  # regression coeffs
  p <- plot_summs(fit_list[[2]],fit_list[[1]],fit_list[[3]],
                  model.names = c('Min. discharge','Mean discharge','Max. discharge'),
                  colors = colorspace::sequential_hcl(4,'Viridis')[1:3],
                  coefs = coefs_names
  ) +
    theme_bw() +
    theme(
      text = element_text(size = 11),
      axis.title = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.direction = "horizontal",
      legend.position = "top",
      legend.title = element_blank(),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-5,-10,-10,-10)
    )
  p
  ggsave(paste0(name,'/coefficients_regression.jpg'), p,width = 150, height = 150, units = 'mm', dpi = 600)
  ggsave(paste0(name,'/coefficients_regression.pdf'), p,width = 150, height = 140, units = 'mm')
  
  export_summs(fit_list,
               model.names = Q_magnitude,
               colors = colorspace::sequential_hcl(4,'Viridis')[1:3],
               coefs = coefs_names
               , to.file = "docx", file.name = paste0(name,'/coefficients_regression.docx'))
  
  # residuals diagnostics
  for(i in 1:length(fit_list)){
    Qv <- Q_magnitude[i]
    mod <- fit_list[[i]]
    
    library(redres)
    dfr <- data.frame(res = compute_redres(mod),
                      var = predict(mod) %>% as.numeric,
                      name = 'predicted')
    for(j in colnames(coefficients(mod)[[1]])[-1]){
      dfr <- rbind(
        dfr,
        data.frame(res = compute_redres(mod),
                   var = mod@frame[,j],
                   name = j)
      )
    }
    
    for(j in 1:length(coefs_names)) dfr$name <- sub(paste0('^',coefs_names[j],'$'),names(coefs_names)[j],dfr$name)
    dfr$name <- gsub("HFP2009","HFI",dfr$name)
    
    library(ggplot2)
    p <- ggplot(dfr) +
      geom_point(aes(x = var, y = res)) +
      facet_wrap('name',ncol = 4,scales = 'free_x') +
      theme_bw()
    ggsave(paste0(name,'/residuals_',Qv,'.jpg'),width = 10, height = 10)
    
    p2 <- ggResidpanel::resid_panel(mod,plots='all')
    ggsave(paste0(name,'/residuals_diagnostics_',Qv,'.jpg'),p2,width = 10, height = 10)
    
  }
  
  # VIF tables
  for(i in 1:length(fit)){
    Qv <- Q_magnitude[i]
    write.csv(
      performance::check_collinearity(fit[[i]]) %>% as.data.frame(),
      paste0(name,'/VIFs_',Qv,'.csv'),
      row.names = F
    )
  }
  
  # goodness of fit
  for(i in 1:length(fit_list)){
    
    Qv <- Q_magnitude[i]
    m <- fit[[i]]
    
    r2_scores <- performance::performance(m)
    
    pg<-ggplot() +
      geom_ribbon(aes(x = c(0,3.5), ymax = c(0.5,4),ymin = c(-0.5,3)),fill = 'gray90') +
      geom_point(data = m@frame,aes(
        x = (SR_tot*attr(scale(log10(tab$SR_tot)),'scaled:scale')+attr(scale(log10(tab$SR_tot)),'scaled:center')), 
        y = (predict(m)*attr(scale(log10(tab$SR_tot)),'scaled:scale')+attr(scale(log10(tab$SR_tot)),'scaled:center'))
      ), 
      color = 'gray20', alpha = 0.5) +
      geom_abline(intercept = 0,slope = 1, color = 'blue') +
      geom_text(aes(x = c(3.4,3.4), y = c(0.2,0.4),label = c(paste0('R^2~marginal == ',round(r2_scores$R2_marginal,2)),
                                                             paste0('R^2~conditional == ',round(r2_scores$R2_conditional,2)))),
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
    
    ggsave(paste0(name,'/goodness_of_fit_',Qv,'.jpg'),pg,
           width = 100,height = 100,dpi = 600,units = 'mm',scale = 1)
    
  }
  
}

# define coefficients correspondence table
coefCorr <- c(
  # streamflow
  "Discharge" = "Q","Discharge seasonality" = "Q_CV",
  
  # anthropogenic
  "Human Footprint Index (HFI)" = "HFP2009", "HFI*Discharge" = "I(HFP2009 * Q)",
  "Fragmentation Status Index (FSI)" = "FSI", "FSI*Discharge" = "I(FSI * Q)",
  
  # habitat area, heterogeneity and isolation
  "Catchment area" = "AREA",
  "Elevation" = "ELEVATION",
  "Topographic Index" = "TI",
  
  # climate
  "Precipitation" = "PREC_PRES",
  "Temperature" = "TEMP_PRES",
  "Latitude" = "LAT",
  
  # quaternary climate stability
  "Precipitation change" = "PREC_DELTA", 
  "Temperature change" = "TEMP_DELTA",
  "Paleo area" = "PALEO_AREA"
  
  
)

# ------------------------------------------------------------------------------
# HISTOGRAMS AND TRANSFORMATIONS -----------------------------------------------

# load dataset
tab <- read.csv('tabs/input_tab.csv')
tab$PALEO_AREA <- (tab$AREA-tab$PALEO_AREA-tab$AREA)/(tab$AREA+tab$PALEO_AREA)
tab$LAT <- abs(tab$LAT)

# make a ranges table
(tab_ranges <- apply(tab %>% select_if(is.numeric),2,
                     function(x) paste0(round(mean(x,na.rm=T),2),' (',round(min(x,na.rm=T),2),' - ',round(max(x,na.rm=T),2), ')') 
) %>% as.data.frame)
write.csv(tab_ranges,'tabs/covariates_range_values.csv',row.names = T)

# check histogram
# tab %>% 
#   select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
#   gather() %>%
#   ggplot(aes(value)) +
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram() +
#   theme_minimal()

# transform highly skewed variables
tab.t <- tab
tab.t$SR_tot <- log10(tab$SR_tot)
tab.t$AREA <- log10(tab$AREA)
tab.t$Q_MAX <- log10(tab$Q_MAX)
tab.t$Q_MIN <- log10(tab$Q_MIN+min(tab$Q_MIN[tab$Q_MIN>0])) # add min to correct for zeroes
tab.t$Q_MEAN <- log10(tab$Q_MEAN)
tab.t$Q_CV <- log10(tab$Q_CV)
tab.t$ELEVATION <- log10(tab$ELEVATION)
tab.t$PREC_PRES <- log10(tab$PREC_PRES)

# check histogram
# tab.t %>% 
#   select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
#   gather() %>%
#   ggplot(aes(value)) +
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram() +
#   theme_minimal()

# scale to 0 mean and 1 std
tab.t <- tab.t %>% select_if(is.numeric) %>% apply(2,scale) %>% as.data.frame()

# overwrite categorical variables that do not need to be scaled
tab.t$KG <- tab$KG
tab.t$REALM <- tab$REALM
tab.t$BAS <- tab$BAS
tab.t$ID <- tab$ID

# check bivariate correlations among covariates
# jpeg('figs/corrplot_manualTrans.jpg',width = 200, height = 200, res = 600, units = 'mm')
# corrplot::corrplot(cor(tab.t %>% select(-c('KG','REALM','BAS','ID')), method = 'spearman', use = "pairwise.complete.obs"), method = 'number', type = 'lower',number.cex = 0.8)
# dev.off()
# # account for random effect in calculation of correlation --> same result between groups
# library(psych)
# stats <- statsBy(tab.t %>% select(-SR_exo,-KG,-REALM,-ID), "BAS",method = 'spearman')
# print(stats, short = FALSE)
# BASED ON CORRPLOT:
# take out LAT and TEMP_DELTA since they are highly correlated with TEMP_PRES, and based on theory
# in the modelling phase, exclude Q_CV in combination with Q_MIN, because they are highly correlated

# ------------------------------------------------------------------------------
# PREDICTORS SELECTION ---------------------------------------------------------

# define covariates and exclude variables with rho > 0.7
covariates_selection <- tab.t %>% 
  select(-BAS,-KG, -REALM, -ID,-SR_tot, -SR_exo,-starts_with('Q_M'), 
         -LAT, -TEMP_DELTA) %>% 
  colnames

# alternatively, run VIF selection
# covariates_selection <- vif_func(in_frame = tab.t[,c(covariates_selection)], thresh = 10)

# ------------------------------------------------------------------------------
# FITTING ----------------------------------------------------------------------

# fitting params
response_selection = 'SR_tot'
random_term <- 'BAS'#/KG
interaction_term <- c('HFP2009','FSI') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX') # fit three different models for each Q variable

# fit the three models and store them in a list
fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab.t
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  csel <- covariates_selection
  if(Qvar == 'Q_MIN') csel <- covariates_selection[-which(covariates_selection == 'Q_CV')]
  
  
  m <- lmer(paste(
    response_selection,"~", # response
    paste(c('Q',csel),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ),
  data = df)
  
  fit[[i]] <- m
  
}

# plot diagnostics in a dedicated folder
plot_and_print(fit_list = fit, name = paste0('proc/finalMod_transManual'),
               coefs_names = coefCorr[-which(coefCorr %in% c('LAT','TEMP_DELTA'))])

# ------------------------------------------------------------------------------
# INTERACTION PLOTS ------------------------------------------------------------

library(sjPlot)
library(sjmisc)
library(ggplot2)

# repeat for three Q variables
for(i in 1:length(fit)){
  
  Qv <- Q_magnitude[i]
  csel <- covariates_selection
  if(Qvar == 'Q_MIN') csel <- covariates_selection[-which(covariates_selection == 'Q_CV')]
  
  df <- tab.t
  dfu <- tab
  
  Qvar <- Q_magnitude[i]
  
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  colnames(dfu) <- gsub(Q_magnitude[i],'Q',colnames(dfu))
  
  
  dfu$BAS <- as.factor(dfu$BAS)
  to_scale <- colnames(dfu %>% select(-SR_tot) %>% purrr::keep(is.numeric))
  for(j in to_scale) dfu[,j] <- scale(dfu[,j]) %>% as.numeric()
  
  fit_ <- fit[[i]]
  
  var = 'FSI'
  xax = 10**seq(-1,6,2)
  yax = 10**(seq(-1.5,2,0.5))
  vec = c(0,5,15,30,66)
  p_fsi <- plot_model(fit_, 
                      type='pred',
                      terms = c(
                        paste0('Q [',paste(round((log10(xax)-attr(scale(log10(tab[,Qvar])),'scaled:center'))/attr(scale(log10(tab[,Qvar])),'scaled:scale'),1),collapse = ','),']'),
                        paste0(var,' [',paste(round((vec-attr(scale(tab[,var]),'scaled:center'))/attr(scale(tab[,var]),'scaled:scale'),1),collapse = ','),']')
                      ),
                      title = '', axis.title = c('Flow','Species richness')) + 
    scale_color_discrete(name = "FSI", labels = vec) + #paste0(vec,' (',c(0,0.25,0.50,0.75,1),')') 
    scale_x_continuous(breaks = round((log10(xax)-attr(scale(log10(tab[,Qvar])),'scaled:center'))/attr(scale(log10(tab[,Qvar])),'scaled:scale'),1), labels = formatC(xax, format = "e", digits = 0)) +
    scale_y_continuous(breaks = log10(yax), labels = formatC(yax, format = "e", digits = 0), limits = log10(range(yax))) +
    theme_bw() +
    theme(legend.direction = 'horizontal',legend.position = 'top',
          panel.grid = element_blank(),
          text = element_text(size = 14))
  # p_fsi
  
  var = 'HFP2009'
  xax = 10**seq(-1,6,2)
  yax = 10**(seq(-1.5,2,0.5))
  vec = c(0,3.5,8,12,36)
  p_hfi <- plot_model(fit_, 
                      type='pred',
                      terms = c(
                        paste0('Q [',paste(round((log10(xax)-attr(scale(log10(tab[,Qvar])),'scaled:center'))/attr(scale(log10(tab[,Qvar])),'scaled:scale'),1),collapse = ','),']'),
                        paste0(var,' [',paste(round((vec-attr(scale(tab[,var]),'scaled:center'))/attr(scale(tab[,var]),'scaled:scale'),1),collapse = ','),']')
                      ),
                      title = '', axis.title = c('Flow','Species richness')) + 
    scale_color_discrete(name = "HFI", labels = vec) + #paste0(vec,' (',c(0,0.25,0.50,0.75,1),')') 
    scale_x_continuous(breaks = round((log10(xax)-attr(scale(log10(tab[,Qvar])),'scaled:center'))/attr(scale(log10(tab[,Qvar])),'scaled:scale'),1), labels = formatC(xax, format = "e", digits = 0)) +
    scale_y_continuous(breaks = log10(yax), labels = formatC(yax, format = "e", digits = 0), limits = log10(range(yax))) +
    theme_bw() +
    theme(legend.direction = 'horizontal',legend.position = 'top',
          panel.grid = element_blank(),
          text = element_text(size = 14))
  # p_hfi
  
  library(ggpubr)
  figure <- ggarrange(p_hfi + theme(legend.title = element_blank()) + rremove("ylab") + rremove("xlab"),
                      p_fsi + theme(legend.title = element_blank()) + rremove("ylab") + rremove("xlab"),
                      labels = c("a) Human Footprint Index", "b) Fragmentation Status Index"),
                      font.label = list(size = 14, color = "black", face = "plain", family = NULL),
                      # labels.x = 0,
                      # labels.y = 5,
                      vjust = 2,
                      hjust = -0.3,
                      align = "hv", 
                      ncol = 2, nrow = 1)
  # figure
  f <- annotate_figure(figure, 
                       left = text_grob("Species richness [-]", size = 14, rot = 90),
                       bottom = text_grob(bquote('Discharge ['*m^3*'/s]'), size = 14))
  # save
  ggsave(paste0('proc/finalMod_transManual/interaction_plots_HFI_FSI',Qv,'.jpg'),f,width = 220,height = 130,dpi = 300,units = 'mm')
  ggsave(paste0('proc/finalMod_transManual/interaction_plots_HFI_FSI',Qv,'.pdf'),f,width = 220,height = 130,units = 'mm')
}


# ------------------------------------------------------------------------------
# SPATIAL RESIDUALS ------------------------------------------------------------

library(sf)
s <- read_sf('spatial/stations_catchments.gpkg')

# base params
crs_custom <- st_crs(4326)
crop_main <- c(-180, 180, -60, 90)
cropping_poly <- st_bbox(raster::extent(crop_main), crs = st_crs(4326)) %>% st_as_sfc(.)

sf_use_s2(FALSE)

# base layers
world <- read_sf('data/naturalearth/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp')[,1] %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)
bb <- read_sf('data/naturalearth/ne_110m_graticules_all/ne_110m_wgs84_bounding_box.shp') %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)
graticules <- read_sf('data/naturalearth/ne_110m_graticules_all/ne_110m_graticules_30.shp') %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)

for(i in 1:length(fit)){
  
  df <- tab.t
  Qv <- Q_magnitude[i]
  df$resid <- resid(fit[[i]])
  
  # link to spatial data for lat long
  sc <- s %>% select(ID = gsim.no) %>% right_join(df[,c('ID','resid')]) %>%
    # calculate centroids of s to calculate distances
    st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>%
    st_centroid
  
  
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
    coord_flip(xlim = c(-60,90), ylim = c(-2.2,2.2),expand = F) +
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
  ggsave(paste0('proc/finalMod_transManual/residuals_map',Qv,'.jpg'), pm,
         width = 200,height = 90,dpi = 300,units = 'mm')
  # ggsave('figs/map_resid.tiff', pm,
  #        width = 200,height = 90,dpi = 600,units = 'mm')
  
}

# ------------------------------------------------------------------------------
# SPATIAL AUTOCORRELATION ------------------------------------------------------

library(dplyr); library(sf);

# preprocess data to get residuals at centroid of catchments -------------------
s <- read_sf('spatial/stations_catchments.gpkg')

for(i in 1:length(fit)){
  
  df <- tab.t
  Qv <- Q_magnitude[i]
  df$resid <- resid(fit[[i]])
  
  # link to spatial data for lat long
  sc <- s %>% select(ID = gsim.no) %>% right_join(df[,c('ID','resid')]) %>%
    # calculate centroids of s to calculate distances
    st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>%
    st_centroid
  
  # inverse of distance among centroids for Moran's I
  dist <- st_distance(sc,sc)
  dist.inv <- 1/dist
  diag(dist.inv) <- 0
  attr(dist.inv,"units") <- NULL
  class(dist.inv) <- setdiff(class(dist.inv),"units")
  
  # check Moran's I --------------------------------------------------------------
  mi <- ape::Moran.I(sc$resid, dist.inv)
  print(mi)
  
  # check correlogram ------------------------------------------------------------
  # sp.cor <- pgirmess::correlog(coords = st_coordinates(sc),z = sc$resid)
  # jpeg(paste0('proc/finalMod_transManual/correlogram_',Qv,'_pgirmessPackage.jpg'),width = 250, height = 175, res = 600, units = 'mm')
  # plot(sp.cor,title='')
  # dev.off()
  
}

# ------------------------------------------------------------------------------
# UNTRANSFORMED ----------------------------------------------------------------
# refit the model without transforming the covariates

# fit the three models and store them in a list
fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  csel <- covariates_selection
  if(Qvar == 'Q_MIN') csel <- covariates_selection[-which(covariates_selection == 'Q_CV')]
  
  
  m <- lmer(paste(
    response_selection,"~", # response
    paste(c('Q',csel),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ),
  data = df)
  
  fit[[i]] <- m
  
}

# plot diagnostics
plot_and_print(fit_list = fit, name = paste0('proc/finalMod_untransformed'),
               coefs_names = coefCorr[-which(coefCorr %in% c('LAT','TEMP_DELTA'))])


# ------------------------------------------------------------------------------
# NESTED -----------------------------------------------------------------------
# refit the model with ID nested in REALM

random_term <- 'BAS/REALM'

# fit the three models and store them in a list
fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab.t
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  csel <- covariates_selection
  if(Qvar == 'Q_MIN') csel <- covariates_selection[-which(covariates_selection == 'Q_CV')]
  
  
  m <- lmer(paste(
    response_selection,"~", # response
    paste(c('Q',csel),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ),
  data = df)
  
  fit[[i]] <- m
  
}

# plot diagnostics
plot_and_print(fit_list = fit, name = paste0('proc/finalMod_transManual_nestedREALM'),
               coefs_names = coefCorr[-which(coefCorr %in% c('LAT','TEMP_DELTA'))])


# ------------------------------------------------------------------------------
# MAIN BASINS ------------------------------------------------------------------
# refit the model with main basins only

# major basins
tab <- read.csv('tabs/input_tab.csv') %>% split(.$BAS) %>% lapply(.,function(x) x[x$AREA == max(x$AREA),]) %>% do.call('rbind',.) %>% as_tibble()
tab$PALEO_AREA <- (tab$AREA-tab$PALEO_AREA-tab$AREA)/(tab$AREA+tab$PALEO_AREA)
tab$LAT <- abs(tab$LAT)

# check histogram
# tab %>% 
#   select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
#   gather() %>%
#   ggplot(aes(value)) +
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram() +
#   theme_minimal()

# transform highly skewed variables
tab.t <- tab
tab.t$SR_tot <- log10(tab$SR_tot)
tab.t$AREA <- log10(tab$AREA)
tab.t$Q_MAX <- log10(tab$Q_MAX)
tab.t$Q_MIN <- log10(tab$Q_MIN+min(tab$Q_MIN[tab$Q_MIN>0])) # add min to correct for zeroes
tab.t$Q_MEAN <- log10(tab$Q_MEAN)
tab.t$Q_CV <- log10(tab$Q_CV)
tab.t$ELEVATION <- log10(tab$ELEVATION)
tab.t$PREC_PRES <- log10(tab$PREC_PRES)

# check histogram
# tab.t %>% 
#   select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
#   gather() %>%
#   ggplot(aes(value)) +
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram() +
#   theme_minimal()

# scale to 0 mean and 1 std
tab.t <- tab.t %>% select_if(is.numeric) %>% apply(2,scale) %>% as.data.frame()
# overwrite variables that do not need to be scaled
tab.t$KG <- tab$KG
tab.t$REALM <- tab$REALM
tab.t$BAS <- tab$BAS
tab.t$ID <- tab$ID

# ------------------------------------------------------------------------------
# PREDICTORS SELECTION ---------------------------------------------------------

covariates_selection <- tab.t %>% 
  select(-BAS,-KG, -REALM, -ID,-SR_tot, -SR_exo,-starts_with('Q_M'), 
         -LAT, -TEMP_DELTA) %>% 
  colnames

# ------------------------------------------------------------------------------
# FITTING ----------------------------------------------------------------------

response_selection = 'SR_tot'
interaction_term <- c('HFP2009','FSI') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX') # fit three different models for each Q variable

# fit the three models and store them in a list
fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab.t
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  csel <- covariates_selection
  # if(Qvar == 'Q_MIN') csel <- covariates_selection[-which(covariates_selection == 'Q_CV')]
  if(Qvar %in% c('Q_MEAN','Q_MAX')) csel <- covariates_selection[-which(covariates_selection == 'AREA')]
  
  m <- glm(paste(
    response_selection,"~", # response
    paste(c('Q',csel),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*','Q',')',collapse = ' + ') # interaction terms with Qvar
  ),
  data = df)
  
  fit[[i]] <- m
  
}

# plot diagnostics
plot_and_print(fit_list = fit, name = paste0('proc/finalMod_transManual_majorBasins'),
               coefs_names = coefCorr[-which(coefCorr %in% c('LAT','TEMP_DELTA'))])
