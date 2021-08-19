# set wd
knitr::opts_knit$set(root.dir = 'D:/projects/SDRs')

#' load assigned variables from MASTER.R
source('R/MASTER.R')

# what has been loaded
as.list(.GlobalEnv)

#' load packages needed
valerioUtils::libinv(c('dplyr','tidyr','ggplot2','foreach'))

# function to normalize variables
normalize_vars <- function(t,figs_name = 'figs/histograms', BN_name = 'proc/BN'){
  
  valerioUtils::libinv(c('ggplot2','purrr','tidyr'))
  
  
  # names of the variables that need to be transformed to reach a more normalized distribution
  variables_to_transform <- t %>% colnames
  
  library(bestNormalize)
  
  BN <- list()
  t.t <- foreach(i = 1:length(variables_to_transform),.combine = 'cbind') %do% {
    
    BN[[i]] <- bestNormalize(as.data.frame(t)[,variables_to_transform[i]],allow_orderNorm = FALSE,standardize = TRUE)
    d <- as.data.frame(BN[[i]]$x.t)
    colnames(d) <- variables_to_transform[i] 
    return( d )
    
  }
  names(BN) <- variables_to_transform
  
  p <- t %>% gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_minimal()
  
  p.t <- t.t %>%
    gather() %>% 
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_minimal()
  
  if(!is.na(figs_name)) ggsave(paste0(figs_name,'.jpg'), p, scale = 2)
  
  if(!is.na(figs_name)) ggsave(paste0(figs_name,'_normalized.jpg'), p.t, scale = 2)
  
  saveRDS(BN,paste0(BN_name,'.rds'))
  
  return(as_tibble(t.t))
  
}


#' # Table formatting
#' 
#' Select the variables used in the modelling
#' 
#' First, check the correlation among the flow indices provided by GSIM

#' put together an overall table with attributes, flow metrics and SR
tab_all <- full_join(
  read.csv('tabs/stations_attributes.csv') %>% as_tibble(),
  read.csv('tabs/stations_indices.csv') %>% as_tibble()
) %>%
  inner_join(.,sf::read_sf('spatial/stations_catchments.gpkg') %>%
               as_tibble() %>% dplyr::select(-geom,-area.est)) %>%
  full_join(.,read.csv('tabs/stations_SR.csv'),by='gsim.no') %>%
  as_tibble()

write.csv(tab_all,'tabs/stations_all_attributes.csv',row.names = F)

# compute the correlation matrix for flow indices

# raw
flow_ind <- tab_all %>% select(colnames(tab_all)[which(colnames(tab_all) == 'MEAN'):which(colnames(tab_all) == 'DOYMAX7')])
# normalized
flow_ind.t <- normalize_vars(flow_ind,figs_name = 'figs/flow_indices_hist', BN_name = 'proc/flow_indices_BN')

# corr matrix
jpeg('figs/corrplot_initial_flow_indices.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cor(flow_ind, use = "pairwise.complete.obs"), method = 'number', type = 'lower',number.cex = 0.8)
dev.off()

jpeg('figs/corrplot_initial_flow_indices_normalized.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cor(flow_ind.t, use = "pairwise.complete.obs"), method = 'number', type = 'lower',number.cex = 0.8)
dev.off()

# # check dams GOOD2 vs GRanD
# plot(log10(tab_all$dams_good2_no),log10(tab_all$no.dams))
# abline(0,1) 

# interesting, some catchments have more grand dams than good2
# true for
# sum((tab_all$dams_good2_no-tab_all$no.dams) < 0) # catchments

#' Select flow indices and additional covariates based on informed choice and theoreical reasoning (see paper)
#' 
#' **Discharge covariates**
#' 
#' -	Flow variability: IQR/p50 
#' 
#' -	Minimum flow: MIN7 
#' 
#' -	Maximum flow: MAX7 
#' 
#' -	Mean (to link up with previous studies and for applied purposes) 
#' 
#' 
#' **Additional covariates**
#' 
#' *Ecosystem productivity:*
#' 
#' - climate zone (major Koppen-Geiger climate zone A-E) 
#' 
#' - latitude 
#' 
#' *Evolutionary diversification potential:*
#' 
#' -	Basin area (log-transformed) 
#' 
#' -	Mean topographic index of the catchment ln(CA/tan(slope_gradient)). Ranges from -1 to xx where high values indicate small catchments with high terrain roughness.
#' 
#' -	Elevation heterogeneity: use (q3-q1)/q2
#' 
#' -	Drainage network density
#' 
#' -	Quaternary climate stability

(tab <- tab_all %>%
    # variables that need to be computed
    mutate(
      PREC_DELTA = abs(prec_cur_mean - prec_hist_mean)/prec_cur_mean,
      TEMP_DELTA = abs(temp_cur_mean - temp_hist_mean)/temp_cur_mean,
      CROP_2015 = cropland_2015_sum/cropland_2015_count,
      CROP_1992 = cropland_1992_sum/cropland_1992_count,
      URB = nl.mean*area.est
    ) %>%
    mutate(
      CROP_PRES = (CROP_1992+CROP_2015)/2 * area.est
    ) %>%
    # select the variables and rename them
    select(
      # ID variables
      ID = gsim.no, 
      BAS = MAIN_BAS, 
      QUALITY = quality,
      
      # Discharge covariates
      Q_MEAN = MEAN, 
      Q_MIN = MIN7, 
      Q_MAX = MAX7, 
      Q_CV = CV, 
      Q_DOYMIN = DOYMIN7, 
      Q_DOYMAX = DOYMAX7,
      
      # Ecosystem productivity
      PREC_PRES = prec_cur_mean, 
      TEMP_PRES = temp_cur_mean, 
      
      # Quaternary climate stabolity
      PREC_DELTA, 
      TEMP_DELTA,
      
      # Habitat area and heterogeneity
      AREA = area.est, 
      TI = tp.mean, 
      ELEVATION = ele.mean, 
      SLOPE = slp.mean,
      
      # Anthropogenic
      POP = pop.count, 
      DAMS = dams_good2_no, 
      URB,
      CROP_PRES,
      HFP2009 = HFP2009_mean,
      HFP1993 = HFP1993_mean,
      CSI = csi,
      # Response
      SR_tot
    ) %>%
    mutate(
      FSI = round(100 - CSI,4)
    ))

#' kick-out quality level 'caution' (meaning catchment area estimate is uncertain)
(tab <- tab %>%
    filter(QUALITY !="Caution"))

#' check for NAs and flow/precipitation metrics < 0 and remove those records

# NAs
apply(tab,2,function(x) sum(is.na(x)))

# looks like URB has NAs instead of zeroes
tab$URB[is.na(tab$URB)] <- 0

# Q < 0
apply(tab %>% select(starts_with('Q'),-QUALITY,starts_with('PREC')),2,function(x) sum(x < 0,na.rm=T))

# adjust NA and Q<0 values
(tab <- tab %>%
    tidyr::drop_na() %>% # we are also dropping 351 basins with NAs in DRAINAGE
    filter(Q_MIN >= 0))

#'  save the final table
write.csv(tab,'tabs/stations_filtered.csv',row.names = F)

# make a table with densities instead ########################################################################################

(tab <- tab_all %>%
   # variables that need to be computed
   mutate(
     PREC_DELTA = abs(prec_cur_mean - prec_hist_mean)/prec_cur_mean,
     TEMP_DELTA = abs(temp_cur_mean - temp_hist_mean)/temp_cur_mean,
     CROP_2015 = cropland_2015_sum/cropland_2015_count,
     CROP_1992 = cropland_1992_sum/cropland_1992_count,
     URB = nl.mean*area.est
   ) %>%
   mutate(
     CROP_PRES = (CROP_1992+CROP_2015)/2
   ) %>%
   # select the variables and rename them
   select(
     # ID variables
     ID = gsim.no, 
     BAS = MAIN_BAS, 
     QUALITY = quality,
     
     # Discharge covariates
     Q_MEAN = MEAN, 
     Q_MIN = MIN7, 
     Q_MAX = MAX7, 
     Q_CV = CV, 
     Q_DOYMIN = DOYMIN7, 
     Q_DOYMAX = DOYMAX7,
     
     # Ecosystem productivity
     PREC_PRES = prec_cur_mean, 
     TEMP_PRES = temp_cur_mean, 
     
     # Quaternary climate stabolity
     PREC_DELTA, 
     TEMP_DELTA,
     
     # Habitat area and heterogeneity
     AREA = area.est, 
     TI = tp.mean, 
     ELEVATION = ele.mean, 
     SLOPE = slp.mean,
     
     # Anthropogenic
     POP = pop.count, 
     DAMS = dams_good2_no, 
     URB,
     CROP_PRES,
     HFP2009 = HFP2009_mean,
     HFP1993 = HFP1993_mean,
     CSI = csi,
     # Response
     SR_tot
   ) %>%
   mutate(
     FSI = round(100 - CSI,4),
     POP = POP/AREA,
     DAMS = DAMS/AREA,
     URB = URB/AREA,
     CROP_PRES = CROP_PRES/AREA
   ))

#' kick-out quality level 'caution' (meaning catchment area estimate is uncertain)
(tab <- tab %>%
    filter(QUALITY !="Caution"))

#' check for NAs and flow/precipitation metrics < 0 and remove those records

# NAs
apply(tab,2,function(x) sum(is.na(x)))

# looks like URB has NAs instead of zeroes
tab$URB[is.na(tab$URB)] <- 0

# Q < 0
apply(tab %>% select(starts_with('Q'),-QUALITY,starts_with('PREC')),2,function(x) sum(x < 0,na.rm=T))

# adjust NA and Q<0 values
(tab <- tab %>%
    tidyr::drop_na() %>% # we are also dropping 351 basins with NAs in DRAINAGE
    filter(Q_MIN >= 0))

#'  save the final table
write.csv(tab,'tabs/stations_filtered_divAREA.csv',row.names = F)

#------------------------------------------------------------------------------------------------------------------------------
#' ## NORMALIZE VARIABLES
#' 

source("R/HighstatLibV10.R") # For VIFs

#'  load the preprocessed table with all attributes
(tab <- read.csv('tabs/stations_filtered_divAREA.csv') %>%
    as_tibble() %>%
    select(-QUALITY)  )

#' filter and convert to absolute values, e.g., total cropland area etc.

# make a table with mean and in brackets min - max
(tab_ranges <- apply(tab %>% select_if(is.numeric),2,
                     function(x) paste0(round(mean(x,na.rm=T),2),' (',round(min(x,na.rm=T),2),' - ',round(max(x,na.rm=T),2), ')') 
) %>% as.data.frame)
write.csv(tab_ranges,'tabs/covariates_range_values_divAREA.csv',row.names = T)

#'***********************************************************************************************

# check correlation covariates
# raw
covariates <- tab %>%  
  select(-ID,-BAS, -Q_DOYMAX, -Q_DOYMIN, -DAMS, -HFP1993, -POP, -URB, -CROP_PRES, -SR_tot, -CSI) %>%
  rename(
    # streamflow
    "Flow mean" = "Q_MEAN", "Flow max" = "Q_MAX","Flow min" = "Q_MIN","Flow seasonality" = "Q_CV",
    # habitat area, heterogeneity and isolation
    "Catchment area" = "AREA", "Topographic Index" = "TI", "Elevation" = "ELEVATION","Slope" = "SLOPE",
    # climate
    "Precipitation" = "PREC_PRES","Temperature" = "TEMP_PRES", 
    # quaternary climate stability
    "Precipitation change" = "PREC_DELTA", "Temperature change" = "TEMP_DELTA",
    # anthropogenic
    "Human Footprint Index" = "HFP2009", 
    "Fragmentation Status Index" = "FSI"
  )

# normalized
covariates.t <- normalize_vars(covariates,figs_name = 'figs/covariates_hist', BN_name = 'proc/covariates_BN')
#predict(readRDS('proc/covariates_BN.rds')$CSI,v_est2,inverse = T)
# corr matrix
jpeg('figs/corrplot_covariates.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cor(covariates, use = "pairwise.complete.obs"), method = 'number', type = 'lower',number.cex = 0.8)
dev.off()

jpeg('figs/corrplot_covariates_normalized.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cor(covariates.t, use = "pairwise.complete.obs"), method = 'number', type = 'lower',number.cex = 0.8)
dev.off()

# 
# # Then VIFS
# corvif(tab.t %>% select(-starts_with('SR'), -HFP1993))
# 
# corvif(tab.t %>% select(-starts_with('SR'),-HFP1993,-Q_MIN, -Q_MEAN))
# 
# corvif(tab.t %>% select(-starts_with('SR'),-HFP1993,-Q_MAX, -Q_MEAN))
# 
# corvif(tab.t %>% select(-starts_with('SR'),-HFP1993,-Q_MIN, -Q_MAX))
# 
# corvif(tab.t %>% select(-starts_with('SR'),-HFP1993,-Q_MAX, -Q_MIN, -PREC_PRES))

covariates <- tab %>%  
  select(-ID,-BAS, -Q_DOYMAX, -Q_DOYMIN, -DAMS, -HFP1993, -POP, -URB, -CROP_PRES, -SR_tot, -CSI) %>%
  keep(is.numeric)
# normalized
covariates.t <- normalize_vars(covariates,figs_name = NA, BN_name = 'proc/covariates_BN')


covariates$SR_tot <- tab$SR_tot
covariates$BAS <- tab$BAS
covariates$ID <- tab$ID

covariates.t$SR_tot <- log10(tab$SR_tot) #use log10 transf for SR
covariates.t$BAS <- tab$BAS
covariates.t$ID <- tab$ID

write.csv(covariates,'tabs/input_tab_divAREA.csv',row.names=F)
write.csv(covariates.t,'tabs/input_tab_transformed_divAREA.csv',row.names=F)

