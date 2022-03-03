# set wd
# knitr::opts_knit$set(root.dir = 'C:/Users/barbarossav/Documents/projects/SDRs')

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

source("R/HighstatLibV10.R") # For VIFs

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

tab <- read.csv('tabs/input_tab.csv') %>% split(.$BAS) %>% lapply(.,function(x) x[x$AREA == max(x$AREA),]) %>% do.call('rbind',.) %>% as_tibble()
# tab.t <- read.csv('tabs/input_tab_transformed_divAREA.csv') %>% filter(ID %in% tab$ID)

covariates <- tab %>%  
  select(-ID,-BAS, -SR_tot, -SR_exo) %>%
  rename(
    # streamflow
    "Flow mean" = "Q_MEAN", "Flow max" = "Q_MAX","Flow min" = "Q_MIN","Flow seasonality" = "Q_CV",
    # habitat area, heterogeneity and isolation
    "Catchment area" = "AREA", "Topographic Index" = "TI", "Elevation" = "ELEVATION",
    # climate
    "Precipitation" = "PREC_PRES","Temperature" = "TEMP_PRES", 
    # quaternary climate stability
    "Precipitation change" = "PREC_DELTA", "Temperature change" = "TEMP_DELTA", "Paleo area" = "PALEO_AREA",
    # anthropogenic
    "Human Footprint Index" = "HFP2009",
    "Fragmentation Status Index" = "FSI"
  )
# normalized
covariates.t <- normalize_vars(covariates,figs_name = 'figs/covariates_hist_majorbasins', BN_name = 'proc/covariates_BN_majorbasins')

# corr matrix
jpeg('figs/corrplot_covariates_majorbasins.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cor(covariates, use = "pairwise.complete.obs"), method = 'number', type = 'lower',number.cex = 0.8)
dev.off()

jpeg('figs/corrplot_covariates__majorbasins_normalized.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cor(covariates.t, use = "pairwise.complete.obs"), method = 'number', type = 'lower',number.cex = 0.8)
dev.off()

covariates <- tab %>%  
  select(-ID,-BAS, -SR_tot, -SR_exo)
# normalized
covariates.t <- normalize_vars(covariates,figs_name = NA, BN_name = 'proc/covariates_BN_majorbasins')

covariates$SR_tot <- tab$SR_tot
covariates$BAS <- tab$BAS
covariates$ID <- tab$ID

covariates.t$SR_tot <- log10(tab$SR_tot) #use log10 transf for SR
covariates.t$BAS <- tab$BAS
covariates.t$ID <- tab$ID

tab <- covariates
tab.t <- covariates.t

# define variables
covariates_selection <- tab.t %>% select(-BAS,-ID,-REALM,-KG,-starts_with('SR'), -starts_with('Q_M'), -starts_with('Q_DOY')) %>% colnames
# covariates_selection <- tab.t %>% select(starts_with('Q'),-starts_with('Q_M'),TEMP_PRES,ELEVATION) %>% colnames

response_selection = 'SR_tot'
# random_term <- 'BAS'
interaction_term <- c('HFP2009','FSI') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')

for(Qvar in Q_magnitude){
  
  df <- tab.t
  
  fit <- glm(paste(
    response_selection,"~", # response
    paste(c(Qvar,covariates_selection),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*',Qvar,')',collapse = ' + ') # interaction terms with Qvar
  ),
  data = df)
  for(nv in 1:(length(covariates_selection)+length(interaction_term)+1)){
    
    dred <- glmulti::glmulti(formula(fit),
                             data=df, method="h",
                             conseq=3, crit=BIC,
                             includeobjects = TRUE,
                             confsetsize=1000, # The number of models to be looked for, i.e. the size of the returned confidence set.
                             popsize = 100, # The population size for the genetic algorithm
                             mutrate = 10^-3, # The per locus (i.e. per term) mutation rate for genetic algorithm, between 0 and 1
                             sexrate = 0.1, # The rate of sexual reproduction for the genetic algorithm, between 0 and 1
                             imm = 0.3, # The rate of immigration for the genetic algorithm, between 0 and 1
                             deltaM = 2,
                             minsize = nv,
                             maxsize = nv,
                             # chunk = i,
                             # chunks = 8,
                             # report = FALSE,
                             plotty = FALSE,
                             intercept=TRUE, marginality=FALSE, level=1)
    
    
    coef_tab <- data.frame(
      matrix(ncol = 1 + 
               length(covariates_selection) + 
               length(interaction_term) + 
               length(Q_magnitude),
             nrow = 1))
    if(length(interaction_term) > 0){ 
      colnames(coef_tab) <- c('(Intercept)',Q_magnitude,covariates_selection,paste0('I(',interaction_term,' * ',Qvar,')')) 
    }else{
      colnames(coef_tab) <- c('(Intercept)',Q_magnitude,covariates_selection)
    }
    
    
    t_res <- foreach(i = seq_along(dred@objects),.combine = 'rbind') %do% {
      
      rsq <- MuMIn::r.squaredGLMM(dred@objects[[i]])
      
      t <- cbind(
        data.frame(
          response = dred@formulas[[i]] %>% as.character() %>% .[2],
          model = dred@formulas[[i]] %>% as.character() %>% .[3],
          BIC = dred@crits[i],
          R2_marginal = rsq[1,1],
          R2_conditional = rsq[1,2]
        ),
        coef_tab
      )
      
      coef <- coef(dred@objects[[i]])
      
      for(j in seq_along(coef)) t[1,names(coef[j])] <- as.numeric(coef[j])
      
      row.names(t) <- NULL
      return(t)
    }
    write.csv(t_res,paste0(valerioUtils::dir_('tabs/dredging_majorbasins/'),'dredge_coefficients_no',nv,'_',response_selection,'_',Qvar,'.csv'),row.names = F)
    
  }
  
}

# make a table with the best models with decreasing number of variables
# terms mutually exclusive given the high pearson's r (>= 0.7)
# Elevation - Slope
# Temp_pres - Temp_delta
# Qmin - Qcv

res_filtered <- foreach(Qvar = Q_magnitude) %do% {
  
  d <- foreach(nv = 1:(length(covariates_selection)+length(interaction_term)+1),.combine = 'rbind') %do% read.csv(paste0('tabs/dredging_majorbasins/dredge_coefficients_no',nv,'_',response_selection,'_',Qvar,'.csv'))
  
  # filter out rows with correlated terms
  rows_to_filter <- numeric()
  for(j in 1:nrow(d)){
    # if(sum(as.integer(is.na(d[j,c('ELEVATION','SLOPE')]))) == 0) rows_to_filter <- c(rows_to_filter,j)
    if(sum(as.integer(is.na(d[j,c('TEMP_PRES','TEMP_DELTA')]))) == 0) rows_to_filter <- c(rows_to_filter,j)
    # if(sum(as.integer(is.na(d[j,c('Q_MIN','Q_CV')]))) == 0) rows_to_filter <- c(rows_to_filter,j)
    if(sum(as.integer(is.na(d[j,c('AREA','Q_MEAN')]))) == 0) rows_to_filter <- c(rows_to_filter,j)
    if(sum(as.integer(is.na(d[j,c('AREA','Q_MAX')]))) == 0) rows_to_filter <- c(rows_to_filter,j)
  }
  rows_to_filter <- unique(rows_to_filter)
  
  # filter out
  d <- d[-rows_to_filter,]
  
  # create count predictors column
  no_pred_tot <- ncol(d[,(which(colnames(d) == 'X.Intercept.')+1):ncol(d)])
  d$no_pred <- apply(d[,(which(colnames(d) == 'X.Intercept.')+1):ncol(d)],
                     1,function(x) no_pred_tot - sum(is.na(x)))
  
  dfilt <- lapply(split(d,d$no_pred),function(x) x[which(x$BIC == min(x$BIC)),])%>%
    do.call('rbind',.) %>%
    arrange(desc(no_pred))
  
  write.csv(dfilt,paste0('tabs/dredge_coefficients_',response_selection,'_',Qvar,'_FILTERED_majorbasins.csv'),row.names = F)
  
  return(dfilt)
  
}

# no_pred_sel <- c(10,10,11)
fit <- list()
for(i in 1:3){
  Qvar <- Q_magnitude[i]
  
  mod <- paste0('SR_tot ~',
                res_filtered[[i]] %>% filter(BIC == min(BIC)) %>% dplyr::select(model) %>% pull %>% as.character)
  mod <- gsub(Q_magnitude[i],'Q',mod)
  
  df <- tab.t
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  fit[[i]] <- glm(mod, data = df)
  
}

# df$resid <- resid(fit)

# fit.res <- glm(resid ~ HFP2009,data = df)
# fit.res %>% summary
# 
# library(ggplot2)
# ggplot(df,aes(x = HFP2009, y = resid)) +
#   geom_point(color = 'gray60', size = 3, alpha = 0.5) +
#   stat_smooth(method = 'lm', color = 'black', size = 1.5) +
#   theme_minimal()
# 
# 
# plot(df$HFP2009,df$resid)
# plot(df$HFP2009,df$POP)

library(jtools); library(ggplot2)

# check common vars
ll <- character()
for(j in 1:3) ll <- c(ll,coefficients(fit[[j]]) %>% names())
ll %>% unique

p <- plot_summs(fit[[2]],fit[[1]],fit[[3]],
           model.names = c('Minimum flow','Mean flow','Maximum flow'),
           colors = colorspace::sequential_hcl(4,'Viridis')[1:3],
           coefs = c(
             
             # streamflow
             "Flow" = "Q","Flow seasonality" = "Q_CV",
             
             # anthropogenic
             # "No. exotic species" = "SR_exo", "Exotic*Flow" = "I(SR_exo * Q)",
             "Human Footprint Index (HFI)" = "HFP2009",
             "Fragmentation Status Index (FSI)" = "FSI", "FSI*Flow" = "I(FSI * Q)",
             
             # habitat area, heterogeneity and isolation
             "Catchment area" = "AREA", 
             "Elevation" = "ELEVATION",
             
             # climate
             "Precipitation" = "PREC_PRES", 
             "Temperature" = "TEMP_PRES",
             
             # quaternary climate stability
              "Paleo area" = "PALEO_AREA"
           )) +
  theme_bw() +
  theme(
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
ggsave('figs/coefficients_regression_majorbasins.jpg', p,width = 150, height = 150, units = 'mm', dpi = 600)

export_summs(fit,
           model.names = Q_magnitude,
           colors = colorspace::sequential_hcl(4,'Viridis')[1:3],
           coefs = c(
             # streamflow
             "Flow" = "Q","Flow seasonality" = "Q_CV",
             
             # anthropogenic
             # "No. exotic species" = "SR_exo", "Exotic*Flow" = "I(SR_exo * Q)",
             "Human Footprint Index (HFI)" = "HFP2009",
             "Fragmentation Status Index (FSI)" = "FSI", "FSI*Flow" = "I(FSI * Q)",
             
             # habitat area, heterogeneity and isolation
             "Catchment area" = "AREA", 
             "Elevation" = "ELEVATION",
             
             # climate
             "Precipitation" = "PREC_PRES", 
             "Temperature" = "TEMP_PRES",
             
             # quaternary climate stability
             "Paleo area" = "PALEO_AREA"
             
           ), to.file = "docx", file.name = 'tabs/coefficients_regression_majorbasins.docx')

# plot(tab.t$HFP2009,tab$SR_tot/(tab$AREA**1) %>% log10)
# install.packages('effects')
# library(effects)
# 
# all <- allEffects(fit[[1]])
# plot(all)


# check Moran's I


