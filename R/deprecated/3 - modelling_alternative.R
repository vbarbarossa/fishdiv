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

#' assign general model name
mod_name <- 'VIFsel_manualTrans'

#' set VIF threshold for model selection
VIF_th <- 5

# ------------------------------------------------------------------------------
# CUSTOM FUNCTIONS -------------------------------------------------------------

# function to plot and print table of coeffs
plot_and_print <- function(fit_list,name, coefs_names = NULL){
  
  if(!dir.exists(name)) dir.create(name,recursive = T)
  
  library(jtools); library(ggplot2)
  
  # histograms
  ph1 <- tab %>% 
    select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
    gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_minimal()
  ph2 <- tab.t %>% 
    select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
    gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_minimal()
  
  ggsave(paste0(name,'/hist.jpg'),ph1,width = 200,height = 200,units = 'mm')
  ggsave(paste0(name,'/hist_trans.jpg'),ph2,width = 200,height = 200,units = 'mm')
  
  
  # bivariate plots
  pb1 <- tab %>%
    select(-SR_exo, -KG, -REALM,-BAS,-ID) %>%
    gather(-SR_tot, key = "var", value = "value") %>% 
    ggplot(aes(x = value, y = SR_tot)) +
    geom_point(alpha = 0.2) +
    facet_wrap(~ var, scales = "free",ncol = 4) +
    theme_bw()
  pb2 <- tab.t %>%
    select(-SR_exo, -KG, -REALM,-BAS,-ID) %>%
    gather(-SR_tot, key = "var", value = "value") %>% 
    ggplot(aes(x = value, y = SR_tot)) +
    geom_point(alpha = 0.2) +
    facet_wrap(~ var, scales = "free",ncol = 4) +
    theme_bw()
  ggsave(paste0(name,'/bivariate.jpg'),pb1,width = 200,height = 200,units = 'mm')
  ggsave(paste0(name,'/bivariate_trans.jpg'),pb2,width = 200,height = 200,units = 'mm')
  
  # regression coeffs
  p <- plot_summs(fit_list[[2]],fit_list[[1]],fit_list[[3]],
                  model.names = c('Min. discharge','Mean discharge','Max. discharge'),
                  colors = colorspace::sequential_hcl(4,'Viridis')[1:3],
                  coefs = coefs_names
  ) +
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
  ggsave(paste0(name,'/coefficients_regression.jpg'), p,width = 150, height = 150, units = 'mm', dpi = 600)
  
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
    for(i in colnames(coefficients(mod)[[1]])[-1]){
      dfr <- rbind(
        dfr,
        data.frame(res = compute_redres(mod),
                   var = mod@frame[,i],
                   name = i)
      )
    }
    
    library(ggplot2)
    p <- ggplot(dfr) +
      geom_point(aes(x = var, y = res)) +
      facet_wrap('name',ncol = 4,scales = 'free_x') +
      theme_bw()
    ggsave(paste0(name,'/residuals_',Qv,'.jpg'),width = 10, height = 10)
    
    p2 <- ggResidpanel::resid_panel(mod,plots='all')
    ggsave(paste0(name,'/residuals_diagnostics_',Qv,'.jpg'),p2,width = 10, height = 10)
    
  }
  
}

# VIF selection function
vif_func<-function(in_frame,thresh=10,trace=T,...){
  library(fmsb)
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    in_dat<-in_frame
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      vif_vals<-NULL
      var_names <- names(in_dat)
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      vif_max<-as.numeric(vif_vals[max_row,2])
      if(vif_max<thresh) break
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
    }
    return(names(in_dat))
  }
}

# run the models for different VIF thresholds
for(VIF_th in c(2.5,5,7.5)){
  
  # ------------------------------------------------------------------------------
  # TRANSFORMATIONS --------------------------------------------------------------
  
  # load dataset
  tab <- read.csv('tabs/input_tab.csv')
  
  # check histogram
  tab %>% 
    select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
    gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_minimal()
  
  # transform highly skewed variables
  tab.t <- tab
  tab.t$SR_tot <- log10(tab$SR_tot)
  tab.t$AREA <- log10(tab$AREA)
  tab.t$Q_MAX <- log10(tab$Q_MAX)
  tab.t$Q_MIN <- log10(tab$Q_MIN+min(tab$Q_MIN[tab$Q_MIN>0])) # add min to correct for zeroes
  tab.t$Q_MEAN <- log10(tab$Q_MEAN)
  # tab.t$PALEO_AREA <- log10(tab$PALEO_AREA + min(tab$PALEO_AREA[tab$PALEO_AREA != 0]))
  # tab.t$ELEVATION <- log10(tab$ELEVATION)
  # tab.t$FSI <- log10(tab$FSI)
  
  # check histogram
  tab.t %>% 
    select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
    gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_minimal()
  
  # scale to 0 mean and 1 std
  tab.t <- tab.t %>% select_if(is.numeric) %>% apply(2,scale) %>% as.data.frame()
  tab.t$KG <- tab$KG
  tab.t$REALM <- tab$REALM
  tab.t$BAS <- tab$BAS
  tab.t$ID <- tab$ID
  
  
  # ------------------------------------------------------------------------------
  # PREDICTORS SELECTION ---------------------------------------------------------
  
  # define variables
  covariates_selection <- tab.t %>% select(-BAS,-KG, -REALM, -ID,-SR_tot, -SR_exo,-starts_with('Q_M')) %>% colnames
  
  # run VIF selection
  covariates_selection <- vif_func(in_frame = tab.t[,covariates_selection], thresh = VIF_th)
  
  
  # ------------------------------------------------------------------------------
  # FITTING ----------------------------------------------------------------------
  
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
    
    m <- lmer(paste(
      response_selection,"~", # response
      paste(c('Q',covariates_selection),collapse=" + "), # fixed terms
      '+',
      paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
      '+',
      paste0("(1|",random_term,")") # random term
    ),
    data = df)
    
    fit[[i]] <- m
    
  }
  
  # check common vars
  ll <- character()
  for(j in 1:3) ll <- c(ll,coefficients(fit[[j]])[[1]] %>% colnames())
  ll %>% unique
  
  # # define coefficients correspondence table
  # cc <- c(
  #   # streamflow
  #   "Discharge" = "Q","Discharge seasonality" = "Q_CV",
  #   
  #   # anthropogenic
  #   "Human Footprint Index (HFI)" = "HFP2009", "HFI*Discharge" = "I(HFP2009 * Q)",
  #   "Fragmentation Status Index (FSI)" = "FSI", "FSI*Discharge" = "I(FSI * Q)",
  #   
  #   # habitat area, heterogeneity and isolation
  #   "Catchment area" = "AREA", 
  #   "Elevation" = "ELEVATION",
  #   "Topographic Index" = "TI", 
  #   
  #   # climate
  #   "Precipitation" = "PREC_PRES", 
  #   # "Temperature" = "TEMP_PRES",
  #   "Latitude" = "LAT",
  #   
  #   # quaternary climate stability
  #   "Precipitation change" = "PREC_DELTA", "Temperature change" = "TEMP_DELTA",
  #   "Paleo area" = "PALEO_AREA"
  #   
  #   
  # )
  
  # plot diagnostics
  plot_and_print(fit_list = fit, name = paste0('proc/fullMod_vif',VIF_th,'_transManual'))
  
}

# do the same for automatically transformed variables
for(VIF_th in c(2.5,5,7.5)){
  
  # ------------------------------------------------------------------------------
  # TRANSFORMATIONS --------------------------------------------------------------
  
  # load dataset
  tab <- read.csv('tabs/input_tab.csv')
  tab.t <- read.csv('tabs/input_tab_transformed.csv')
  
  # ------------------------------------------------------------------------------
  # PREDICTORS SELECTION ---------------------------------------------------------
  
  # define variables
  covariates_selection <- tab.t %>% select(-BAS,-KG, -REALM, -ID,-SR_tot, -SR_exo,-starts_with('Q_M')) %>% colnames
  
  # run VIF selection
  covariates_selection <- vif_func(in_frame = tab.t[,covariates_selection], thresh = VIF_th)
  
  
  # ------------------------------------------------------------------------------
  # FITTING ----------------------------------------------------------------------
  
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
    
    m <- lmer(paste(
      response_selection,"~", # response
      paste(c('Q',covariates_selection),collapse=" + "), # fixed terms
      '+',
      paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
      '+',
      paste0("(1|",random_term,")") # random term
    ),
    data = df)
    
    fit[[i]] <- m
    
  }
  
  # check common vars
  ll <- character()
  for(j in 1:3) ll <- c(ll,coefficients(fit[[j]])[[1]] %>% colnames())
  ll %>% unique
  
  # plot diagnostics
  plot_and_print(fit_list = fit, name = paste0('proc/fullMod_vif',VIF_th,'_transAuto'))
  
}




# ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
# ------------------------------------------------------------------------------
# FINAL MODEL ------------------------------------------------------------------
# ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

# MANUAL ----------------------
# adjust labels for final model
# VIF_th <- 10
# load dataset
tab <- read.csv('tabs/input_tab.csv')

# check histogram
tab %>% 
  select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal()

# transform highly skewed variables
tab.t <- tab
tab.t$SR_tot <- log10(tab$SR_tot)
tab.t$AREA <- log10(tab$AREA)
tab.t$Q_MAX <- log10(tab$Q_MAX)
tab.t$Q_MIN <- log10(tab$Q_MIN+min(tab$Q_MIN[tab$Q_MIN>0])) # add min to correct for zeroes
tab.t$Q_MEAN <- log10(tab$Q_MEAN)
tab.t$Q_CV <- log10(tab$Q_CV)
tab.t$ELEVATION <- log10(tab$ELEVATION)
tab.t$PALEO_AREA <- (tab$AREA-tab$PALEO_AREA-tab$AREA)/(tab$AREA+tab$PALEO_AREA)
tab$PREC_PRES <- log10(tab$PREC_PRES)
tab.t$LAT <- abs(tab$LAT)
# tab.t$FSI <- log10(tab$FSI)

# check histogram
tab.t %>% 
  select_if(is.numeric) %>% select(-BAS,-KG,-REALM,-SR_exo) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal()

# scale to 0 mean and 1 std
tab.t <- tab.t %>% select_if(is.numeric) %>% apply(2,scale) %>% as.data.frame()
tab.t$KG <- tab$KG
tab.t$REALM <- tab$REALM
tab.t$BAS <- tab$BAS
tab.t$ID <- tab$ID

# check bivariate correlations among covariates
jpeg('figs/corrplot_manualTrans.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cor(tab.t %>% select(-c('KG','REALM','BAS','ID')), method = 'spearman', use = "pairwise.complete.obs"), method = 'number', type = 'lower',number.cex = 0.8)
dev.off()
# BASED ON CORRPLOT:
# take out LAT and TEMP_DELTA since they are highly correlated with TEMP_PRES, and based on theory
# in the modelling phase, exclude Q_CV in combination with Q_MIN, because they are highly correlated

# ------------------------------------------------------------------------------
# PREDICTORS SELECTION ---------------------------------------------------------

# define variables
covariates_selection <- tab.t %>% 
  select(-BAS,-KG, -REALM, -ID,-SR_tot, -SR_exo,-starts_with('Q_M'), -LAT, -TEMP_DELTA) %>% colnames

# run VIF selection
# covariates_selection <- vif_func(in_frame = tab.t[,c(covariates_selection)], thresh = 10)


# ------------------------------------------------------------------------------
# FITTING ----------------------------------------------------------------------

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
  
  if(Qvar == 'Q_MIN') covariates_selection <- covariates_selection[-which(covariates_selection == 'Q_CV')]


  m <- lmer(paste(
    response_selection,"~", # response
    paste(c('Q',covariates_selection),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ),
  data = df)

  fit[[i]] <- m

}

# check common vars
ll <- character()
for(j in 1:3) ll <- c(ll,coefficients(fit[[j]])[[1]] %>% colnames())
ll %>% unique

# define coefficients correspondence table
cc <- c(
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
  # "Latitude" = "LAT",

  # quaternary climate stability
  "Precipitation change" = "PREC_DELTA", 
  # "Temperature change" = "TEMP_DELTA",
  "Paleo area" = "PALEO_AREA"


)

# plot diagnostics
# plot_and_print(fit_list = fit, name = paste0('proc/finalMod_vif',VIF_th,'_transManual'),coefs_names = cc)
plot_and_print(fit_list = fit, name = paste0('proc/finalMod_transManual_aprioriSel'),coefs_names = cc)


# AUTO ------------------------
# adjust labels for final model
VIF_th <- 7.5
# load dataset
tab <- read.csv('tabs/input_tab.csv')
tab.t <- read.csv('tabs/input_tab_transformed.csv')

# ------------------------------------------------------------------------------
# PREDICTORS SELECTION ---------------------------------------------------------

# define variables
covariates_selection <- tab.t %>% select(-BAS,-KG, -REALM, -ID,-SR_tot, -SR_exo,-starts_with('Q_M')) %>% colnames

# run VIF selection
covariates_selection <- vif_func(in_frame = tab.t[,covariates_selection], thresh = VIF_th)


# ------------------------------------------------------------------------------
# FITTING ----------------------------------------------------------------------

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
  
  m <- lmer(paste(
    response_selection,"~", # response
    paste(c('Q',covariates_selection),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ),
  data = df)
  
  fit[[i]] <- m
  
}

# check common vars
ll <- character()
for(j in 1:3) ll <- c(ll,coefficients(fit[[j]])[[1]] %>% colnames())
ll %>% unique

# define coefficients correspondence table
cc <- c(
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
  "Precipitation change" = "PREC_DELTA", "Temperature change" = "TEMP_DELTA",
  "Paleo area" = "PALEO_AREA"
  
  
)

# plot diagnostics
plot_and_print(fit_list = fit, name = paste0('proc/finalMod_vif',VIF_th,'_transAuto'),coefs_names = cc)


