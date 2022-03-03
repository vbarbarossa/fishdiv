# check residuals

library(dplyr); library(foreach)

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

# function to plot and print table of coeffs
plot_and_print <- function(fit_list,name){
  
  if(!dir.exists(name)) dir.create(name,recursive = T)
  
  library(jtools); library(ggplot2)
  
  p <- plot_summs(fit_list[[2]],fit_list[[1]],fit_list[[3]],
                  model.names = c('Min. discharge','Mean discharge','Max. discharge'),
                  colors = colorspace::sequential_hcl(4,'Viridis')[1:3],
                  coefs = c(
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
               coefs = c(
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


tab <- read.csv('tabs/input_tab.csv')
tab.t <- read.csv('tabs/input_tab_transformed.csv')

# check bivariate plots
p <- tab %>%
  select(-SR_exo, -KG, -REALM,-BAS,-ID) %>%
  gather(-SR_tot, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = SR_tot)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~ var, scales = "free",ncol = 4) +
  theme_bw()
ggsave('figs/check_bivariate.pdf',p,width = 200,height = 200,units = 'mm')

tab %>% 
  select_if(is.numeric) %>% gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal()

p <- tab.t %>%
  select(-SR_exo, -KG, -REALM,-BAS,-ID) %>%
  gather(-SR_tot, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = SR_tot)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~ var, scales = "free",ncol = 4) +
  theme_bw()
ggsave('figs/check_bivariate_normalized.pdf',p,width = 200,height = 200,units = 'mm')

# define variables
covariates_selection <- tab.t %>% select(-BAS,-KG, -REALM, -ID,-SR_tot, -SR_exo,-starts_with('Q_M'), -starts_with('Q_DOY')) %>% colnames

# run VIF selection
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

covariates_selection_vif <- vif_func(in_frame = tab.t[,covariates_selection], thresh = 5)


response_selection = 'SR_tot'
random_term <- 'BAS'#/KG
# interaction_term <- c('TEMP_PRES','ELEVATION') # interactions with Q_magnitude variables
interaction_term <- c('HFP2009','FSI') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')

# standardized and normalized --------------------------------------------------
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

plot_and_print(fit_list = fit, name = 'proc/fullmod_normalized_standardized')

# cooksd <- cooks.distance(fit[[1]])
# sample_size <- nrow(df)
# plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
# abline(h = 3*mean(cooksd), col="red")  # add cutoff line
# 
# df_ <- df[-which(cooksd > 4/sample_size),]
# m <- lmer(paste(
#   response_selection,"~", # response
#   paste(c('Q',covariates_selection),collapse=" + "), # fixed terms
#   '+',
#   paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
#   '+',
#   paste0("(1|",random_term,")") # random term
# ),
# data = df_)
# summ(m)
# ggResidpanel::resid_panel(m,plots='all')
# ggResidpanel::resid_panel(fit[[1]],plots='all')

# only Q and human vars --------------------------------------------------------
covariates_selection_alt <- c("HFP2009","FSI")
tab.s <- tab %>% select(starts_with('Q_'),SR_tot)
for(i in 1:ncol(tab.s)) tab.s[,i] <- tab.s[,i]/tab$AREA
tab.s$

fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  m <- lmer(paste(
    response_selection,"~", # response
    paste(c('Q',covariates_selection_alt),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ),
  data = df)
  
  fit[[i]] <- m
  
}

plot_and_print(fit_list = fit, name = 'proc/Q_HFI_FSI_normalized_standardized')


# raw --------------------------------------------------------------------------
fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab
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

plot_and_print(fit_list = fit, name = 'fullmod')

# custom normalization ---------------------------------------------------------
# raw --------------------------------------------------------------------------

tab %>% 
  select(covariates_selection,Q_magnitude) %>%
  select_if(is.numeric) %>% gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal()

tab.s <- tab
tab.s$SR_tot <- log10(tab$SR_tot)
tab.s$AREA <- log10(tab$AREA)
# tab.s$PALEO_AREA <- log10(tab$PALEO_AREA + min(tab$PALEO_AREA[tab$PALEO_AREA != 0]))
# tab.s$ELEVATION <- log10(tab$ELEVATION)
# tab.s$FSI <- log10(tab$FSI)
tab.s$Q_MAX <- tab.t$Q_MAX
tab.s$Q_MIN <- tab.t$Q_MIN
tab.s$Q_MEAN <- tab.t$Q_MEAN

tab.s <- tab.s %>% select_if(is.numeric) %>% apply(2,scale) %>% as.data.frame()
tab.s$KG <- tab$KG
tab.s$REALM <- tab$REALM
tab.s$BAS <- tab$BAS
tab.s$ID <- tab$ID

tab.s %>% 
  select(covariates_selection,Q_magnitude) %>%
  select_if(is.numeric) %>% gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal()


fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab.s
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


plot_and_print(fit_list = fit, name = 'fullmod_custom_trans')



# scaled only ------------------------------------------------------------------
tab.s <- tab %>% select_if(is.numeric) %>% apply(2,scale) %>% as.data.frame()
tab.s$SR_tot <- tab$SR_tot
tab.s$KG <- tab$KG
tab.s$REALM <- tab$REALM
tab.s$BAS <- tab$BAS
tab.s$ID <- tab$ID

fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab.s
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

plot_and_print(fit_list = fit, name = 'fullmod_std')

# scaled logSR -----------------------------------------------------------------
tab.s <- tab %>% select_if(is.numeric) %>% apply(2,scale) %>% as.data.frame()
tab.s$SR_tot <- log10(tab$SR_tot)
tab.s$KG <- tab$KG
tab.s$REALM <- tab$REALM
tab.s$BAS <- tab$BAS
tab.s$ID <- tab$ID

fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab.s
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

plot_and_print(fit_list = fit, name = 'fullmod_std_logSR')


# glmer with loglink ------------------------------------------------------------------
tab.s <- tab %>% select_if(is.numeric) %>% apply(2,scale) %>% as.data.frame()
tab.s$SR_tot <- tab$SR_tot
tab.s$KG <- tab$KG
tab.s$REALM <- tab$REALM
tab.s$BAS <- tab$BAS
tab.s$ID <- tab$ID

fit <- list()
for(i in 1:length(Q_magnitude)){
  Qvar = Q_magnitude[i]
  df <- tab.s
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  m <- glmer(paste(
    response_selection,"~", # response
    paste(c('Q',covariates_selection),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*','Q',')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ),
  data = df, family = 'poisson')
  
  fit[[i]] <- m
  
}

plot_and_print(fit_list = fit, name = 'fullmod_std_glmmPoisson')

