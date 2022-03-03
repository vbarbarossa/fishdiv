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

tab <- read.csv('tabs/input_tab.csv')
tab.t <- read.csv('tabs/input_tab_transformed.csv')

# # first fit with Q only
# (fit <- lmer('SR_tot ~ Q_MEAN + (1 + 1|BAS)',data = tab.t))
# MuMIn::r.squaredGLMM(fit)
# 
# plot(residuals(fit),tab.t$POP)
# df <- tab.t %>% mutate(resid = residuals(fit))
# 
# (fit2 <- lm('resid ~ POP + DAMS + URB + CROP_PRES',data = df)) %>% summary
# MuMIn::r.squaredGLMM(fit2)

# modify lmer call for dredging
lmer.glmulti<-function(formula,data,random="",...) {
  newf <- formula
  newf[[3]] <- substitute(f+r,
                          list(f=newf[[3]],
                               r=reformulate(random)[[2]]))
  lme4::lmer(newf,data=data,
             REML=FALSE,...)
}


# define variables
covariates_selection <- tab.t %>% select(-BAS,-KG, -REALM, -ID,-SR_tot, -SR_exo,-starts_with('Q_M'), -starts_with('Q_DOY')) %>% colnames

response_selection = 'SR_tot'
random_term <- 'BAS'
# interaction_term <- c('TEMP_PRES','ELEVATION') # interactions with Q_magnitude variables
interaction_term <- c('HFP2009','FSI') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')


res_filtered <- foreach(Qvar = Q_magnitude) %do% {
  
  d <- foreach(nv = 1:(length(covariates_selection)+length(interaction_term)+1),.combine = 'rbind') %do% read.csv(paste0('tabs/dredging/dredge_coefficients_no',nv,'_',response_selection,'_',Qvar,'.csv'))
  
  # filter out rows with correlated terms
  rows_to_filter <- numeric()
  for(j in 1:nrow(d)){
    # if(sum(as.integer(is.na(d[j,c('ELEVATION','SLOPE')]))) == 0) rows_to_filter <- c(rows_to_filter,j)
    if(sum(as.integer(is.na(d[j,c('TEMP_PRES','TEMP_DELTA')]))) == 0) rows_to_filter <- c(rows_to_filter,j)
    if(sum(as.integer(is.na(d[j,c('Q_MIN','Q_CV')]))) == 0) rows_to_filter <- c(rows_to_filter,j)
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
  
  # write.csv(dfilt,paste0('tabs/dredge_coefficients_',response_selection,'_',Qvar,'_FILTERED.csv'),row.names = F)
  
  return(dfilt)
  
}

# no_pred_sel <- c(10,10,11)
fit <- list()
for(i in 1:3){
  Qvar <- Q_magnitude[i]
  
  mod <- paste0('SR_tot ~',
                res_filtered[[i]] %>% filter(BIC == min(BIC)) %>% dplyr::select(model) %>% pull %>% as.character,
                ' + ',
                paste0("(1|",random_term,")"))
  mod <- gsub(Q_magnitude[i],'Q',mod)
  
  df <- tab.t
  colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))
  
  fit[[i]] <- lmer(mod, data = df)
  
}

mod <- fit[[1]]

library(redres)
# plot(mod)
# plot_redres(mod, type = 'std_mar')

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
ggsave('figs/residuals_Qmean.jpg',width = 10, height = 10)

p2 <- ggResidpanel::resid_panel(mod,plots='all')
ggsave('figs/residuals_diagnostics_Qmean.jpg',p2,width = 10, height = 10)

# refit residuals based on unstransformed vars --------------------------------------
i=1
Qvar <- Q_magnitude[i]
mod <- paste0('SR_tot ~',
              res_filtered[[i]] %>% filter(BIC == min(BIC)) %>% dplyr::select(model) %>% pull %>% as.character,
              ' + ',
              paste0("(1|",random_term,")"))
mod <- gsub(Q_magnitude[i],'Q',mod)

df <- tab
colnames(df) <- gsub(Q_magnitude[i],'Q',colnames(df))

mod <- lmer(mod, data = df)

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
ggsave('figs/residuals_untransformed_Qmean.jpg',p,width = 10, height = 10)

p2 <- ggResidpanel::resid_panel(mod,plots='all')
ggsave('figs/residuals_diagnostics_untransformed_Qmean.jpg',p2,width = 10, height = 10)

