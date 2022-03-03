# set wd
# knitr::opts_knit$set(root.dir = '~/projects/fishdiv')

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
# covariates_selection <- tab.t %>% select(starts_with('Q'),-starts_with('Q_M'),TEMP_PRES,ELEVATION) %>% colnames

response_selection = 'SR_tot'
random_term <- 'BAS/REALM'#/KG
# interaction_term <- c('TEMP_PRES','ELEVATION') # interactions with Q_magnitude variables
interaction_term <- c('HFP2009','FSI') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')

for(Qvar in Q_magnitude){
  
  df <- tab.t
  
  fit <- lmer(paste(
    response_selection,"~", # response
    paste(c(Qvar,covariates_selection),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*',Qvar,')',collapse = ' + '), # interaction terms with Qvar
    '+',
    # paste0('I(',interaction_term,'*AREA)',collapse = ' + '), # interaction terms with Area
    # '+',
    paste0("(1|",random_term,")") # random term
  ),
  data = df)
  for(nv in 1:(length(covariates_selection)+length(interaction_term)+1)){
    
    dred <- glmulti::glmulti(formula(fit,
                                     fixed.only=TRUE),
                             random=paste0("(1|",random_term,")"),
                             data=df, method="h",
                             conseq=3, crit=BIC,
                             fitfunc=lmer.glmulti,
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
      
      
      # rsq <- MuMIn::r.squaredGLMM(dred@objects[[i]]) #<< not working
      rsq <- piecewiseSEM::rsquared(dred@objects[[i]])[,5:6]
      
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
      
      coef <- fixef(dred@objects[[i]])
      
      for(j in seq_along(coef)) t[1,names(coef[j])] <- as.numeric(coef[j])
      
      row.names(t) <- NULL
      return(t)
    }
    write.csv(t_res,paste0(valerioUtils::dir_('tabs/dredging_nested/'),'dredge_coefficients_no',nv,'_',response_selection,'_',Qvar,'.csv'),row.names = F)
    
  }
  
}

# make a table with the best models with decreasing number of variables
# terms mutually exclusive given the high pearson's r (>= 0.7)
# Elevation - Slope
# Temp_pres - Temp_delta
# Qmin - Qcv

res_filtered <- foreach(Qvar = Q_magnitude) %do% {
  
  d <- foreach(nv = 1:(length(covariates_selection)+length(interaction_term)+1),.combine = 'rbind') %do% read.csv(paste0('tabs/dredging_nested/dredge_coefficients_no',nv,'_',response_selection,'_',Qvar,'.csv'))
  
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
  
  write.csv(dfilt,paste0('tabs/dredge_coefficients_nested_',response_selection,'_',Qvar,'_FILTERED.csv'),row.names = F)
  
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
for(j in 1:3) ll <- c(ll,coefficients(fit[[j]])[[1]] %>% colnames())
ll %>% unique

p <- plot_summs(fit[[2]],fit[[1]],fit[[3]],
                model.names = c('Minimum flow','Mean flow','Maximum flow'),
                colors = colorspace::sequential_hcl(4,'Viridis')[1:3],
                coefs = c(
                  # streamflow
                  "Flow" = "Q","Flow seasonality" = "Q_CV",
                  
                  # habitat area, heterogeneity and isolation
                  "Catchment area" = "AREA", "Topographic Index" = "TI", "Elevation" = "ELEVATION", 
                  
                  # climate
                  "Precipitation" = "PREC_PRES","Temperature" = "TEMP_PRES", 
                  
                  # quaternary climate stability
                  "Precipitation change" = "PREC_DELTA",
                  "Paleo area" = "PALEO_AREA",
                  
                  # anthropogenic
                  # "No. exotic species" = "SR_exo", "Exotic*Flow" = "I(SR_exo * Q)",
                  "Human Footprint Index (HFI)" = "HFP2009",
                  "Fragmentation Status Index (FSI)" = "FSI", "FSI*Flow" = "I(FSI * Q)"
                  
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
ggsave('figs/coefficients_regression_nested.jpg', p,width = 150, height = 150, units = 'mm', dpi = 600)

export_summs(fit,
             model.names = Q_magnitude,
             colors = colorspace::sequential_hcl(4,'Viridis')[1:3],
             coefs = c(
               # streamflow
               "Flow" = "Q","Flow seasonality" = "Q_CV",
               
               # habitat area, heterogeneity and isolation
               "Catchment area" = "AREA", "Topographic Index" = "TI", "Elevation" = "ELEVATION", 
               
               # climate
               "Precipitation" = "PREC_PRES","Temperature" = "TEMP_PRES", 
               
               # quaternary climate stability
               "Precipitation change" = "PREC_DELTA",
               "Paleo area" = "PALEO_AREA",
               
               # anthropogenic
               # "No. exotic species" = "SR_exo", "Exotic*Flow" = "I(SR_exo * Q)",
               "Human Footprint Index (HFI)" = "HFP2009",
               "Fragmentation Status Index (FSI)" = "FSI", "FSI*Flow" = "I(FSI * Q)"
               
             ), to.file = "docx", file.name = 'tabs/coefficients_regression_nested.docx')



