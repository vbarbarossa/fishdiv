# set wd
knitr::opts_knit$set(root.dir = 'C:/Users/barbarossav/Documents/projects/SDRs')

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

#'  load the preprocessed table with all attributes
(tab <- read.csv('tabs/stations_filtered.csv') %>%
  as_tibble() %>%
    select(-QUALITY,-SR_end)  )

#' filter and convert to absolute values, e.g., total cropland area etc.

# make a table with mean and in brackets min - max
(tab_ranges <- apply(tab %>% select_if(is.numeric),2,
       function(x) paste0(round(mean(x,na.rm=T),2),' (',round(min(x,na.rm=T),2),' - ',round(max(x,na.rm=T),2), ')') 
      ) %>% as.data.frame)
write.csv(tab_ranges,'tabs/covariates_range_values.csv',row.names = T)

#'***********************************************************************************************
#' # Modelling
#' 
#' ## Some diagnostics first
#' 
#' ### Variables distribution 
#'  
#' Check distributions 
valerioUtils::libinv(c('ggplot2','purrr','tidyr'))
tab %>%
  select(-ID,-BAS) %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal()

# names of the variables that need to be transformed to reach a more normalized distribution
variables_to_transform <- tab %>% select(-ID,-BAS) %>% keep(is.numeric) %>% colnames

library(bestNormalize)

BN <- list()
tab.t <- foreach(i = 1:length(variables_to_transform),.combine = 'cbind') %do% {
  
  BN[[i]] <- list()
  BN[[i]][[1]] <- bestNormalize(as.data.frame(tab)[,variables_to_transform[i]],allow_orderNorm = FALSE)
  
  val <- scale(BN[[i]][[1]]$x.t)
  BN[[i]][[2]] <- attr(val,"scaled:center")
  BN[[i]][[3]] <- attr(val,"scaled:scale")
  
  names(BN[[i]]) <- c('BN','center','scale')
  
  d <- as.data.frame(as.numeric(val))
  colnames(d) <- variables_to_transform[i] 
  
  return( d )
  
}
names(BN) <- variables_to_transform

# and see what happened with the transformations
BN

# and check again the distribution
tab.t %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() +
  theme_minimal()

#' ## Multicollinearity
#' 
#' First bivariate correlations
cm <- cor(tab.t %>% select(-starts_with('SR')), use = "pairwise.complete.obs")
# and visualize
corrplot::corrplot(cm, method = 'number', type = 'lower', number.cex = 0.8)

# and save
jpeg('figs/corrplot_covariates_transformed.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cm, method = 'number', type = 'lower',number.cex = 0.8)
dev.off()


# Then VIFS
corvif(tab.t %>% select(-starts_with('SR')))

corvif(tab.t %>% select(-starts_with('SR'),-Q_MIN, -Q_MEAN))

corvif(tab.t %>% select(-starts_with('SR'),-Q_MAX, -Q_MEAN))

corvif(tab.t %>% select(-starts_with('SR'),-Q_MIN, -Q_MAX))

corvif(tab.t %>% select(-starts_with('SR'),-Q_MAX, -Q_MIN, -PREC_PRES))

#' * One model per flow variables.
#'  
#' * Excluding the central tendency of DOY decreased the VIF for both DOYMIN and DOYMAX.
#'  
#' * Excluding either drainage, TI or slope lowers the VIFs for those variables as well. I would kick out drainage as it makes us lose about 300 observations in the earlier stages of the modelling

#' ## RUN MODELS
#' 
#' 
#' ### Dredging the models

# add the BAS id to the transformed variables table
tab.t$BAS <- tab$BAS

write.csv(tab,'tabs/input_tab.csv',row.names=F)
write.csv(tab.t,'tabs/input_tab_transformed.csv',row.names=F)

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
covariates_selection <- c(tab.t %>% select(-BAS,-starts_with('SR'), -starts_with('Q_M'), -starts_with('Q_DOY')) %>% colnames,
                          'I(Q_DOYMIN * PREC_PRES * TEMP_PRES)','I(Q_DOYMAX * PREC_PRES * TEMP_PRES)')
# covariates_selection <- tab.t %>% select(starts_with('Q'),-starts_with('Q_M'),TEMP_PRES,ELEVATION) %>% colnames

response_selection = 'SR_tot'
random_term <- 'BAS'
interaction_term <- c('TEMP_PRES','ELEVATION') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')

res <- foreach(Qvar = Q_magnitude) %do% {
  
  df <- tab.t
  
  fit <- lmer(paste(
    response_selection,"~", # response
    paste(c(Qvar,covariates_selection),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*',Qvar,')',collapse = ' + '), # interaction terms with Qvar
    
    '+',
    paste0("(1 + ",Qvar,"|",random_term,")") # random term
  ),
  data = df)
  
  dred <- glmulti::glmulti(formula(fit,
                                   fixed.only=TRUE),
                           random=paste0("(1|",random_term,")"),
                           data=df, method="g",
                           conseq=3, crit=BIC,
                           fitfunc=lmer.glmulti,
                           includeobjects = TRUE,
                           confsetsize=10000, # The number of models to be looked for, i.e. the size of the returned confidence set.
                           popsize = 100, # The population size for the genetic algorithm
                           mutrate = 10^-3, # The per locus (i.e. per term) mutation rate for genetic algorithm, between 0 and 1
                           sexrate = 0.1, # The rate of sexual reproduction for the genetic algorithm, between 0 and 1
                           imm = 0.3, # The rate of immigration for the genetic algorithm, between 0 and 1
                           deltaM = 2,
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
    
    coef <- fixef(dred@objects[[i]])
    
    for(j in seq_along(coef)) t[1,names(coef[j])] <- as.numeric(coef[j])
    
    row.names(t) <- NULL
    return(t)
  }
  write.csv(t_res,paste0('tabs/dredge_coefficients_ranslope_',response_selection,'_',Qvar,'.csv'),row.names = F)
  
  return(t_res)
}
names(res) <- Q_magnitude

# make a table with the best models with decreasing number of variables
res_filtered <- foreach(i = seq_along(res)) %do% {
  
  d <- res[[i]]
  
  # create count redictors column
  no_pred_tot <- ncol(d[,(which(colnames(d) == '(Intercept)')+1):ncol(d)])
  d$no_pred <- apply(d[,(which(colnames(d) == '(Intercept)')+1):ncol(d)],
                     1,function(x) no_pred_tot - sum(is.na(x)))
  
  dfilt <- lapply(split(d,d$no_pred),function(x) x[which(x$BIC == min(x$BIC)),])%>%
    do.call('rbind',.) %>%
    arrange(desc(no_pred))
  
  write.csv(dfilt,paste0('tabs/dredge_coefficients_ranslope_',response_selection,'_',names(res)[i],'_FILTERED.csv'),row.names = F)
  
  return(dfilt)
  
}




