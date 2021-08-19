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

tab <- read.csv('tabs/input_tab_divAREA.csv')
tab.t <- read.csv('tabs/input_tab_transformed_divAREA.csv')

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
covariates_selection <- tab.t %>% select(-BAS,-POP,-DAMS,-URB,-CROP_PRES,-HFP1993,-starts_with('SR'), -starts_with('Q_M'), -starts_with('Q_DOY')) %>% colnames
# covariates_selection <- tab.t %>% select(starts_with('Q'),-starts_with('Q_M'),TEMP_PRES,ELEVATION) %>% colnames

response_selection = 'SR_tot'
random_term <- 'BAS'
interaction_term <- c('TEMP_PRES','ELEVATION') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')

for(Qvar in Q_magnitude){
  
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
  for(nv in 1:(length(covariates_selection)+length(interaction_term)+1)){
    
    dred <- glmulti::glmulti(formula(fit,
                                     fixed.only=TRUE),
                             random=paste0("(1 + ",Qvar,"|",random_term,")"),
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
    write.csv(t_res,paste0('tabs/dredge_coefficients_ranslope_no',nv,'_',response_selection,'_',Qvar,'.csv'),row.names = F)
  
  }
  
}

# make a table with the best models with decreasing number of variables
res_filtered <- foreach(Qvar = Q_magnitude) %do% {
  
  
  d <- foreach(nv = 1:(length(covariates_selection)+length(interaction_term)+1),.combine = 'rbind') %do% read.csv(paste0('tabs/dredge_coefficients_no',nv,'_',response_selection,'_',Qvar,'.csv'))
  
  # create count redictors column
  no_pred_tot <- ncol(d[,(which(colnames(d) == 'X.Intercept.')+1):ncol(d)])
  d$no_pred <- apply(d[,(which(colnames(d) == 'X.Intercept.')+1):ncol(d)],
                     1,function(x) no_pred_tot - sum(is.na(x)))
  
  dfilt <- lapply(split(d,d$no_pred),function(x) x[which(x$BIC == min(x$BIC)),])%>%
    do.call('rbind',.) %>%
    arrange(desc(no_pred))
  
  write.csv(dfilt,paste0('tabs/dredge_coefficients_ranslope_',response_selection,'_',Qvar,'_FILTERED.csv'),row.names = F)
  
  return(dfilt)
  
}

no_pred_sel <- c(10,10,10)
i = 1
Qvar <- Q_magnitude[i]

df <- tab.t

fit <- lmer(paste0('SR_tot ~',
                   res_filtered[[i]] %>% filter(no_pred == no_pred_sel[i]) %>% select(model) %>% pull %>% as.character,
                   ' + ',
                   paste0("(1|",random_term,")")),
            data = df)

df$resid <- resid(fit)

fit.res <- glm(resid ~ HFP2009,data = df)
fit.res %>% summary

library(ggplot2)
ggplot(df,aes(x = HFP2009, y = resid)) +
  geom_point(color = 'gray60', size = 3, alpha = 0.5) +
  stat_smooth(method = 'lm', color = 'black', size = 1.5) +
  theme_minimal()


plot(df$HFP2009,df$resid)
plot(df$HFP2009,df$POP)

# plot(tab.t$HFP2009,tab$SR_tot/(tab$AREA**1) %>% log10)
# install.packages('effects')
library(effects)

all <- allEffects(fit)
plot(all)
