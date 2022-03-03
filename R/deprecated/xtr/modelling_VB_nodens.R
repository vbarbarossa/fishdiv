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

# # compute the correlation matrix 
# cm <- cor(tab_all %>% select(colnames(tab_all)[92:which(colnames(tab_all) == 'DOYMAX7')]),
#           use = "pairwise.complete.obs")
# 
# # and visualize
# corrplot::corrplot(cm, method = 'number', type = 'lower',number.cex = 0.8)
# 
# # and save
# jpeg('figs/corrplot_initial_flow_indices.jpg',width = 200, height = 200, res = 600, units = 'mm')
# corrplot::corrplot(cm, method = 'number', type = 'lower',number.cex = 0.8)
# dev.off()

#'  load the preprocessed table with all attributes
tab <- read.csv('tabs/stations_filtered.csv') %>%
  as_tibble()

#' consider two hypotheses
#' 
#' 1 - use absolute values of predictors
#' 
#' 2 - use densities

# absolutes
tab_d <- tab %>%
  mutate(
    # ELEVATION = ELEVATION/AREA,
    # SLOPE = SLOPE/AREA,
    DAMS = NO_DAMS,
    CROP_PRES = CROP_PRES*AREA,
    CROP_DELTA = CROP_DELTA*AREA,
    URB = URB*AREA
  ) %>%
  # variables not needed
  select(-QUALITY,-PREC_HIST,-TEMP_HIST,-CROP_PAST,-POP_DENS,-NO_DAMS,-SR_end)

# # # densities
# tab_d <- tab %>%
#   mutate(
#     # ELEVATION = ELEVATION/AREA,
#     # SLOPE = SLOPE/AREA,
#     DAMS = NO_DAMS/AREA,
#     Q_MAX = Q_MAX/AREA,
#     Q_MEAN = Q_MEAN/AREA,
#     Q_MIN = Q_MIN/AREA,
#     SR_tot_dens = SR_tot/AREA,
#     SR_end_dens = SR_end/AREA
#   ) %>%
#   # variables not needed
#   select(-QUALITY,-PREC_HIST,-TEMP_HIST,-CROP_PAST,-POP,-NO_DAMS,-SR_end)

#'***********************************************************************************************
#' # Modelling
#' 
#' ## Some diagnostics first
#' 
#' ### Variables distribution 
#'  
#' Check distributions 
valerioUtils::libinv(c('ggplot2','purrr','tidyr'))
tab_d %>%
  select(-ID,-BAS) %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

# names of the variables that need to be transformed to reach a more normalized distribution
variables_to_transform <- tab_d %>% select(-ID,-BAS) %>% keep(is.numeric) %>% colnames

library(bestNormalize)

BN <- list()
tab_d.t <- foreach(i = 1:length(variables_to_transform),.combine = 'cbind') %do% {
  
  BN[[i]] <- list()
  BN[[i]][[1]] <- bestNormalize(as.data.frame(tab_d)[,variables_to_transform[i]],allow_orderNorm = FALSE)
  
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
tab_d.t %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

#' ## Multicollinearity
#' 
#' First bivariate correlations
cm <- cor(tab_d.t %>% select(-starts_with('SR')), use = "pairwise.complete.obs")
# and visualize
corrplot::corrplot(cm, method = 'number', type = 'lower', number.cex = 0.8)

# Then VIFS
corvif(tab_d.t %>% select(-starts_with('SR')))

corvif(tab_d.t %>% select(-starts_with('SR'),-Q_DOYCT))

corvif(tab_d.t %>% select(-starts_with('SR'),-Q_MEAN, -Q_MAX, -Q_DOYCT, -PREC_PRES))

corvif(tab_d.t %>% select(-starts_with('SR'),-Q_MAX, -Q_MIN, -Q_DOYCT,  -PREC_PRES))

corvif(tab_d.t %>% select(-starts_with('SR'),-Q_MIN, -Q_MEAN, -Q_DOYCT, -PREC_PRES))

corvif(tab_d.t %>% select(-starts_with('SR'),-Q_MAX, -Q_MIN, -Q_DOYCT,  -PREC_PRES, -TI))

corvif(tab_d.t %>% select(-starts_with('SR'),-Q_MAX, -Q_MIN, -Q_DOYCT,  -PREC_PRES, -DRAINAGE))

corvif(tab_d.t %>% select(-starts_with('SR'),-Q_MAX, -Q_MIN, -Q_DOYCT,  -PREC_PRES, -SLOPE))

#' * One model per flow variables.
#'  
#' * Excluding the central tendency of DOY decreased the VIF for both DOYMIN and DOYMAX.
#'  
#' * Excluding either drainage, TI or slope lowers the VIFs for those variables as well. I would kick out drainage as it makes us lose about 300 observations in the earlier stages of the modelling

#' ## RUN MODELS
#' 
#' 
#' ### Dredging the models

# modify lmer call for dredging
lmer.glmulti<-function(formula,data,random="",...) {
  newf <- formula
  newf[[3]] <- substitute(f+r,
                          list(f=newf[[3]],
                               r=reformulate(random)[[2]]))
  lme4::lmer(newf,data=data,
             REML=FALSE,...)
}

# add the BAS id to the transformed variables table
tab_d.t$BAS <- tab_d$BAS

# define variables
covariates_selection <- c(tab_d.t %>% select(-BAS,-starts_with('SR'), -starts_with('Q_M'), -Q_DOYCT, -DRAINAGE, -PREC_PRES) %>% colnames)
# covariates_selection <- tab_d.t %>% select(starts_with('Q'),-starts_with('Q_M'),TEMP_PRES,ELEVATION) %>% colnames

response_selection = 'SR_tot'
random_term <- 'BAS'
interaction_term <- c('TEMP_PRES','ELEVATION') # interactions with Q_magnitude variables
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')

res <- foreach(Qvar = Q_magnitude) %do% {
  
  df <- tab_d.t %>% select(all_of(response_selection),all_of(covariates_selection),all_of(Qvar),all_of(random_term))
  
  # df <- full_join(df,
  # df %>%
  #   group_by(BAS) %>%
  #   summarise(no_bas = n()))
  
  # GGally::ggpairs(df %>% select(-BAS),progress = FALSE) + theme_bw() + 
  #   theme(panel.grid =element_blank(),strip.background = element_blank())
  
  
  fit <- lmer(paste(
    response_selection,"~", # response
    paste(c(Qvar,covariates_selection),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*',Qvar,')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ), 
  data = df)
  
  
  # 
  # floc <- lm(SR_tot ~ Q_MEAN + Q_DOYMIN + ELEVATION + PREC_DELTA + Q_MEAN*ELEVATION, data = df %>% filter(BAS == 7120047060))
  # summary(floc)
  # plot(floc$model$SR_tot,predict(floc))
  # abline(0,1)
  
  # fit_lm <- glm(paste(
  #   
  #   response_selection,"~", # response
  #   
  #   paste(c(Qvar,covariates_selection),collapse=" + "), # fixed terms
  #    # interaction terms with Qvar
  #   collapse = ' + ')
  #   ,data = df)
  
  
  
  dred <- glmulti::glmulti(formula(fit,
                                   fixed.only=TRUE),
                           random=paste0("(1|",random_term,")"),
                           data=df, method="g",
                           conseq=3, crit=BIC,
                           fitfunc=lmer.glmulti,
                           includeobjects = TRUE,
                           confsetsize=1000, # The number of models to be looked for, i.e. the size of the returned confidence set.
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
  colnames(coef_tab) <- c('(Intercept)',Q_magnitude,covariates_selection,paste0('I(',interaction_term,' * ',Qvar,')'))
  
  
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
  write.csv(t_res,paste0('tabs/dredge_coefficients_',response_selection,'_',Qvar,'.csv'),row.names = F)
  
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
  
  write.csv(dfilt,paste0('tabs/dredge_coefficients_',response_selection,'_',names(res)[i],'_FILTERED.csv'),row.names = F)
  
  return(dfilt)
  
}






summary(dred)
plot(dred, type="r")
print(dred@objects[[1]])
MuMIn::r.squaredGLMM(dred@objects[[1]]) # R squared (marginal and conditional)

# to try later
library(snow)
cl<-makeCluster(8,type="SOCK")
clusterExport(cl, list("fit", "df", "lmer.glmulti","random_term","i"))

dred <- foreach(i = 1:8,.packages = 'glmulti') %dopar% {
  
  
  glmulti::glmulti(formula(fit,
                           fixed.only=TRUE),
                   random=paste0("(1|",random_term,")"),
                   data=df, method="h",
                   conseq=3, crit=BIC,
                   fitfunc=lmer.glmulti,
                   includeobjects = TRUE,
                   
                   confsetsize=100, # The number of models to be looked for, i.e. the size of the returned confidence set.
                   popsize = 100, # The population size for the genetic algorithm
                   # mutrate = 0.05,
                   chunk = i,
                   chunks = 8,
                   plotty = FALSE,
                   intercept=TRUE, marginality=FALSE, level=1)
  
}
Sys.time() - s
stopCluster(cl)



#### PRELIMINARY STUFF ----

## Some trials based on mean annual flow only

# 1a. Linear, log(SR) ~ log(MAF)
mod.1a <- glm(log(SR_tot) ~ log(Q_MEAN), data = tab, family = "gaussian") 
mod.1a.res <- resid(mod.1a)
plot(log(tab$Q_MEAN), mod.1a.res) 
abline(0,0) # Hmmm, doesn't look good

# 1b. Linear, mixed effect: log(SR) ~ log(MAF) + (1|basin)
mod.1b <- lmer(log(SR_tot) ~ log(Q_MEAN) + (1|BAS), data = tab) 
mod.1b.res <- resid(mod.1b)
plot(log(tab$Q_MEAN), mod.1b.res, ylim = c(-4,4)) 
abline(0,0) # 
qqnorm(mod.1b.res, xlim = c(-4,4), ylim = c(-4,4))
abline(0,1)
ks.test(mod.1b.res, "pnorm")

# 2a. Poisson, SR ~ log(MAF)
mod.2a <- glm(SR_tot ~ log(Q_MEAN), data = tab, family = "poisson") 
mod.2a.res <- resid(mod.2a)
plot(log(tab$Q_MEAN), mod.2a.res) 
abline(0,0) # Not good either

# 2b. Poisson, mixed effect: SR ~ log(MAF) + (1|basin)
mod.2b <- glmer(SR_tot ~ log(Q_MEAN) + (1|BAS), data = tab, family = "poisson")
mod.2b.res <- resid(mod.2b)
plot(log(tab$Q_MEAN), mod.2b.res) 
abline(0,0) # 
summary(mod.2b)

# Check possible overdisperion
hist(log(tab$SR_tot))
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(mod.2b) # Indeed overdispersed (assuming H0 represents no overdispersion)

# 3. Negative binomial to account for overdispersion, mixed effect: SR ~ log(MAF) + 1|basin
mod.3 <- glmer.nb(SR_tot ~ log(Q_MEAN) + (1|BAS), data = tab)
mod.3.res <- resid(mod.3)
plot(log(tab$Q_MEAN), mod.3.res, ylim = c(-4,4)) 
abline(0,0) # Heteroscedasticity gone, but we get error message (need to check)
summary(mod.3)
qqnorm(mod.3.res, xlim = c(-4,4), ylim = c(-4,4))
abline(0,1)
ks.test(mod.3.res, "pnorm")

## Conclusion: we go with a simple linear mixed effect model for now (with BAS as random effect).



