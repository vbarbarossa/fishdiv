# set wd
knitr::opts_knit$set(root.dir = 'D:/projects/SDRs')

#' load assigned variables from MASTER.R
source('R/MASTER.R')

# what has been loaded
as.list(.GlobalEnv)

#' load packages needed
valerioUtils::libinv(c('dplyr','tidyr','ggplot2','foreach'))

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
  full_join(.,sf::read_sf('spatial/stations_catchments2.gpkg') %>% 
              as_tibble() %>% dplyr::select(-geom,-area.est)) %>%
  full_join(.,read.csv('tabs/stations_SR.csv'),by='gsim.no') %>%
  as_tibble()

write.csv(tab_all,'tabs/stations_all_attributes.csv',row.names = F)

# compute the correlation matrix 
cm <- cor(tab_all %>% select(colnames(tab_all)[92:which(colnames(tab_all) == 'DOYMAX7')]),
          use = "pairwise.complete.obs")

# and visualize
corrplot::corrplot(cm, method = 'number', type = 'lower',number.cex = 0.8)

# and save
jpeg('figs/corrplot_initial_flow_indices.jpg',width = 200, height = 200, res = 600, units = 'mm')
corrplot::corrplot(cm, method = 'number', type = 'lower',number.cex = 0.8)
dev.off()

# check dams GOOD2 vs GRanD
plot(log10(tab_all$dams_good2_no),log10(tab_all$no.dams))
abline(0,1) 

# interesting, some catchments have more grand dams than good2
# true for
sum((tab_all$dams_good2_no-tab_all$no.dams) < 0) # catchments

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

tab <- tab_all %>%
  # variables that need to be computed
  mutate(
    # Q_var = IQR/P50
    PREC_DELTA = abs(prec_cur_mean - prec_hist_mean),
    TEMP_DELTA = abs(temp_cur_mean - temp_hist_mean),
    LATITUDE = abs(LAT),
    CROP_PAST = cropland_1992_sum/cropland_1992_count,
    CROP_PRES = cropland_2015_sum/cropland_2015_count
  ) %>%
  mutate(
    CROP_DELTA = CROP_PRES - CROP_PAST
  ) %>%
  # select the variables and rename them
  select(
    # ID variables
    ID = gsim.no, BAS = MAIN_BAS, QUALITY = quality,
    # Discharge covariates
    Q_MEAN = MEAN, Q_MIN = MIN7, Q_MAX = MAX7, Q_CV = CV, Q_DOYMIN = DOYMIN7, Q_DOYMAX = DOYMAX7, Q_DOYCT = CT,
    # Ecosystem productivity
    CLIMATE = climate.type,
    PREC_PRES = prec_cur_mean, PREC_HIST = prec_hist_mean, PREC_DELTA,
    TEMP_PRES = temp_cur_mean, TEMP_HIST = temp_hist_mean, TEMP_DELTA,
    # Evolutionary diversification potential
    AREA = area.est, TI = tp.mean, ELEVATION = ele.mean, DRAINAGE = dr.mean, SLOPE = slp.mean,
    # Anthropogenic
    POP = pop.count, POP_DENS = pd.mean, NO_DAMS = dams_good2_no, URB = nl.mean,
    CROP_PAST,CROP_PRES,CROP_DELTA,
    # Response
    SR_tot, SR_end
  )

#' kick-out quality level 'caution' (meaning catchment area estimate is uncertain)
tab <- tab %>%
  filter(QUALITY !="Caution")

#' check for NAs and flow/precipitation metrics < 0 and remove those records

# NAs
apply(tab,2,function(x) sum(is.na(x)))

# looks like URB has NAs instead of zeroes
tab$URB[is.na(tab$URB)] <- 0

# Q < 0
apply(tab %>% select(starts_with('Q'),-QUALITY,starts_with('PREC')),2,function(x) sum(x < 0,na.rm=T))

# remove

tab <- tab %>%
  drop_na() %>% # we are also dropping 351 basins with NAs in DRAINAGE
  filter(Q_MIN >= 0)
tab


#' make main basin ID as factor
tab$BAS <- as.factor(tab$BAS)

#'  save the final table
write.csv(tab,'tabs/stations_filtered.csv',row.names = F)


#'***********************************************************************************************
#' # Modelling
#' 
#' ## Some diagnostics first
#' 
#' ### Check endemic vs non-endemic respose variables
p <- ggplot(tab) +
  geom_point(aes(x = log10(SR_tot+1), y = log10(SR_end+1))) +
  xlab('Log-SR') +
  ylab('Log-SR Endemic Only') +
  theme_bw()
p
ggsave('figs/bivariate_SRtot_vs_SRendemic.jpg',dpi = 600)

#' Load packages
library(lme4) # for mixed effect modelling
library(Hmisc)
library(glmulti) # for dredging models
source("R/HighstatLibV10.R") # For VIFs

#' ### Variables distribution 
#' 
#' Select only strictly needed variables
tab <- tab %>%
  select(-QUALITY,-PREC_HIST,-TEMP_HIST)

#' Check distributions and log-tranform if skewed
valerioUtils::libinv(c('ggplot2','purrr','tidyr'))
tab %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

tab[,c("Q_MEAN", "Q_MAX", "AREA", "ELEVATION","SLOPE")] <- log10(tab[,c("Q_MEAN", "Q_MAX", "AREA", "ELEVATION","SLOPE")])
tab$Q_MIN <- log10(tab$Q_MIN + 0.001) # approx the minimum value of Q mean is ~ 10^-3
tab$SR_tot <- log10(tab$SR_tot + 1) # no zeros but for consistency
tab$SR_end <- log10(tab$SR_end + 1) # zero values in here

# Scale continuous variables because of clearly different units (gives issues in model fitting later on)
tab[,colnames(tab)[-which(colnames(tab) %in% c('ID','BAS','CLIMATE'))]] <- 
  scale(tab[,colnames(tab)[-which(colnames(tab) %in% c('ID','BAS','CLIMATE'))]])
summary(tab)

# and check again the distribution
tab %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

#' ## Multicollinearity
#' 
#' First bivariate correlations
cm <- cor(tab %>% select(-ID,-BAS,-CLIMATE,-starts_with('SR')),
          use = "pairwise.complete.obs")
# and visualize
corrplot::corrplot(cm, method = 'number', type = 'lower', number.cex = 0.8)
# Then VIFS
corvif(tab %>% select(-ID,-BAS,-starts_with('SR')))

#' **High VIFS (larger than 5) for the flow varaibles -> one model per flow metric**

#' ## RUN MODELS
#' 

#' ## Dredging the models

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
covariates_selection <- c('TEMP_PRES','TEMP_DELTA','PREC_PRES','PREC_DELTA','AREA','TI','ELEVATION','Q_CV','Q_DOYMIN','Q_DOYMAX')
# covariates_selection <- c('TEMP_PRES','AREA','ELEVATION','SLOPE')

response_selection = 'SR_tot'
random_term <- 'BAS'
interaction_term <- c('TEMP_PRES','ELEVATION')
Q_magnitude <- c('Q_MEAN','Q_MIN','Q_MAX')

res <- foreach(Qvar = Q_magnitude) %do% {
  
  df <- tab %>% select(response_selection,covariates_selection,Qvar,random_term)
  
  fit <- lmer(paste(
    response_selection,"~", # response
    paste(c(Qvar,covariates_selection),collapse=" + "), # fixed terms
    '+',
    paste0('I(',interaction_term,'*',Qvar,')',collapse = ' + '), # interaction terms with Qvar
    '+',
    paste0("(1|",random_term,")") # random term
  ), 
  data = df)
  
  dred <- glmulti::glmulti(formula(fit,
                                   fixed.only=TRUE),
                           random=paste0("(1|",random_term,")"),
                           data=df, method="h",
                           conseq=3, crit=BIC,
                           fitfunc=lmer.glmulti,
                           includeobjects = TRUE,
                           
                           confsetsize=100, # The number of models to be looked for, i.e. the size of the returned confidence set.
                           popsize = 100, # The population size for the genetic algorithm
                           # chunk = i,
                           # chunks = 8,
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



