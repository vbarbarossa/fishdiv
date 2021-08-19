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
tab_d <- read.csv('tabs/stations_filtered.csv') %>%
  as_tibble() %>%
  select(ID,BAS,SR_tot, Q_MEAN)

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

#' ## test mixed effect models
#' 

# add the BAS id to the transformed variables table
tab_d.t$BAS <- tab_d$BAS

df <- tab_d.t

# add number of observations per main basin
df <- full_join(df,
                df %>%
                  group_by(BAS) %>%
                  summarise(no_bas = n()))

#' # Case 1: entire dataset
d1 <- df

# number of data points
nrow(d1)

# number of basins
length(unique(d1$BAS))

# number of basins with only one observation?
length(unique(d1$BAS[d1$no_bas == 1]))

#' Overall regression
fit <- lm(SR_tot ~ Q_MEAN, data = d1)
summary(fit)

d1$BAS <- as.factor(d1$BAS)
(mm_plot <- ggplot(d1, aes(x = Q_MEAN, y = SR_tot, colour = BAS)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(data = cbind(d1, pred = predict(fit)), aes(y = pred), size = 1, alpha = 0.3) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

#' ### HP1: different number of sampling points per waterhsed influences the SDR validity
#' 
#' Run lmer with only random intercept
fit1 <- lmer(SR_tot ~ Q_MEAN + (1 | BAS), REML = F,
             data = d1)
summary(fit1)
#' looks like BAS explains a lot of variance 0.58/(0.37+0.58) * 100 = ~61%
#' 
#' check the R2
# marginal and conditional
MuMIn::r.squaredGLMM(fit1)
# common way
valerioUtils::r.squared(d1$SR_tot,predict(fit1))

#' visualize
d1$BAS <- as.factor(d1$BAS)
(mm_plot <- ggplot(d1, aes(x = Q_MEAN, y = SR_tot, colour = BAS)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(data = cbind(d1, pred = predict(fit1)), aes(y = pred), size = 1, alpha = 0.3) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

#' **Q1: are the intercept values assigned to the basins with one observation different?**
#' 
# intercept values for basins with only one observation
boxplot(
  coef(fit1)[[1]][row.names(coef(fit1)[[1]]) %in% unique(d1$BAS[d1$no_bas == 1]),'(Intercept)'],
  ylab = 'Intercept value'
)
#' **R1: Yes, one different intercept value for each basin**
#' 
#' 
#' ### HP2: different number of sampling points per watershed influences the SDR and different watersheds have different increase rates of SDRs
#' 
#' Run lmer with random intercep and slope
fit1 <- lmer(SR_tot ~ Q_MEAN + (1 + Q_MEAN | BAS), REML = F,
             data = d1)
summary(fit1)
#' random intercept explains 0.55/(0.30 + 0.09 + 0.55) * 100 = ~58% of variance
#' 
#' random slope explains 0.09/(0.30 + 0.09 + 0.55) * 100 = ~10% of variance 
#' 
#' check the R2
# marginal and conditional
MuMIn::r.squaredGLMM(fit1)
# common way
valerioUtils::r.squared(d1$SR_tot,predict(fit1))

#' visualize
d1$BAS <- as.factor(d1$BAS)
(mm_plot <- ggplot(d1, aes(x = Q_MEAN, y = SR_tot, colour = BAS)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(data = cbind(d1, pred = predict(fit1)), aes(y = pred), size = 1, alpha = 0.3) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

#' the mixed effect model uses one intercept and slope per basin. 
#' My doubt is whether this might be inflating the explained variance 
#' of the random terms because we have many basins with only one observation
#' 
#' 
#' # Case 2: exclude basins with less than 20 observations
#' 
#' 
#' Let's repeat what has been done but using basins with only 20 or more sub-basins

d1 <- df %>% filter(no_bas >= 20)

# number of data points
nrow(d1)

# number of basins
length(unique(d1$BAS))

#' Overall regression
fit <- lm(SR_tot ~ Q_MEAN, data = d1)
summary(fit)

d1$BAS <- as.factor(d1$BAS)
(mm_plot <- ggplot(d1, aes(x = Q_MEAN, y = SR_tot, colour = BAS)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(data = cbind(d1, pred = predict(fit)), aes(y = pred), size = 1, alpha = 0.3) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

#' ### HP1: different number of sampling points per waterhsed influences the SDR validity
#' 
#' Run lmer with only random intercept
fit1 <- lmer(SR_tot ~ Q_MEAN + (1 | BAS), REML = F,
             data = d1)
summary(fit1)
#' random intercept explains 0.47/(0.47+0.43) * 100 = ~52% of variance
#' 
#' check the R2
# marginal and conditional
MuMIn::r.squaredGLMM(fit1)
# common way
valerioUtils::r.squared(d1$SR_tot,predict(fit1))

#' visualize
d1$BAS <- as.factor(d1$BAS)
(mm_plot <- ggplot(d1, aes(x = Q_MEAN, y = SR_tot, colour = BAS)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(data = cbind(d1, pred = predict(fit1)), aes(y = pred), size = 1, alpha = 0.3) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

#' ### HP2: different number of sampling points per watershed influences the SDR and different watersheds have different increase rates of SDRs
#' 
#' Run lmer with random intercep and slope
fit1 <- lmer(SR_tot ~ Q_MEAN + (1 + Q_MEAN | BAS), REML = F,
             data = d1)
summary(fit1)
#' random intercept explains 0.55/(0.55 + 0.10 + 0.36) * 100 = ~54% of variance
#' 
#' random slope explains 0.10/(0.55 + 0.10 + 0.36) * 100 = ~10% of variance 
#' 
#' check the R2
# marginal and conditional
MuMIn::r.squaredGLMM(fit1)
# common way
valerioUtils::r.squared(d1$SR_tot,predict(fit1))

#' visualize
d1$BAS <- as.factor(d1$BAS)
(mm_plot <- ggplot(d1, aes(x = Q_MEAN, y = SR_tot, colour = BAS)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(data = cbind(d1, pred = predict(fit1)), aes(y = pred), size = 1, alpha = 0.3) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

#' multipanel
#+ fig.width=10, fig.height=16
ggplot(d1, aes(x = Q_MEAN, y = SR_tot, colour = BAS)) +
  geom_point(alpha = 0.5) +
  facet_wrap('BAS',ncol = 5) +
  theme_bw() +
  geom_line(data = cbind(d1, pred = predict(fit1)), aes(y = pred), size = 1, alpha = 0.3) +  # adding predicted line from mixed model
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"))  # adding space between panels


#' **Q2: is this last model (random slope and intercept) equivalent to fitting one lm per main basin?**

# one lm for each main basin
(multi_lm <- foreach(id = row.names(coef(fit1)[[1]]),.combine = 'rbind') %do% {
  
  dt <- d1 %>% filter(BAS == id)
  res <- coef(lm(SR_tot ~ Q_MEAN, data = dt)) %>%
    as.data.frame() %>% t() %>% as.data.frame()# %>% select('Q_MEAN','(Intercept)')
  row.names(res) <- id
  
  return(res)
  
})

# results from the lmer with random slope and intercept
(single_lmer <- coef(fit1)[[1]])

# % differences
round((multi_lm - single_lmer)/single_lmer * 100,1)

#' **R2: it does not seem to be the same thing, somehow the lmer still keeps into account for the other data points when calculating slope and intercept, nice!**
#' 