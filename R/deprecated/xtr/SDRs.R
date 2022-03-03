## SDRs - Valerio and Aafke
## January 2020

setwd("C:/Users/aafke/Documents/Projects/Valerio/SDR")

## Load packages
library(lme4) # for mixed effect modelling
library(ggplot2)
library(Hmisc)
library(glmulti) # for dredging models
library(dplyr)
#install.packages("MuMIn")
source("./HighstatLibV10.R") # For VIFs

## Read data
input <- read.csv("stations_filtered.csv", head = TRUE)
head(input)

## PREPARE DATA ----
input <- input[input$QUALITY !="Caution", ] # kick-out quality level 'caution' (meaning catchment area estimate is uncertain)
input$BAS <- as.factor(input$BAS)
input <- na.omit(input)

# Check distributions and log-tranform if skewed
hist.data.frame(input[,c("Q_MEAN", "Q_MIN", "Q_MAX", "Q_CV", "LATITUDE", "AREA", "TI", "ELEVATION", "DRAINAGE", "SR_tot", "SR_end")])
input[,c("Q_MEAN", "Q_MAX", "Q_CV", "AREA", "ELEVATION", "SR_tot")] <- log10(input[,c("Q_MEAN", "Q_MAX", "Q_CV", "AREA", "ELEVATION", "SR_tot")])
input$Q_MIN <- log10(input$Q_MIN + 1)
input$SR_tot <- log10(input$SR_tot + 1) # no zeros but for consistency
input$SR_end <- log10(input$SR_end + 1)

# Scale continuous variables because of clearly different units (gives issues in model fitting later on)
input[,c("Q_MEAN", "Q_MIN", "Q_MAX", "Q_CV", "LATITUDE", "AREA", "TI", "ELEVATION", "DRAINAGE", "SR_tot", "SR_end")] <- 
  scale(input[,c("Q_MEAN", "Q_MIN", "Q_MAX", "Q_CV", "LATITUDE", "AREA", "TI", "ELEVATION", "DRAINAGE", "SR_tot", "SR_end")])
data <- input
summary(data)

## Analysis of multicollinearity
# First bivariate correlations
cm <- cor(data[,c("Q_MEAN", "Q_MIN", "Q_MAX", "Q_CV", "LATITUDE", "TI", "ELEVATION", "DRAINAGE")],
          use = "pairwise.complete.obs")
# and visualize
corrplot::corrplot(cm, method = 'number', type = 'lower', number.cex = 0.8)
# Then VIFS
corvif(data[,c("Q_MEAN", "Q_MIN", "Q_MAX", "Q_CV", "CLIMATE", "LATITUDE", "TI", "ELEVATION", "DRAINAGE")])
## High VIFS (larger than 5) for the flow varaibles -> one model per flow metric


## MODELLING PART ----
# Full models (eight in total; one per flow metric and for all species and endemics separately)
full.mod.Q_MEAN.all <- lmer(paste(colnames(data)[14],"~", 
                       paste(colnames(data)[c(4,8:13)],collapse=" + "), 
                       "+ (1|BAS)"), data = data)
full.mod.Q_MEAN.end <- lmer(paste(colnames(data)[15],"~", 
                       paste(colnames(data)[c(4,8:13)],collapse=" + "), 
                        "+ (1|BAS)"), data = data)
full.mod.Q_MIN.all <- lmer(paste(colnames(data)[14],"~", 
                       paste(colnames(data)[c(5,8:13)],collapse=" + "), 
                       "+ (1|BAS)"), data = data)
full.mod.Q_MIN.end <- lmer(paste(colnames(data)[15],"~", 
                      paste(colnames(data)[c(5,8:13)],collapse=" + "), 
                      "+ (1|BAS)"), data = data)
full.mod.Q_MAX.all <- lmer(paste(colnames(data)[14],"~", 
                      paste(colnames(data)[c(6,8:13)],collapse=" + "), 
                      "+ (1|BAS)"), data = data)
full.mod.Q_MAX.end <- lmer(paste(colnames(data)[15],"~", 
                      paste(colnames(data)[c(6,8:13)],collapse=" + "), 
                      "+ (1|BAS)"), data = data)
full.mod.Q_CV.all <- lmer(paste(colnames(data)[14],"~", 
                      paste(colnames(data)[c(7,8:13)],collapse=" + "), 
                      "+ (1|BAS)"), data = data)
full.mod.Q_CV.end <- lmer(paste(colnames(data)[15],"~", 
                      paste(colnames(data)[c(7,8:13)],collapse=" + "), 
                      "+ (1|BAS)"), data = data)

# Dredging the models
lmer.glmulti<-function(formula,data,random="",...) {
  newf <- formula
  newf[[3]] <- substitute(f+r,
                          list(f=newf[[3]],
                               r=reformulate(random)[[2]]))
  lmer(newf,data=data,
       REML=FALSE,...)
}
random <- "(1|BAS)" 
all.mod.Q_MEAN.all <- glmulti(formula(full.mod.Q_MEAN.all,
                      fixed.only=TRUE),
                      random=random,
                      data=data, method="g",
                      conseq=3, crit=BIC,
                      fitfunc=lmer.glmulti,
                      includeobjects = TRUE,
                      confsetsize=100, # max number of models
                      intercept=TRUE, marginality=FALSE, level=1)
summary(all.mod.Q_MEAN.all)
plot(all.mod.Q_MEAN.all, type="r")
print(all.mod.Q_MEAN.all@objects[[2]])
MuMIn::r.squaredGLMM(all.mod.Q_MEAN.all@objects[[1]]) # R squared (marginal and conditional)


#### PRELIMINARY STUFF ----

## Some trials based on mean annual flow only

# 1a. Linear, log(SR) ~ log(MAF)
mod.1a <- glm(log(SR_tot) ~ log(Q_MEAN), data = input, family = "gaussian") 
mod.1a.res <- resid(mod.1a)
plot(log(input$Q_MEAN), mod.1a.res) 
abline(0,0) # Hmmm, doesn't look good

# 1b. Linear, mixed effect: log(SR) ~ log(MAF) + (1|basin)
mod.1b <- lmer(log(SR_tot) ~ log(Q_MEAN) + (1|BAS), data = input) 
mod.1b.res <- resid(mod.1b)
plot(log(input$Q_MEAN), mod.1b.res, ylim = c(-4,4)) 
abline(0,0) # 
qqnorm(mod.1b.res, xlim = c(-4,4), ylim = c(-4,4))
abline(0,1)
ks.test(mod.1b.res, "pnorm")

# 2a. Poisson, SR ~ log(MAF)
mod.2a <- glm(SR_tot ~ log(Q_MEAN), data = input, family = "poisson") 
mod.2a.res <- resid(mod.2a)
plot(log(input$Q_MEAN), mod.2a.res) 
abline(0,0) # Not good either

# 2b. Poisson, mixed effect: SR ~ log(MAF) + (1|basin)
mod.2b <- glmer(SR_tot ~ log(Q_MEAN) + (1|BAS), data = input, family = "poisson")
mod.2b.res <- resid(mod.2b)
plot(log(input$Q_MEAN), mod.2b.res) 
abline(0,0) # 
summary(mod.2b)

# Check possible overdisperion
hist(log(input$SR_tot))
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
mod.3 <- glmer.nb(SR_tot ~ log(Q_MEAN) + (1|BAS), data = input)
mod.3.res <- resid(mod.3)
plot(log(input$Q_MEAN), mod.3.res, ylim = c(-4,4)) 
abline(0,0) # Heteroscedasticity gone, but we get error message (need to check)
summary(mod.3)
qqnorm(mod.3.res, xlim = c(-4,4), ylim = c(-4,4))
abline(0,1)
ks.test(mod.3.res, "pnorm")

## Conclusion: we go with a simple linear mixed effect model for now (with BAS as random effect).


