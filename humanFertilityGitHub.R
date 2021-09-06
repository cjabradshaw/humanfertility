################################################################################################################
## R code to accompany paper:
## Perry, C, CJA Bradshaw, CM Saraswati, M Judge, J Heyworth, PN Le Souëf. In review. Lower infant mortality
## and access to contraception reduce fertility in low- and middle-income nations. Lancet Global Health
## September 2021
################################################################################################################

# libraries
library(lme4)
library(Hmisc)
library(ggplot2)
library(plotly)
library(nlme)
library(car)
library(dismo)
library(gbm)
library(rgeos)
library(rworldmap)
library(rworldxtra)
library(rcompanion)
library(SpatialEpi)
library(ggplot2)
library(ggridges)
library(dplyr)
library(ggpubr)
library(plyr)
library(fields)
library(ncf)
library(AICcmodavg)
library(modEvA)

# source files
source("new_lmer_AIC_tables3.R")
source("r.squared.R")

# set functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

## functions
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

# import data
dat <- read.csv("basedata.csv", header=T)
dim(dat)
head(dat)
str(dat)

dat.tfr.sort <- dat[order(dat[,5],decreasing=T), ]
dat.tfr.sort[, 1:5]

# tfr = total fertility rate
hist(dat$tfr,border="black",col="grey")
ggqqplot(dat, 'tfr')
shapiro.test(dat$tfr)


#####################
## PHASE 1: THEMES ##
#####################

#####################################################
# THEME 1: AVAILABILITY OF FAMILY PLANNING SERVICES
#####################################################
## availability of & access to family planning
hist(dat$access,border="black",col="grey")
ggqqplot(dat, 'access')
shapiro.test(dat$access)
access.sc <- scale(logit(dat$access/100), scale=T, center=T)
hist(access.sc,border="black",col="grey")
ggqqplot(access.sc)

# privsect = involvement of private-sector agencies/groups
hist(dat$privsect,border="black",col="grey")
ggqqplot(dat, 'privsect')
shapiro.test(dat$privsect)
privsect.sc <- scale(logit(dat$privsect/100), scale=T, center=T)
hist(privsect.sc,border="black",col="grey")
ggqqplot(access.sc)

# cbdistr = community-based distribution
hist(dat$cbdistr,border="black",col="grey")
ggqqplot(dat, 'cbdistr')
shapiro.test(dat$cbdistr)
cbdistr.sc <- scale(logit(dat$cbdistr/100), scale=T, center=T)
hist(cbdistr.sc,border="black",col="grey")
ggqqplot(cbdistr.sc)

# socmarkt = social marketing programs for subsidised contraceptive sales
hist(dat$socmarkt,border="black",col="grey")
ggqqplot(dat, 'socmarkt')
shapiro.test(dat$socmarkt)
socmarkt.sc <- scale(logit(dat$socmarkt/100), scale=T, center=T)
hist(socmarkt.sc,border="black",col="grey")
ggqqplot(socmarkt.sc)

# homevst = home-visiting workers for rural areas
hist(dat$homevst,border="black",col="grey")
ggqqplot(dat, 'homevst')
shapiro.test(dat$homevst)
homevst.sc <- scale(logit(dat$homevst/100), scale=T, center=T)
hist(homevst.sc,border="black",col="grey")
ggqqplot(homevst.sc)

# logtrans = logistics and transport for family planning services
hist(dat$logtrans,border="black",col="grey")
ggqqplot(dat, 'logtrans')
shapiro.test(dat$logtrans)
logtrans.sc <- scale(logit(dat$logtrans/100), scale=T, center=T)
hist(logtrans.sc,border="black",col="grey")
ggqqplot(logtrans.sc)
shapiro.test(logtrans.sc)

# correlation matrix
avail.dat <- data.frame(dat$tfr, access.sc,privsect.sc,cbdistr.sc,socmarkt.sc,homevst.sc,logtrans.sc)
colnames(avail.dat) <- c("fert","access","privsect","cbdist","socmrkt","homevis","logtrans")
cormat.avail <- cor(na.omit(avail.dat[,-c(1)]), method="kendall")
cormat.avail

# general linear model

# model set
m1 <- "fert ~ access + privsect + cbdist + socmrkt + homevis + logtrans"
m2 <- "fert ~ access + cbdist + socmrkt + homevis + logtrans"
m3 <- "fert ~ access + cbdist + socmrkt + logtrans"
m4 <- "fert ~ access + privsect"
m5 <- "fert ~ access + cbdist"
m6 <- "fert ~ access + socmrkt"
m7 <- "fert ~ access + homevis"
m8 <- "fert ~ access + logtrans"
m9 <- "fert ~ privsect + cbdist"
m10 <- "fert ~ privsect + socmrkt"
m11 <- "fert ~ privsect + homevis"
m12 <- "fert ~ privsect + logtrans"
m13 <- "fert ~ cbdist + socmrkt"
m14 <- "fert ~ cbdist + homevis"
m15 <- "fert ~ cbdist + logtrans"
m16 <- "fert ~ socmrkt + homevis"
m17 <- "fert ~ socmrkt + logtrans"
m18 <- "fert ~ homevis + logtrans"
m19 <- "fert ~ access"
m20 <- "fert ~ privsect"
m21 <- "fert ~ cbdist"
m22 <- "fert ~ socmrkt"
m23 <- "fert ~ homevis"
m24 <- "fert ~ logtrans"
m25 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=avail.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=avail.dat, na.action=na.omit)
plot(fit)

## boosted regression tree
brt.avail <- gbm.step(avail.dat, gbm.x = attr(avail.dat, "names")[c(2:7)], gbm.y = attr(avail.dat, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.avail)
gbm.plot(brt.avail)
tmp <- gbm.plot.fits(brt.avail, v=0)
D2 <- 100 * (brt.avail$cv.statistics$deviance.mean - brt.avail$self.statistics$mean.resid) / brt.avail$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.avail$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.avail$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

avail.dat.raw <- data.frame(dat$tfr, dat$access,dat$privsect,dat$cbdistr,dat$socmarkt,dat$homevst,dat$logtrans)
colnames(avail.dat.raw) <- c("fert","access","privsect","cbdist","socmrkt","homevis","logtrans")
brt.availraw <- gbm.step(avail.dat.raw, gbm.x = attr(avail.dat.raw, "names")[c(2:7)], gbm.y = attr(avail.dat.raw, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.availraw)
gbm.plot(brt.availraw)
tmp <- gbm.plot.fits(brt.availraw, v=0)
D2 <- 100 * (brt.availraw$cv.statistics$deviance.mean - brt.availraw$self.statistics$mean.resid) / brt.availraw$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.availraw$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.availraw$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

## ** keep homevis **



#####################################################
# THEME 2: QUALITY OF FAMILY PLANNING SERVICES
#####################################################

# quality = quality of family planning services 
hist(dat$quality,border="black",col="grey")
ggqqplot(dat, 'quality')
shapiro.test(dat$quality)
quality.sc <- scale(logit(dat$quality/100), scale=T, center=T)
hist(quality.sc,border="black",col="grey")
ggqqplot(quality.sc)

# anycontr = any contraceptive use
hist(dat$anycontr,border="black",col="grey")
ggqqplot(dat, 'anycontr')
shapiro.test(dat$anycontr)
anycontr.sc <- scale(logit(dat$anycontr/100), scale=T, center=T)
hist(anycontr.sc,border="black",col="grey")
ggqqplot(anycontr.sc)

# modcontr = modern contraceptive use
hist(dat$modcontr,border="black",col="grey")
ggqqplot(dat, 'modcontr')
shapiro.test(dat$modcontr)
modcontr.sc <- scale(logit(dat$modcontr/100), scale=T, center=T)
hist(modcontr.sc,border="black",col="grey")
ggqqplot(modcontr.sc)

# tradcontr = traditional contraceptive use
hist(dat$tradcontr,border="black",col="grey")
ggqqplot(dat, 'tradcontr')
shapiro.test(dat$tradcontr)
tradcontr.sc <- scale(logit(dat$tradcontr/100), scale=T, center=T)
hist(tradcontr.sc,border="black",col="grey")
ggqqplot(tradcontr.sc)

# folkcontr = folk method contraceptive use (REMOVE)
hist(dat$folkcontr,border="black",col="grey")
ggqqplot(dat, 'folkcontr')
shapiro.test(dat$folkcontr)
folkcontr.sc <- scale(logit(dat$folkcontr/100), scale=T, center=T)
hist(folkcontr.sc,border="black",col="grey")
ggqqplot(folkcontr.sc)

# nocontr = no contraceptive use (REMOVE)
hist(dat$nocontr,border="black",col="grey")
ggqqplot(dat, 'nocontr')
shapiro.test(dat$nocontr)
nocontr.sc <- scale(logit(dat$nocontr/100), scale=T, center=T)
hist(nocontr.sc,border="black",col="grey")
ggqqplot(nocontr.sc)

# correlation matrix
qual.dat <- data.frame(dat$tfr, quality.sc, anycontr.sc,modcontr.sc,tradcontr.sc,folkcontr.sc,nocontr.sc)
colnames(qual.dat) <- c("fert", "quality","anyContr","modContr","tradContr","folkContr","noContr")
cormat.qual <- cor(na.omit(qual.dat[,-c(1)]), method="kendall")
cormat.qual

# general linear model
# model set
m1 <- "fert ~ quality + anyContr + modContr + tradContr + noContr"
m2 <- "fert ~ quality + anyContr"
m3 <- "fert ~ quality + modContr"
m4 <- "fert ~ quality + tradContr"
m5 <- "fert ~ quality + noContr"
m6 <- "fert ~ anyContr + modContr"
m7 <- "fert ~ anyContr + tradContr"
m8 <- "fert ~ anyContr + noContr"
m9 <- "fert ~ modContr + tradContr"
m10 <- "fert ~ modContr + noContr"
m11 <- "fert ~ tradContr + noContr"
m12 <- "fert ~ quality"
m13 <- "fert ~ anyContr"
m14 <- "fert ~ modContr"
m15 <- "fert ~ tradContr"
m16 <- "fert ~ noContr"
m17 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=qual.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=qual.dat, na.action=na.omit)
plot(fit)

## boosted regression tree
# transformed values
brt.qual <- gbm.step(qual.dat, gbm.x = attr(qual.dat, "names")[c(2:5,7)], gbm.y = attr(qual.dat, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.qual)
gbm.plot(brt.qual)
tmp <- gbm.plot.fits(brt.qual, v=0)
D2 <- 100 * (brt.qual$cv.statistics$deviance.mean - brt.qual$self.statistics$mean.resid) / brt.qual$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.qual$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.qual$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

# raw values
qual.dat.raw <- data.frame(dat$tfr, dat$quality, dat$anycontr, dat$modcontr, dat$tradcontr, dat$nocontr)
colnames(qual.dat.raw) <- c("fert", "quality","anyContr","modContr","tradContr","noContr")
brt.qualraw <- gbm.step(qual.dat.raw, gbm.x = attr(qual.dat.raw, "names")[c(2:6)], gbm.y = attr(qual.dat.raw, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.qualraw)
gbm.plot(brt.qualraw)
tmp <- gbm.plot.fits(brt.qualraw, v=0)
D2 <- 100 * (brt.qualraw$cv.statistics$deviance.mean - brt.qualraw$self.statistics$mean.resid) / brt.qualraw$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.qualraw$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.qualraw$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

## ** keep anyContr **


#####################################################
# THEME 3: MATERNAL EDUCATION
#####################################################

## education
# primwom = female primary education completion
hist(dat$primwom,border="black",col="grey")
ggqqplot(dat, 'primwom')
shapiro.test(dat$primwom)
primwom.sc <- scale(logit(dat$primwom/100), scale=T, center=T)
hist(primwom.sc,border="black",col="grey")
ggqqplot(primwom.sc)

# secwom = female secondary education completion
hist(dat$secwom,border="black",col="grey")
ggqqplot(dat, 'secwom')
shapiro.test(dat$secwom)
secwom.sc <- scale(logit(dat$secwom/100), scale=T, center=T)
hist(secwom.sc,border="black",col="grey")
ggqqplot(secwom.sc)

# correlation matrix
ed.dat <- data.frame(dat$tfr, primwom.sc,secwom.sc)
colnames(ed.dat) <- c("fert","FprimEd","FsecEd")
cormat.ed <- cor(na.omit(ed.dat[,-c(1)]), method="kendall")
cormat.ed

# general linear model

# model set
m1 <- "fert ~ FprimEd + FsecEd + FprimEd*FsecEd"
m2 <- "fert ~ FprimEd + FsecEd"
m3 <- "fert ~ FprimEd"
m4 <- "fert ~ FsecEd"
m5 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=ed.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=ed.dat, na.action=na.omit)
plot(fit)

## ** keep FsecEd **


#####################################################
# THEME 4: RELIGION
#####################################################

# catholic_percentage = female primary education completion
hist(dat$catholic_percentage,border="black",col="grey")
ggqqplot(dat, 'catholic_percentage')
shapiro.test(dat$catholic_percentage)
cath <- scale(logit(dat$catholic_percentage/100), scale=T, center=T)
hist(cath,border="black",col="grey")
ggqqplot(cath)

## ** keep cath **


#####################################################
# THEME 5: MORTALITY
#####################################################

# imr = infant mortality rate
hist(dat$imr,border="black",col="grey")
ggqqplot(dat, 'imr')
shapiro.test(dat$imr)
imr.sc <- scale(logit(dat$imr/100), scale=T, center=T)
hist(imr.sc,border="black",col="grey")
ggqqplot(imr.sc)

# conflict-related death (crd_avg) (average 2010-2016)
hist(dat$crd_avg,border="black",col="grey")
ggqqplot(dat, 'crd_avg')

# import population data
popN <- read.table("pop.yr.csv", sep=",", header=T) # 2015
popN2015 <- popN[, c(1, dim(popN)[2])]
contreg <- read.table("continent.country2.csv", sep=",", header=T)

dat.mort1 <- data.frame(dat$country, dat$crd_avg)
colnames(dat.mort1) <- c("country","crd")
dim(dat.mort1)
int.cntry <- intersect(dat.mort1$country, contreg$country)
int.cntry
length(int.cntry)

# maternal mortality per 100,000 live births
matmort.dat <- read.csv("matmort.csv")

matmortmr <- rep(NA, dim(matmort.dat)[1])
for (i in 1:dim(matmort.dat)[1]) {
  maxsub <- max(which(is.na(matmort.dat[i,]) == F))
  matmortmr[i] <- ifelse(maxsub == 1, NA, matmort.dat[i,maxsub])
}
matmort <- data.frame(matmort.dat$cntry.code, matmortmr)
colnames(matmort) <- c("cntry.code","matmortmr")
head(matmort)
hist(log10(matmort$matmortmr))

dat.mort2 <- merge(dat.mort1, contreg, by="country")
dat.mort3 <- merge(dat.mort2, popN2015, by="cntry.code")
dat.mort4 <- merge(dat.mort3, matmort, by="cntry.code")

dat.mort <- dat.mort4[, c(2,3,6,7)]
dat.mort$crd1 <- ifelse(dat.mort$crd == 0, 1, dat.mort$crd)
dat.mort$matmortmr1 <- ifelse(dat.mort$matmortmr == 0, 1, dat.mort$matmortmr)
crdpc <- dat.mort$crd1/dat.mort$X2015 # per capita conflict-related deaths

hist(log10(crdpc),border="black",col="grey")
ggqqplot(log10(crdpc))
shapiro.test(na.omit(log10(crdpc)))
crdpc.sc <- scale(log10(crdpc), scale=T, center=T)
hist(crdpc.sc,border="black",col="grey")
ggqqplot(crdpc.sc)

matmort.sc <- scale(log10(dat.mort$matmortmr1), scale=T, center=T)
hist(matmort.sc,border="black",col="grey")
shapiro.test(na.omit(log10(crdpc)))
ggqqplot(matmort.sc)

# correlation matrix
mrt.dat <- data.frame(dat$tfr, imr.sc, matmort.sc, crdpc.sc)
colnames(mrt.dat) <- c("fert","infmort","matmort","crd")
cormat.mrt <- cor(na.omit(mrt.dat[,-c(1)]), method="kendall")
cormat.mrt

# general linear model

# model set
m1 <- "fert ~ infmort + matmort + crd + infmort*crd + infmort*matmort"
m2 <- "fert ~ infmort + matmort + crd + infmort*crd"
m3 <- "fert ~ infmort + matmort + crd + infmort*matmort"
m4 <- "fert ~ infmort + matmort + infmort*matmort"
m5 <- "fert ~ infmort + matmort + crd"
m6 <- "fert ~ infmort + crd + infmort*crd"
m7 <- "fert ~ infmort + crd"
m8 <- "fert ~ infmort + matmort"
m9 <- "fert ~  matmort + crd"
m10 <- "fert ~ infmort"
m11 <- "fert ~ matmort"
m12 <- "fert ~ crd"
m13 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=mrt.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=mrt.dat, na.action=na.omit)
plot(fit)
dim(na.omit(mrt.dat))

# transformed values
brt.mrt <- gbm.step(mrt.dat, gbm.x = attr(mrt.dat, "names")[2:4], gbm.y = attr(mrt.dat, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.mrt)
gbm.plot(brt.mrt)
tmp <- gbm.plot.fits(brt.mrt, v=0)
D2 <- 100 * (brt.mrt$cv.statistics$deviance.mean - brt.mrt$self.statistics$mean.resid) / brt.mrtl$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.mrt$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.mrt$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

## ** keep infmort **



#####################################################
# THEME 6: SOCIO-ECONOMICS
#####################################################

# qlowtotal = proportion of the population residing in low socioeconomic areas 
hist(dat$qlowtotal,border="black",col="grey")
ggqqplot(dat, 'qlowtotal')
shapiro.test(dat$qlowtotal)
qlowtotal.sc <- scale(logit(dat$qlowtotal/100), scale=T, center=T)
hist(qlowtotal.sc,border="black",col="grey")
ggqqplot(qlowtotal.sc)

# qhightotal = proportion of the population residing in high socioeconomic areas 
hist(dat$qhightotal,border="black",col="grey")
ggqqplot(dat, 'qhightotal')
shapiro.test(dat$qhightotal)
qhightotal.sc <- scale(logit(dat$qhightotal/100), scale=T, center=T)
hist(qhightotal.sc,border="black",col="grey")
ggqqplot(qhightotal.sc)

# threegentotal = three generations living in the same household
hist(dat$threegentotal,border="black",col="grey")
ggqqplot(dat, 'threegentotal')
shapiro.test(dat$threegentotal)
threegentotal.sc <- scale(logit(dat$threegentotal/100), scale=T, center=T)
hist(threegentotal.sc,border="black",col="grey")
ggqqplot(threegentotal.sc)

# meanmem = mean members of a household
hist(dat$meanmem,border="black",col="grey")
ggqqplot(dat, 'meanmem')
shapiro.test(dat$meanmem)
meanmem.sc <- scale(logit(dat$meanmem/100), scale=T, center=T)
hist(meanmem.sc,border="black",col="grey")
ggqqplot(meanmem.sc)

# correlation matrix
socioec.dat <- data.frame(dat$tfr, qlowtotal.sc,qhightotal.sc,threegentotal.sc,meanmem.sc)
colnames(socioec.dat) <- c("fert","SElo","SEhi","gen3","housesiz")
cormat.socioec <- cor(na.omit(socioec.dat[,-c(1)]), method="kendall")
cormat.socioec

# general linear model

# model set
m1 <- "fert ~ SElo + SEhi + gen3 + housesiz"
m2 <- "fert ~ SElo + SEhi"
m3 <- "fert ~ SElo + gen3"
m4 <- "fert ~ SElo + housesiz"
m5 <- "fert ~ SEhi + gen3"
m6 <- "fert ~ SEhi + housesiz"
m7 <- "fert ~ gen3 + housesiz"
m8 <- "fert ~ SElo"
m9 <- "fert ~ SEhi"
m10 <- "fert ~ housesiz"
m11 <- "fert ~ gen3"
m12 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=socioec.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=socioec.dat, na.action=na.omit)
plot(fit)

## boosted regression tree
# transformed values
brt.socioec <- gbm.step(socioec.dat, gbm.x = attr(socioec.dat, "names")[c(2:5)], gbm.y = attr(socioec.dat, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.socioec)
gbm.plot(brt.socioec)
tmp <- gbm.plot.fits(brt.socioec, v=0)
D2 <- 100 * (brt.socioec$cv.statistics$deviance.mean - brt.socioec$self.statistics$mean.resid) / brt.socioec$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.socioec$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.socioec$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

### ** keep gen3 **



###########################
## PHASE 2: GLOBAL MODEL ##
###########################

## homevis, FsecEd, cath, infmort, anyContr, gen3
final.dat <- data.frame(dat$tfr, homevst.sc, secwom.sc, cath, imr.sc, anycontr.sc, threegentotal.sc)
colnames(final.dat) <- c("fert", "homevis","FsecEd","cath","infmort","anyContr", "gen3")
cormat.final <- cor(na.omit(final.dat[,-c(1)]), method="kendall")
cormat.final

# general linear model

# model set
m1 <- "fert ~ homevis + FsecEd + cath + infmort + anyContr + gen3"
m2 <- "fert ~ homevis + FsecEd"
m3 <- "fert ~ homevis + cath"
m4 <- "fert ~ homevis + infmort"
m5 <- "fert ~ homevis + anyContr"
m6 <- "fert ~ homevis + gen3"
m7 <- "fert ~ FsecEd + cath"
m8 <- "fert ~ FsecEd + infmort"
m9 <- "fert ~ FsecEd + anyContr"
m10 <- "fert ~ FsecEd + gen3"
m11 <- "fert ~ cath + infmort"
m12 <- "fert ~ cath + anyContr"
m13 <- "fert ~ cath + gen3"
m14 <- "fert ~ infmort + anyContr"
m15 <- "fert ~ infmort + gen3"
m16 <- "fert ~ homevis"
m17 <- "fert ~ FsecEd"
m18 <- "fert ~ cath"
m19 <- "fert ~ infmort"
m20 <- "fert ~ anyContr"
m21 <- "fert ~ gen3"
m22 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=final.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=final.dat, na.action=na.omit)
plot(fit)

## boosted regression tree
# transformed values
brt.final <- gbm.step(final.dat, gbm.x = attr(final.dat, "names")[c(2:7)], gbm.y = attr(final.dat, "names")[1], family="gaussian", tolerance = 0.001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.final)
gbm.plot(brt.final)
tmp <- gbm.plot.fits(brt.final, v=0)
D2 <- 100 * (brt.final$cv.statistics$deviance.mean - brt.final$self.statistics$mean.resid) / brt.final$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.final$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.final$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))


## back-transform infant mortality to look at threshold effect
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
infmort.scale <- as.numeric(imr.sc) * attr(imr.sc, 'scaled:scale')
infmort.raw <- logit2prob(infmort.scale)
infmort.raw
fert.pred.md <- read.csv('BRT.boot.pred.med.FERT.csv')$infmort
fert.pred.up <- read.csv('BRT.boot.pred.up.FERT.csv')$infmort
fert.pred.lo <- read.csv('BRT.boot.pred.lo.FERT.csv')$infmort
infmort.val <- read.csv('BRT.boot.val.med.FERT.csv')$infmort
infmort.scale.fit <- as.numeric(infmort.val) * attr(imr.sc, 'scaled:scale')
infmort.raw.fit <- logit2prob(infmort.scale.fit)

plot(imr.sc, brt.final$fitted, pch=19)
plot(infmort.raw, brt.final$fitted, pch=19, xlab="infant mortality", ylab="fertility")
lines(infmort.raw.fit,fert.pred.md,lty=1)
lines(infmort.raw.fit,fert.pred.lo,lty=2,col="red")
lines(infmort.raw.fit,fert.pred.up,lty=2,col="red")

infmort.out <- data.frame(infmort.raw,brt.final$fitted)
colnames(infmort.out) <- c("infmort","fert.fit")
write.table(infmort.out,file="infmort.out.csv",sep=",", row.names = F, col.names = T)
write.table(infmort.raw.fit,file="infmortrawfit.csv",sep=",", row.names = F, col.names = T)

# raw values
final.dat.raw <- data.frame(dat$tfr, dat$homevst, dat$secwom, dat$catholic_percentage, dat$imr, dat$anycontr, dat$threegentotal)
colnames(final.dat.raw) <- c("fert", "homevis","FsecEd","cath","infmort","anyContr","gen3")
brt.finalraw <- gbm.step(final.dat.raw, gbm.x = attr(final.dat.raw, "names")[c(2:7)], gbm.y = attr(final.dat.raw, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.finalraw)
gbm.plot(brt.finalraw)
tmp <- gbm.plot.fits(brt.finalraw, v=0)
D2 <- 100 * (brt.finalraw$cv.statistics$deviance.mean - brt.finalraw$self.statistics$mean.resid) / brt.finalraw$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.finalraw$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.finalraw$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))


## generalised linear mixed-effects models

# first regroup regions
# ASIAPAC = East Asia/Pacific + South Asia (n = 12)
# EURCASMENA = Europe/Central Asia + Middle East/North Africa (n = 7)
# AFRICA = Sub-Saharan Africa (n = 34)
# LACAR = Latin America/Caribbean (n = 6)

dat$region2 <- ifelse((dat$region == "East Asia/Pacific" | dat$region == "South Asia"), "ASIAPAC", NA)
dat$region2 <- ifelse((dat$region == "Europe/Central Asia" | dat$region == "Middle East/North Africa"), "EURCASMENA", dat$region2)
dat$region2 <- ifelse((dat$region == "Sub-Saharan Africa"), "AFRICA", dat$region2)
dat$region2 <- as.factor(ifelse((dat$region == "Latin America/Caribbean"), "LACAR", dat$region2))
dat$region2
table(dat$region2)

final.dat2 <- data.frame(dat$country, dat$code, dat$region2, dat$tfr, homevst.sc, secwom.sc ,cath, imr.sc,anycontr.sc, threegentotal.sc)
colnames(final.dat2) <- c("country","cntry.code", "region2", "fert", "homevis","FsecEd","cath","infmort","anyContr", "gen3")
final.dat2

# model set
m1 <- "fert ~ homevis + FsecEd + cath + infmort + anyContr + gen3 + (1|region2)"
m2 <- "fert ~ homevis + FsecEd + (1|region2)"
m3 <- "fert ~ homevis + cath + (1|region2)"
m4 <- "fert ~ homevis + infmort + anyContr + (1|region2)"
m5 <- "fert ~ homevis + infmort + (1|region2)"
m6 <- "fert ~ homevis + anyContr + (1|region2)"
m7 <- "fert ~ homevis + gen3 + (1|region2)"
m8<- "fert ~ FsecEd + cath + (1|region2)"
m9 <- "fert ~ FsecEd + infmort + (1|region2)"
m10 <- "fert ~ FsecEd + anyContr + (1|region2)"
m11 <- "fert ~ FsecEd + gen3 + (1|region2)"
m12 <- "fert ~ cath + infmort + (1|region2)"
m13 <- "fert ~ cath + gen3 + (1|region2)"
m14 <- "fert ~ cath + anyContr + (1|region2)"
m15 <- "fert ~ infmort + anyContr + (1|region2)"
m16 <- "fert ~ infmort + gen3 + (1|region2)"
m17 <- "fert ~ homevis + (1|region2)"
m18 <- "fert ~ FsecEd + (1|region2)"
m19 <- "fert ~ cath + (1|region2)"
m20 <- "fert ~ infmort + (1|region2)"
m21 <- "fert ~ anyContr + (1|region2)"
m22 <- "fert ~ gen3 + (1|region2)"
m23 <- "fert ~ 1 + (1|region2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]),data=final.dat2, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,8],decreasing=F),]
summary.table


#######################################################################
## bootstrap iteration loop to estimate prediction confidence interval
iter <- 1000
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = 0, dim = c(eq.sp.points, 6, iter), dimnames=list(paste("x",1:eq.sp.points,sep=""), attr(final.dat, "names")[c(2:7)], paste("b",1:iter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- homevis.ri <- FsecEd.ri <- cath.ri <- infmort.ri <- anyContr.ri <- gen3.ri <- rep(0,iter)

# iterate
for (b in 1:iter) {  # start b loop
  
  # bootstrap data
  #boot.sub <- sort(sample(x = 1:dim(final.dat)[1], size = dim(final.dat)[1], replace=TRUE))
  boot.sub <- sort(sample(x = 1:dim(final.dat)[1], size = 50, replace=TRUE))
  final.boot <- final.dat[boot.sub,]
  
  brt.fit <- gbm.step(final.boot, gbm.x = attr(final.boot, "names")[c(2:7)], gbm.y = attr(final.boot, "names")[1], family="gaussian", tolerance = 0.001, learning.rate = 0.0001, bag.fraction=0.8, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000, silent=TRUE)
  summ.fit <- summary(brt.fit)
  
  homevis.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(2:7)][1])]
  FsecEd.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(2:7)][2])]
  cath.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(2:7)][3])]
  infmort.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(2:7)][4])]
  anyContr.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(2:7)][5])]
  gen3.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(2:7)][6])]
  
  #gbm.plot(brt.fit)
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
  D2.vec[b] <- D2
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.vec[b] <- CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se.vec[b] <- CV.cor.se

  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=6)
  
  ## output average predictions
  for (p in 1:6) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names
  
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  print(b)
}  # end b loop

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 6
pred.update <- pred.arr

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:iter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr, MARGIN=c(1,2), median)

par(mfrow=c(2,3)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←fewer) home-visiting workers for rural areas (more→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←lower) female secondary education completion (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←less) proportion Catholic (more→)")
lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←lower) infant mortality (higher→)" )
lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←less) access to any contraception (more→)" )
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

plot(val.med[,6],pred.med[,6],type="l",ylim=c(min(pred.lo[,6]),max(pred.up[,6])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←less) proportion households 3 gen (more→)" )
lines(val.med[,6], pred.lo[,6], type="l", lty=2, col="red")
lines(val.med[,6], pred.up[,6], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec
CV.cor.update <- CV.cor.vec
CV.cor.se.update <- CV.cor.se.vec
homevis.ri.update <- homevis.ri
FsecEd.ri.update <- FsecEd.ri
cath.ri.update <- cath.ri
infmort.ri.update <- infmort.ri
anyContr.ri.update <- anyContr.ri
gen3.ri.update <- gen3.ri

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  homevis.mean <- mean(homevis.ri.update, na.rm=T); homevis.sd <- sd(homevis.ri.update, na.rm=T)
  FsecEd.mean <- mean(FsecEd.ri.update, na.rm=T); FsecEd.sd <- sd(FsecEd.ri.update, na.rm=T)
  cath.mean <- mean(cath.ri.update, na.rm=T); cath.sd <- sd(cath.ri.update, na.rm=T)
  infmort.mean <- mean(infmort.ri.update, na.rm=T); infmort.sd <- sd(infmort.ri.update, na.rm=T)
  anyContr.mean <- mean(anyContr.ri.update, na.rm=T); anyContr.sd <- sd(anyContr.ri.update, na.rm=T)
  gen3.mean <- mean(gen3.ri.update, na.rm=T); gen3.sd <- sd(gen3.ri.update, na.rm=T)
  
  for (u in 1:iter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    homevis.ri.update[u] <- ifelse((homevis.ri.update[u] < (homevis.mean-kappa*homevis.sd) | homevis.ri.update[u] > (homevis.mean+kappa*homevis.sd)), NA, homevis.ri.update[u])
    FsecEd.ri.update[u] <- ifelse((FsecEd.ri.update[u] < (FsecEd.mean-kappa*FsecEd.sd) | FsecEd.ri.update[u] > (FsecEd.mean+kappa*FsecEd.sd)), NA, FsecEd.ri.update[u])
    cath.ri.update[u] <- ifelse((cath.ri.update[u] < (cath.mean-kappa*cath.sd) | cath.ri.update[u] > (cath.mean+kappa*cath.sd)), NA, cath.ri.update[u])
    infmort.ri.update[u] <- ifelse((infmort.ri.update[u] < (infmort.mean-kappa*infmort.sd) | infmort.ri.update[u] > (infmort.mean+kappa*infmort.sd)), NA, infmort.ri.update[u])
    anyContr.ri.update[u] <- ifelse((anyContr.ri.update[u] < (anyContr.mean-kappa*anyContr.sd) | anyContr.ri.update[u] > (anyContr.mean+kappa*anyContr.sd)), NA, anyContr.ri.update[u])
    gen3.ri.update[u] <- ifelse((gen3.ri.update[u] < (gen3.mean-kappa*anyContr.sd) | gen3.ri.update[u] > (gen3.mean+kappa*anyContr.sd)), NA, gen3.ri.update[u])
  }
  
  print(k)
}


D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

homevis.ri.lo <- quantile(homevis.ri.update, probs=0.025, na.rm=TRUE)
homevis.ri.med <- median(homevis.ri.update, na.rm=TRUE)
homevis.ri.up <- quantile(homevis.ri.update, probs=0.975, na.rm=TRUE)

FsecEd.ri.lo <- quantile(FsecEd.ri.update, probs=0.025, na.rm=TRUE)
FsecEd.ri.med <- median(FsecEd.ri.update, na.rm=TRUE)
FsecEd.ri.up <- quantile(FsecEd.ri.update, probs=0.975, na.rm=TRUE)

cath.ri.lo <- quantile(cath.ri.update, probs=0.025, na.rm=TRUE)
cath.ri.med <- median(cath.ri.update, na.rm=TRUE)
cath.ri.up <- quantile(cath.ri.update, probs=0.975, na.rm=TRUE)

infmort.ri.lo <- quantile(infmort.ri.update, probs=0.025, na.rm=TRUE)
infmort.ri.med <- median(infmort.ri.update, na.rm=TRUE)
infmort.ri.up <- quantile(infmort.ri.update, probs=0.975, na.rm=TRUE)

anyContr.ri.lo <- quantile(anyContr.ri.update, probs=0.025, na.rm=TRUE)
anyContr.ri.med <- median(anyContr.ri.update, na.rm=TRUE)
anyContr.ri.up <- quantile(anyContr.ri.update, probs=0.975, na.rm=TRUE)

gen3.ri.lo <- quantile(gen3.ri.update, probs=0.025, na.rm=TRUE)
gen3.ri.med <- median(gen3.ri.update, na.rm=TRUE)
gen3.ri.up <- quantile(gen3.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(homevis.ri.lo,FsecEd.ri.lo,cath.ri.lo,infmort.ri.lo,anyContr.ri.lo,gen3.ri.lo)
ri.med <- c(homevis.ri.med,FsecEd.ri.med,cath.ri.med,infmort.ri.med,anyContr.ri.med,gen3.ri.med)
ri.up <- c(homevis.ri.up,FsecEd.ri.up,cath.ri.up,infmort.ri.up,anyContr.ri.up,gen3.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.boot, "names")[c(2:7)]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort

write.table(pred.med,file="BRT.boot.pred.med.FERT.csv",sep=",", row.names = T, col.names = T) # median predicted values
write.table(pred.lo,file="BRT.boot.pred.lo.FERT.csv",sep=",", row.names = T, col.names = T) # lower 95% confidence bounds of predicted values
write.table(pred.up,file="BRT.boot.pred.up.FERT.csv",sep=",", row.names = T, col.names = T) # upper 95% confidence bounds of predicted values
write.table(val.med,file="BRT.boot.val.med.FERT.csv",sep=",", row.names = T, col.names = T) # range of values from which predictions derive


###############################
## general least-squares models

# get world map
wmap <- getMap(resolution="high")

# get centroids
centroids <- gCentroid(wmap, byid=TRUE)

# get a data.frame with centroids
centroids.df <- as.data.frame(centroids)
head(centroids.df)
centroids.df$country <- (row.names(centroids.df))
colnames(centroids.df) <- c("lon","lat","country")
head(centroids.df)

# standardise country names for merging
centroids.df$country <- revalue(centroids.df$country, c("Ivory Coast" = "Cote d'Ivoire",
                                              "Democratic Republic of the Congo" = "Congo Democratic Republic",
                                              "Republic of the Congo" = "Congo",
                                              "Kyrgyzstan" = "Kyrgyz Republic",
                                              "East Timor" = "Timor-Leste",
                                              "United Republic of Tanzania" = "Tanzania"))

# merge
final.dat.centroids <- merge(final.dat2, centroids.df, by="country")
dim(final.dat.centroids)

final.dat.c <- na.omit(final.dat.centroids)
final.dat.c$x <- latlong2grid(final.dat.c[,c(10:11)])$x # equidistant coordinates
final.dat.c$y <- latlong2grid(final.dat.c[,c(10:11)])$y # equidistant coordinates
plot(final.dat.c$x,final.dat.c$y,pch=19)

## determine best correlation structure
m1 <- gls(fert ~ homevis + FsecEd + cath + infmort + anyContr + gen3, data = final.dat.c)
vario1 <- Variogram(m1, form = ~ lon + lat, resType = "pearson")
plot(vario1, smooth = TRUE)

m2 <- gls(fert ~ homevis + FsecEd + cath + infmort + anyContr + gen3, correlation = corExp(form = ~ lon + lat, nugget = T), data = final.dat.c)
m3 <- gls(fert ~ homevis + FsecEd + cath + infmort + anyContr + gen3, correlation = corGaus(form = ~ lon + lat, nugget = T), data = final.dat.c)
m4 <- gls(fert ~ homevis + FsecEd + cath + infmort + anyContr + gen3, correlation = corSpher(form = ~ lon + lat, nugget = T), data = final.dat.c)
m5 <- gls(fert ~ homevis + FsecEd + cath + infmort + anyContr + gen3, correlation = corRatio(form = ~ lon + lat, nugget = T), data = final.dat.c)

mod.lab <- c("noCor","Exp","Gaus","Spher","Ratio")
AIC.vec <- c(AICc(m1), AICc(m2), AICc(m3), AICc(m4), AIC(m5))
#AIC.vec <- BIC(m1,m2,m3,m4,m5)$BIC
dAIC.vec <- delta.IC(AIC.vec)
wAIC.vec <- weight.IC(dAIC.vec)
psR2.mcf <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[1])
psR2.cs <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[2])
psR2.cu <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[3])
results.out <- data.frame(mod.lab,AIC.vec,dAIC.vec,wAIC.vec,psR2.mcf,psR2.cs,psR2.cu)
colnames(results.out) <- c("mod","AICc","dAICc","wAICc","psR2mcf","psR2cs","psR2cu")
results.sort <- results.out[order(results.out[,4],decreasing=T),1:7]
results.sort

# percentage of variance explained by geographic coordinates
100*(1 - results.sort[1,5]/(results.sort[4,5])) # pseudo R2 - McFadden
100*(1 - results.sort[1,6]/(results.sort[4,6])) # pseudo R2 - Cox & Snell
100*(1 - results.sort[1,7]/(results.sort[4,7])) # pseudo R2 - Craig & Uhler

vario4 <- Variogram(m2, form = ~ lon + lat, resType = "pearson")
plot(vario4, smooth = TRUE)
vario4.nr <- Variogram(m2, form = ~ lon + lat, resType = "normalized")
plot(vario4.nr, smooth = TRUE)

# run with Spherical spatial autocorrelation
# model set
m1 <- "fert ~ homevis + FsecEd + cath + infmort + anyContr + gen3"
m2 <- "fert ~ homevis + FsecEd"
m3 <- "fert ~ homevis + cath"
m4 <- "fert ~ homevis + infmort"
m5 <- "fert ~ homevis + anyContr"
m6 <- "fert ~ homevis + gen3"
m7 <- "fert ~ FsecEd + cath"
m8 <- "fert ~ FsecEd + infmort"
m9 <- "fert ~ FsecEd + anyContr"
m10 <- "fert ~ FsecEd + gen3"
m11 <- "fert ~ cath + infmort"
m12 <- "fert ~ cath + anyContr"
m13 <- "fert ~ cath + gen3"
m14 <- "fert ~ infmort + anyContr"
m15 <- "fert ~ infmort + gen3"
m16 <- "fert ~ homevis"
m17 <- "fert ~ FsecEd"
m18 <- "fert ~ cath"
m19 <- "fert ~ infmort"
m20 <- "fert ~ anyContr"
m21 <- "fert ~ gen3"
m22 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
SaveCount <- BIC.vec <- AICc.vec <- LL.vec <- k.vec <- psR2mcf <- psR2cs <- psR2cu <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- gls(as.formula(mod.vec[i]), correlation = corSpher(form = ~ lon + lat, nugget = T), method="ML", data = final.dat.c)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- logLik(fit)
  k.vec[i] <-  (AIC(fit) - -2*as.numeric(logLik(fit)))/2
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  psR2mcf[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[1] # McFadden pseudo-R2
  psR2cs[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[2] # Cox & Snell pseudo-R2
  psR2cu[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[3] # Craig & Uhler pseudo-R2
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

# re-interpret pseudo-R2s as % deviance from maximum pseudo-R2
psR2mcf.rel <- 100 - 100*(max(psR2mcf) - psR2mcf)/max(psR2mcf)
psR2cs.rel <- 100 - 100*(max(psR2cs) - psR2cs)/max(psR2cs)
psR2cu.rel <- 100 - 100*(max(psR2cu) - psR2cu)/max(psR2cu)

sumtable <- data.frame(mod.num,k.vec,round(LL.vec,3),round(AICc.vec,3),round(dAICc,3),round(wAICc,3),BIC.vec,round(dBIC,3),round(wBIC,3),round(psR2mcf,3),round(psR2mcf.rel,1),round(psR2cs,3),round(psR2cs.rel,1),round(psR2cu,3),round(psR2cu.rel,1))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","psR2mcf","mcfRel","psR2cs","csRel","psR2cu","cuRel")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,9],decreasing=T),1:15] # order by wBIC
summary.table

top.mod.sub <- summary.table[1,1]
summary(mod.list[[top.mod.sub]])
variotop <- Variogram(mod.list[[top.mod.sub]], form = ~ lon + lat, resType = "pearson")
plot(variotop, smooth = F, pch=19, grid=T)
variotop.nr <- Variogram(mod.list[[top.mod.sub]], form = ~ lon + lat, resType = "normalized")
plot(variotop.nr, smooth = F, pch=19, grid=T)

plot(mod.list[[top.mod.sub]],pch=19)
ggqqplot(mod.list[[top.mod.sub]]$residual)
plot(variotop.nr, smooth = F, pch=19, grid=T)

plot(mod.list[[1]],pch=19)
ggqqplot(mod.list[[1]]$residual)
variosat.nr <- Variogram(mod.list[[1]], form = ~ lon + lat, resType = "normalized")
plot(variosat.nr, smooth = F, pch=19, grid=T)

