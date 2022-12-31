################################################################################################################
## R code to accompany paper:
##
## BRADSHAW, CJA, C PERRY, M JUDGE, CM SARASWATI, J HEYWORTH, PN LE SOUËF. 2023. Lower infant mortality,
## household size, and access to contraception reduce fertility in low- and middle-income nations. PLoS One. In press
##
## Corey J. A. Bradshaw
## Flinders University
## August 2022 / updated December 2022
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
library(ggstatsplot)

# source files
source("new_lmer_AIC_tables3.R") # change path as required
source("r.squared.R") # change path as required

# Set functions
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
dat <- read.csv("basedata.update.csv", header=T) # updated data

dat.tfr.sort <- dat[order(dat[,6],decreasing=T), ]
dat.tfr.sort[, 1:6]


## response(s)
# tfr = total fertility rate
hist(dat$tfr10t20,border="black",col="grey")
ggqqplot(dat, 'tfr10t20')
shapiro.test(dat$tfr10t20)

## availability of & access to family planning
# access = access
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
avail.dat <- data.frame(dat$tfr10t20, access.sc,privsect.sc,cbdistr.sc,socmarkt.sc,homevst.sc,logtrans.sc)
colnames(avail.dat) <- c("fert","access","privsect","cbdist","socmrkt","homevis","logtrans")
cormat.avail <- cor(na.omit(avail.dat[,-c(1)]), method="kendall")
cormat.avail

dim(na.omit(avail.dat))

# GLM
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


## ** keep 'homevis' **


# weighted years of education (females)
hist(dat$educyrsmn,border="black",col="grey")
ggqqplot(dat, 'educyrsmn')
shapiro.test(dat$educyrsmn)
educyrsmn.sc <- scale((dat$educyrsmn), scale=T, center=T)
hist(educyrsmn.sc,border="black",col="grey")
ggqqplot(educyrsmn.sc)

# correlation matrix
ed.dat <- data.frame(dat$tfr10t20, educyrsmn.sc)
colnames(ed.dat) <- c("fert","educyrsmn")
cormat.ed <- cor(na.omit(ed.dat, method="kendall"))
cormat.ed

# GLM
# model set
m1 <- "fert ~ educyrsmn"
m2 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2)

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


## religion (Catholics + Muslims)
hist(dat$ardacathmusl2020,border="black",col="grey")
ggqqplot(dat, 'ardacathmusl2020')
shapiro.test(dat$ardacathmusl2020)
cathmusl.sc <- scale(logit(dat$ardacathmusl2020/100), scale=T, center=T)
hist(cathmusl.sc,border="black",col="grey")
ggqqplot(cathmusl.sc)

## ** keep 'cathmusl' **


# mortality
# imr = infant mortality rate
hist(dat$imr2020,border="black",col="grey")
ggqqplot(dat, 'imr2020')
shapiro.test(dat$imr2020)
imr.sc <- scale(logit(dat$imr2020/100), scale=T, center=T)
hist(imr.sc,border="black",col="grey")
ggqqplot(imr.sc)


# import population data
popN <- read.table("pop.yr.update.csv", sep=",", header=T) # 2021
popN2021 <- popN[, c(1, dim(popN)[2])]
colnames(popN2021)[2] <- c("pop2021")
contreg <- read.table("continent.country2.csv", sep=",", header=T)


# battle-related deaths
brd <- read.table("battledeaths.csv", sep=",", header=T)
brdcols <- dim(brd)[2]
brd$brdAvg5yr <- rowSums(brd[,(brdcols-5):brdcols], na.rm=T)
brdpop <- merge(brd, popN2021, "cntry.code")
brdpop$brdAvg5yrPC <- brdpop$brdAvg5yr/brdpop$pop2021
brdpop$lbrdAvg5yrPC <- ifelse(is.infinite(log10(brdpop$brdAvg5yrPC)) == T, NA, log10(brdpop$brdAvg5yrPC))

brdpopdat <- merge(dat, brdpop, "cntry.code")
hist((brdpopdat$lbrdAvg5yrPC),border="black",col="grey")
ggqqplot((brdpopdat$lbrdAvg5yrPC))
shapiro.test(na.omit((brdpopdat$lbrdAvg5yrPC)))
brd.sc <- scale((brdpopdat$lbrdAvg5yrPC), scale=T, center=T)
hist(brd.sc,border="black",col="grey")
ggqqplot(brd.sc)

# maternal mortality per 100,000 live births
matmort.dat <- read.csv("matmort.update.csv")

matmortmr <- rep(NA, dim(matmort.dat)[1])
for (i in 1:dim(matmort.dat)[1]) {
  maxsub <- max(which(is.na(matmort.dat[i,]) == F))
  matmortmr[i] <- ifelse(maxsub == 1, NA, matmort.dat[i,maxsub])
}
matmort <- data.frame(matmort.dat$cntry.code, matmortmr)
colnames(matmort) <- c("cntry.code","matmortmr")
head(matmort)
hist(log10(matmort$matmortmr))

matmortmrdat <- merge(dat, matmort, "cntry.code")
matmort.sc <- scale(log10(matmortmrdat$matmortmr), scale=T, center=T)
hist(matmort.sc,border="black",col="grey")
ggqqplot(matmort.sc)
shapiro.test(na.omit(matmort.sc))

# correlation matrix
mrt.dat <- data.frame(dat$tfr10t20, imr.sc, matmort.sc, brd.sc)
colnames(mrt.dat) <- c("fert","infmort","matmort","brd")
cormat.mrt <- cor(na.omit(mrt.dat[,-c(1)]), method="kendall")
cormat.mrt

# GLM
# model set
m1 <- "fert ~ infmort + matmort + brd + infmort*brd + infmort*matmort"
m2 <- "fert ~ infmort + matmort + brd + infmort*brd"
m3 <- "fert ~ infmort + matmort + brd + infmort*matmort"
m4 <- "fert ~ infmort + matmort + infmort*matmort"
m5 <- "fert ~ infmort + matmort + brd"
m6 <- "fert ~ infmort + brd + infmort*brd"
m7 <- "fert ~ infmort + brd"
m8 <- "fert ~ infmort + matmort"
m9 <- "fert ~  matmort + brd"
m10 <- "fert ~ infmort"
m11 <- "fert ~ matmort"
m12 <- "fert ~ brd"
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
D2 <- 100 * (brt.mrt$cv.statistics$deviance.mean - brt.mrt$self.statistics$mean.resid) / brt.mrt$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.mrt$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.mrt$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

## ** keep 'infmort' **


## socio-economics
# npwb50 = net personal wealth bottom 50% (neg values = debt)
hist(dat$npwb50,border="black",col="grey")
ggqqplot(dat, 'npwb50')
shapiro.test(dat$npwb50)
npwb50.sc <- scale((dat$npwb50)^3, scale=T, center=T)
hist(npwb50.sc,border="black",col="grey")
ggqqplot(npwb50.sc)
shapiro.test(npwb50.sc)

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
socioec.dat <- data.frame(dat$tfr10t20, npwb50.sc, threegentotal.sc,meanmem.sc)
colnames(socioec.dat) <- c("fert","inequal","gen3","housesiz")
cormat.socioec <- cor(na.omit(socioec.dat[,-c(1)]), method="kendall")
cormat.socioec

# GLM
# model set
m1 <- "fert ~ inequal + gen3 + housesiz"
m2 <- "fert ~ inequal"
m3 <- "fert ~ inequal + gen3"
m4 <- "fert ~ inequal + housesiz"
m5 <- "fert ~ gen3 + housesiz"
m6 <- "fert ~ housesiz"
m7 <- "fert ~ gen3"
m8 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8)

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
brt.socioec <- gbm.step(socioec.dat, gbm.x = attr(socioec.dat, "names")[c(2:4)], gbm.y = attr(socioec.dat, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.socioec)
gbm.plot(brt.socioec)
tmp <- gbm.plot.fits(brt.socioec, v=0)
D2 <- 100 * (brt.socioec$cv.statistics$deviance.mean - brt.socioec$self.statistics$mean.resid) / brt.socioec$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.socioec$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.socioec$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

### keep housesize


## quality
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

# nocontr = no contraceptive use (REMOVE)
hist(dat$nocontr,border="black",col="grey")
ggqqplot(dat, 'nocontr')
shapiro.test(dat$nocontr)
nocontr.sc <- scale(logit(dat$nocontr/100), scale=T, center=T)
hist(nocontr.sc,border="black",col="grey")
ggqqplot(nocontr.sc)

# correlation matrix
qual.dat <- data.frame(dat$tfr, quality.sc, anycontr.sc,modcontr.sc,tradcontr.sc,nocontr.sc)
colnames(qual.dat) <- c("fert", "quality","anyContr","modContr","tradContr","noContr")
cormat.qual <- cor(na.omit(qual.dat[,-c(1)]), method="kendall")
cormat.qual


# GLM
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
dim(na.omit(qual.dat))

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=qual.dat, na.action=na.omit)
plot(fit)

## boosted regression tree
# transformed values
brt.qual <- gbm.step(qual.dat, gbm.x = attr(qual.dat, "names")[c(2:6)], gbm.y = attr(qual.dat, "names")[1], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
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

## ** keep 'anyContr' + 'quality'**


## multiple imputation using mice
final.dat.raw <- data.frame(dat$homevst, dat$educyrsmn, dat$ardacathmusl2020, dat$imr2020, dat$anycontr, dat$quality, dat$meanmem)
colnames(final.dat.raw) <- c("homevis","Feduc","CathMusl","infmort","anyContr","quality", "housesiz")
dim(final.dat.raw)
dim(na.omit(final.dat.raw))

library(VIM)
aggr_plot <- aggr(final.dat.raw, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(final.dat.raw), cex.axis=.7, gap=3, ylab=c("histogram of missing data","pattern"))
marginplot(final.dat.raw[c(1,2)])

library(mice)
md.pattern(final.dat.raw)
final.dat.imp <- mice(final.dat.raw, m=7, maxit=50,method="pmm", seed=101)
summary(final.dat.imp)

final.dat.imp$imp$homevis
final.dat.imp$imp$Feduc
final.dat.imp$imp$CathMusl
final.dat.imp$imp$infmort
final.dat.imp$imp$anyContr
final.dat.imp$imp$quality
final.dat.imp$imp$housesiz


xyplot(final.dat.imp, homevis ~ Feduc+CathMusl+infmort+anyContr+quality+housesiz,pch=18,cex=1)
densityplot(final.dat.imp)
stripplot(final.dat.imp, pch = 20, cex = 1.2)


final.dat.compl1 <- complete(final.dat.imp, 1)
final.dat.compl2 <- complete(final.dat.imp, 2)
final.dat.compl3 <- complete(final.dat.imp, 3)
final.dat.compl4 <- complete(final.dat.imp, 4)
final.dat.compl5 <- complete(final.dat.imp, 5)
final.dat.compl6 <- complete(final.dat.imp, 6)
final.dat.compl7 <- complete(final.dat.imp, 7)

final.dat.compl.mn <- final.dat.raw
sub1 <- which(is.na(final.dat.raw[,1]) == T)
dat1 <- data.frame(final.dat.compl1[sub1, 1], final.dat.compl2[sub1, 1], final.dat.compl3[sub1, 1], final.dat.compl4[sub1, 1], final.dat.compl5[sub1, 1], final.dat.compl6[sub1, 1], final.dat.compl7[sub1, 1])
final.dat.compl.mn[sub1, 1] <- apply(dat1, MARGIN=1, mean)

sub2 <- which(is.na(final.dat.raw[,2]) == T)
dat2 <- data.frame(final.dat.compl1[sub2, 2], final.dat.compl2[sub2, 2], final.dat.compl3[sub2, 2], final.dat.compl4[sub2, 2], final.dat.compl5[sub2, 2], final.dat.compl6[sub2, 2], final.dat.compl7[sub2, 2])
final.dat.compl.mn[sub2, 2] <- apply(dat2, MARGIN=1, mean)

sub3 <- which(is.na(final.dat.raw[,3]) == T)
dat3 <- data.frame(final.dat.compl1[sub3, 3], final.dat.compl2[sub3, 3], final.dat.compl3[sub3, 3], final.dat.compl4[sub3, 3], final.dat.compl5[sub3, 3], final.dat.compl6[sub3, 3], final.dat.compl7[sub3, 3])
final.dat.compl.mn[sub3, 3] <- apply(dat3, MARGIN=1, mean)

sub4 <- which(is.na(final.dat.raw[,4]) == T)
dat4 <- data.frame(final.dat.compl1[sub4, 4], final.dat.compl2[sub4, 4], final.dat.compl3[sub4, 4], final.dat.compl4[sub4, 4], final.dat.compl5[sub4, 4], final.dat.compl6[sub4, 4], final.dat.compl7[sub4, 4])
final.dat.compl.mn[sub4, 4] <- apply(dat4, MARGIN=1, mean)

sub5 <- which(is.na(final.dat.raw[,5]) == T)
dat5 <- data.frame(final.dat.compl1[sub5, 5], final.dat.compl2[sub5, 5], final.dat.compl3[sub5, 5], final.dat.compl4[sub5, 5], final.dat.compl5[sub5, 5], final.dat.compl6[sub5, 5], final.dat.compl7[sub5, 5])
final.dat.compl.mn[sub5, 5] <- apply(dat5, MARGIN=1, mean)

sub6 <- which(is.na(final.dat.raw[,6]) == T)
dat6 <- data.frame(final.dat.compl1[sub6, 6], final.dat.compl2[sub6, 6], final.dat.compl3[sub6, 6], final.dat.compl4[sub6, 6], final.dat.compl5[sub6, 6], final.dat.compl6[sub6, 6], final.dat.compl7[sub6, 6])
final.dat.compl.mn[sub6, 6] <- apply(dat6, MARGIN=1, mean)

sub7 <- which(is.na(final.dat.raw[,7]) == T)
dat7 <- data.frame(final.dat.compl1[sub7, 7], final.dat.compl2[sub7, 7], final.dat.compl3[sub7, 7], final.dat.compl4[sub7, 7], final.dat.compl5[sub7, 7], final.dat.compl6[sub7, 7], final.dat.compl7[sub7, 7])
final.dat.compl.mn[sub7, 7] <- apply(dat7, MARGIN=1, mean)



homevis.imp.sc <- scale(logit(final.dat.compl.mn$homevis/100), scale=T, center=T)
Feduc.imp.sc <- scale(logit(final.dat.compl.mn$Feduc/100), scale=T, center=T)
CathMusl.imp.sc <- scale(logit(final.dat.compl.mn$CathMusl/100), scale=T, center=T)
infmort.imp.sc <- scale(logit(final.dat.compl.mn$infmort/100), scale=T, center=T)
anyContr.imp.sc <- scale(logit(final.dat.compl.mn$anyContr/100), scale=T, center=T)
quality.imp.sc <- scale(logit(final.dat.compl.mn$quality/100), scale=T, center=T)
housesiz.imp.sc <- scale(logit(final.dat.compl.mn$housesiz/100), scale=T, center=T)


final.dat.imputed <- data.frame(dat$cntry.code, dat$WBregion, dat$tfr10t20, homevis.imp.sc, Feduc.imp.sc, CathMusl.imp.sc, infmort.imp.sc, anyContr.imp.sc, quality.imp.sc, housesiz.imp.sc)
colnames(final.dat.imputed) <- c("cntry.code", "reg", "fert", "homevis","Feduc","CathMusl","infmort","anyContr", "quality", "housesiz")
final.datNA <- data.frame(dat$cntry.code,dat$WBregion, dat$tfr10t20, homevst.sc, educyrsmn.sc, cathmusl.sc, imr.sc, anycontr.sc, quality.sc, meanmem.sc)
colnames(final.datNA) <- c("cntry.code","reg", "fert", "homevis","Feduc","CathMusl","infmort","anyContr", "quality", "housesiz")

## global model
## homevis, Fsece, cath, infmort, anyContr, gen3

final.dat <- final.dat.imputed
dim(final.dat)

cormat.final <- cor(na.omit(final.dat[,-c(1:3)]), method="kendall")
cormat.final
ggcorrmat(data=final.dat, cor.vars = homevis:housesiz)


# GLM
# model set
m1 <- "fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz"
m2 <- "fert ~ homevis + Feduc"
m3 <- "fert ~ homevis + CathMusl"
m4 <- "fert ~ homevis + infmort"
m5 <- "fert ~ homevis + anyContr"
m6 <- "fert ~ homevis + quality"
m7 <- "fert ~ homevis + housesiz"
m8 <- "fert ~ Feduc + infmort"
m9 <- "fert ~ Feduc + anyContr"
m10 <- "fert ~ Feduc + quality"
m11 <- "fert ~ Feduc + housesiz"
m12 <- "fert ~ CathMusl + infmort"
m13 <- "fert ~ CathMusl + anyContr"
m14 <- "fert ~ CathMusl + quality"
m15 <- "fert ~ CathMusl + housesiz"
m16 <- "fert ~ infmort + anyContr"
m17 <- "fert ~ infmort + quality"
m18 <- "fert ~ infmort + housesiz"
m19 <- "fert ~ homevis"
m20 <- "fert ~ Feduc"
m21 <- "fert ~ CathMusl"
m22 <- "fert ~ infmort"
m23 <- "fert ~ anyContr"
m24 <- "fert ~ quality"
m25 <- "fert ~ housesiz"
m26 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26)

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
brt.final <- gbm.step(final.dat, gbm.x = attr(final.dat, "names")[c(4:10)], gbm.y = attr(final.dat, "names")[3], family="gaussian", tolerance = 0.001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.final)
gbm.plot(brt.final)
tmp <- gbm.plot.fits(brt.final, v=0)
D2 <- 100 * (brt.final$cv.statistics$deviance.mean - brt.final$self.statistics$mean.resid) / brt.final$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.final$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.final$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))





## GLMM

# model set
m1 <- "fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz + (1|reg)"
m2 <- "fert ~ homevis + Feduc + (1|reg)"
m3 <- "fert ~ homevis + CathMusl + (1|reg)"
m4 <- "fert ~ homevis + infmort + (1|reg)"
m5 <- "fert ~ homevis + anyContr + (1|reg)"
m6 <- "fert ~ homevis + quality + (1|reg)"
m7 <- "fert ~ homevis + housesiz + (1|reg)"
m8 <- "fert ~ Feduc + infmort + (1|reg)"
m9 <- "fert ~ Feduc + anyContr + (1|reg)"
m10 <- "fert ~ Feduc + quality + (1|reg)"
m11 <- "fert ~ Feduc + housesiz + (1|reg)"
m12 <- "fert ~ CathMusl + infmort + (1|reg)"
m13 <- "fert ~ CathMusl + anyContr + (1|reg)"
m14 <- "fert ~ CathMusl + quality + (1|reg)"
m15 <- "fert ~ CathMusl + housesiz + (1|reg)"
m16 <- "fert ~ infmort + anyContr + (1|reg)"
m17 <- "fert ~ infmort + quality + (1|reg)"
m18 <- "fert ~ infmort + housesiz + (1|reg)"
m19 <- "fert ~ homevis + (1|reg)"
m20 <- "fert ~ Feduc + (1|reg)"
m21 <- "fert ~ CathMusl + (1|reg)"
m22 <- "fert ~ infmort + (1|reg)"
m23 <- "fert ~ anyContr + (1|reg)"
m24 <- "fert ~ quality + (1|reg)"
m25 <- "fert ~ housesiz + (1|reg)"
m26 <- "fert ~ 1 + (1|reg)"


## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]),data=final.dat, na.action=na.omit)
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


#################################
## bootstrap iteration loop to estimate prediction confidence interval
iter <- 1000
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = 0, dim = c(eq.sp.points, 7, iter), dimnames=list(paste("x",1:eq.sp.points,sep=""), attr(final.dat, "names")[c(4:10)], paste("b",1:iter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- homevis.ri <- Feduc.ri <- CathMusl.ri <- infmort.ri <- anyContr.ri <- quality.ri <- housesiz.ri <- rep(0,iter)

# iterate
for (b in 1:iter) {  # start b loop
  
  # bootstrap data
  #boot.sub <- sort(sample(x = 1:dim(final.dat)[1], size = dim(final.dat)[1], replace=TRUE))
  boot.sub <- sort(sample(x = 1:dim(final.dat)[1], size = 50, replace=TRUE))
  final.boot <- final.dat[boot.sub,]
  
  brt.fit <- gbm.step(final.boot, gbm.x = attr(final.boot, "names")[c(4:10)], gbm.y = attr(final.boot, "names")[3], family="gaussian", tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.8, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000, silent=T)
  summ.fit <- summary(brt.fit)
  
  homevis.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(4:10)][1])]
  Feduc.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(4:10)][2])]
  CathMusl.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(4:10)][3])]
  infmort.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(4:10)][4])]
  anyContr.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(4:10)][5])]
  quality.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(4:10)][6])]
  housesiz.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(final.boot, "names")[c(4:10)][7])]
  
  #gbm.plot(brt.fit)
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
  D2.vec[b] <- D2
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.vec[b] <- CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se.vec[b] <- CV.cor.se
  
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=7)
  
  ## output average predictions
  for (p in 1:7) {
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
kappa.n <- 5
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

par(mfrow=c(2,4)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←fewer) home-visiting workers (more→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←lower) female education (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←less) % Cath/Musl more→)")
lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red") 

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←lower) infant mortality (higher→)" )
lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←less) access to any contraception (more→)" )
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

plot(val.med[,6],pred.med[,6],type="l",ylim=c(min(pred.lo[,6]),max(pred.up[,6])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←less) family-planning quality (more→)" )
lines(val.med[,6], pred.lo[,6], type="l", lty=2, col="red")
lines(val.med[,6], pred.up[,6], type="l", lty=2, col="red")

plot(val.med[,7],pred.med[,7],type="l",ylim=c(min(pred.lo[,7]),max(pred.up[,7])), lwd=2, ylab="(←lower) fertility (higher →)", xlab="(←lower) household size (higher→)" )
lines(val.med[,7], pred.lo[,7], type="l", lty=2, col="red")
lines(val.med[,7], pred.up[,7], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec
CV.cor.update <- CV.cor.vec
CV.cor.se.update <- CV.cor.se.vec
homevis.ri.update <- homevis.ri
Feduc.ri.update <- Feduc.ri
CathMusl.ri.update <- CathMusl.ri
infmort.ri.update <- infmort.ri
anyContr.ri.update <- anyContr.ri
quality.ri.update <- quality.ri
housesiz.ri.update <- housesiz.ri

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  homevis.mean <- mean(homevis.ri.update, na.rm=T); homevis.sd <- sd(homevis.ri.update, na.rm=T)
  Feduc.mean <- mean(Feduc.ri.update, na.rm=T); Feduc.sd <- sd(Feduc.ri.update, na.rm=T)
  CathMusl.mean <- mean(CathMusl.ri.update, na.rm=T); CathMusl.sd <- sd(CathMusl.ri.update, na.rm=T)
  infmort.mean <- mean(infmort.ri.update, na.rm=T); infmort.sd <- sd(infmort.ri.update, na.rm=T)
  anyContr.mean <- mean(anyContr.ri.update, na.rm=T); anyContr.sd <- sd(anyContr.ri.update, na.rm=T)
  quality.mean <- mean(quality.ri.update, na.rm=T); quality.sd <- sd(quality.ri.update, na.rm=T)
  housesiz.mean <- mean(housesiz.ri.update, na.rm=T); housesiz.sd <- sd(housesiz.ri.update, na.rm=T)
  
  for (u in 1:iter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    homevis.ri.update[u] <- ifelse((homevis.ri.update[u] < (homevis.mean-kappa*homevis.sd) | homevis.ri.update[u] > (homevis.mean+kappa*homevis.sd)), NA, homevis.ri.update[u])
    Feduc.ri.update[u] <- ifelse((Feduc.ri.update[u] < (Feduc.mean-kappa*Feduc.sd) | Feduc.ri.update[u] > (Feduc.mean+kappa*Feduc.sd)), NA, Feduc.ri.update[u])
    CathMusl.ri.update[u] <- ifelse((CathMusl.ri.update[u] < (CathMusl.mean-kappa*CathMusl.sd) | CathMusl.ri.update[u] > (CathMusl.mean+kappa*CathMusl.sd)), NA, CathMusl.ri.update[u])
    infmort.ri.update[u] <- ifelse((infmort.ri.update[u] < (infmort.mean-kappa*infmort.sd) | infmort.ri.update[u] > (infmort.mean+kappa*infmort.sd)), NA, infmort.ri.update[u])
    anyContr.ri.update[u] <- ifelse((anyContr.ri.update[u] < (anyContr.mean-kappa*anyContr.sd) | anyContr.ri.update[u] > (anyContr.mean+kappa*anyContr.sd)), NA, anyContr.ri.update[u])
    quality.ri.update[u] <- ifelse((quality.ri.update[u] < (quality.mean-kappa*quality.sd) | quality.ri.update[u] > (quality.mean+kappa*quality.sd)), NA, quality.ri.update[u])
    housesiz.ri.update[u] <- ifelse((housesiz.ri.update[u] < (housesiz.mean-kappa*housesiz.sd) | housesiz.ri.update[u] > (housesiz.mean+kappa*housesiz.sd)), NA, housesiz.ri.update[u])
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

Feduc.ri.lo <- quantile(Feduc.ri.update, probs=0.025, na.rm=TRUE)
Feduc.ri.med <- median(Feduc.ri.update, na.rm=TRUE)
Feduc.ri.up <- quantile(Feduc.ri.update, probs=0.975, na.rm=TRUE)

CathMusl.ri.lo <- quantile(CathMusl.ri.update, probs=0.025, na.rm=TRUE)
CathMusl.ri.med <- median(CathMusl.ri.update, na.rm=TRUE)
CathMusl.ri.up <- quantile(CathMusl.ri.update, probs=0.975, na.rm=TRUE)

infmort.ri.lo <- quantile(infmort.ri.update, probs=0.025, na.rm=TRUE)
infmort.ri.med <- median(infmort.ri.update, na.rm=TRUE)
infmort.ri.up <- quantile(infmort.ri.update, probs=0.975, na.rm=TRUE)

anyContr.ri.lo <- quantile(anyContr.ri.update, probs=0.025, na.rm=TRUE)
anyContr.ri.med <- median(anyContr.ri.update, na.rm=TRUE)
anyContr.ri.up <- quantile(anyContr.ri.update, probs=0.975, na.rm=TRUE)

quality.ri.lo <- quantile(quality.ri.update, probs=0.025, na.rm=TRUE)
quality.ri.med <- median(quality.ri.update, na.rm=TRUE)
quality.ri.up <- quantile(quality.ri.update, probs=0.975, na.rm=TRUE)

housesiz.ri.lo <- quantile(housesiz.ri.update, probs=0.025, na.rm=TRUE)
housesiz.ri.med <- median(housesiz.ri.update, na.rm=TRUE)
housesiz.ri.up <- quantile(housesiz.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(homevis.ri.lo,Feduc.ri.lo,CathMusl.ri.lo,infmort.ri.lo,anyContr.ri.lo,quality.ri.lo,housesiz.ri.lo)
ri.med <- c(homevis.ri.med,Feduc.ri.med,CathMusl.ri.med,infmort.ri.med,anyContr.ri.med,quality.ri.med,housesiz.ri.med)
ri.up <- c(homevis.ri.up,Feduc.ri.up,CathMusl.ri.up,infmort.ri.up,anyContr.ri.up,quality.ri.up,housesiz.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.boot, "names")[c(3:9)]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort


# raw values
final.dat.raw <- data.frame(dat$cntry.code, dat$tfr10t20, dat$homevst, dat$educyrsmn, dat$ardacathmusl2020, dat$imr2020, dat$anycontr, dat$quality, dat$meanmem)
colnames(final.dat.raw) <- c("cntry.code","fert", "homevis","Feduc","CathMusl","infmort","anyContr","quality","housesiz")
brt.finalraw <- gbm.step(final.dat.raw, gbm.x = attr(final.dat.raw, "names")[c(3:9)], gbm.y = attr(final.dat.raw, "names")[2], family="gaussian", tolerance = 0.00001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
summary(brt.finalraw)
gbm.plot(brt.finalraw)
tmp <- gbm.plot.fits(brt.finalraw, v=0)
D2 <- 100 * (brt.finalraw$cv.statistics$deviance.mean - brt.finalraw$self.statistics$mean.resid) / brt.finalraw$cv.statistics$deviance.mean
D2 # % deviance explained
CV.cor <- 100 * brt.finalraw$cv.statistics$correlation.mean
CV.cor.se <- 100 *brt.finalraw$cv.statistics$correlation.se
print(c(CV.cor, CV.cor.se))

## back-transform infant mortality to look at threshold effect
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}


infmort.scale <- as.numeric(infmort.imp.sc) * attr(infmort.imp.sc, 'scaled:scale') + attr(infmort.imp.sc, 'scaled:center')
infmort.raw <- logit2prob(infmort.scale)
infmort.raw

infmort.test.out <- data.frame(final.dat.compl.mn$cntry.code, final.dat.compl.mn$infmort, final.dat.compl.mn$infmort/100,
                               infmort.imp.sc, infmort.scale, infmort.raw, final.dat.compl.mn$fert)
colnames(infmort.test.out) <- c("cntry.code", "raw", "prop", "imp.sc", "scale", "rawout","fert")
infmort.test.out

rawfert.cntry.out <- data.frame(final.dat$cntry.code, infmort.raw/10, final.dat$fert)
colnames(rawfert.cntry.out) <- c("cntry.code","infmort.raw", "fert")
head(rawfert.cntry.out)
rawfert.cntry.out

head(final.dat)


fert.pred.md <- pred.med[,4]
fert.pred.up <- pred.up[,4]
fert.pred.lo <- pred.lo[,4]
infmort.val <- val.med[,4]
infmort.scale.fit <- as.numeric(infmort.val) * attr(infmort.imp.sc, 'scaled:scale') + attr(infmort.imp.sc, 'scaled:center')
infmort.raw.fit <- logit2prob(infmort.scale.fit)

rawinfmortfertout <- data.frame(infmort.raw.fit, fert.pred.md, fert.pred.up, fert.pred.lo)


plot(imr.sc, brt.final$fitted, pch=19)
plot(infmort.raw, brt.final$fitted, pch=19, xlab="infant mortality", ylab="fertility")
lines(infmort.raw.fit,fert.pred.md,lty=1)
lines(infmort.raw.fit,fert.pred.lo,lty=2,col="red")
lines(infmort.raw.fit,fert.pred.up,lty=2,col="red")

infmort.out <- data.frame(infmort.raw,brt.final$fitted)
colnames(infmort.out) <- c("infmort","fert.fit")



##################
## GLS

# get world map
wmap <- getMap(resolution="high")

# get centroids
centroids <- gCentroid(wmap, byid=TRUE)

# get a data.frame with centroids
centroids.df <- as.data.frame(centroids)
head(centroids.df)
centroids.df$country <- (row.names(centroids.df))
colnames(centroids.df) <- c("lon","lat","cntry")
head(centroids.df)

final.dat2 <- data.frame(dat$cntry,dat$cntry.code,final.dat)
colnames(final.dat2)[1:2] <- c("cntry","cntry.code")
head(final.dat2)

# standardise country names for merging
centroids.df$cntry <- revalue(centroids.df$cntry, c("Laos" = "Lao PDR", "East Timor" = "TimorLeste", "Kyrgyzstan" = "Kyrgyz Republic",
                                                    "Turkey" = "Turkiye", "Yemen" = "Yemen Rep", "Egypt" = "Egypt Arab Rep",
                                                    "Democratic Republic of the Congo" = "Congo Dem Rep",
                                                    "Republic of the Congo" = "Congo Rep", "Ivory Coast" = "Cote dIvoire",
                                                    "Gambia" = "Gambia The", "United Republic of Tanzania" = "Tanzania"))
final.dat2

# merge
final.dat.centroids <- merge(final.dat2, centroids.df, by="cntry")
dim(final.dat.centroids)
head(final.dat.centroids)

final.dat.c <- na.omit(final.dat.centroids)
final.dat.c$x <- latlong2grid(final.dat.c[,c(10:11)])$x # equidistant coordinates
final.dat.c$y <- latlong2grid(final.dat.c[,c(10:11)])$y # equidistant coordinates
plot(final.dat.c$x,final.dat.c$y,pch=19)

## determine best correlation structure
m1 <- gls(fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz, data = final.dat.c)
vario1 <- Variogram(m1, form = ~ lon + lat, resType = "pearson")
plot(vario1, smooth = TRUE)

m2 <- gls(fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz, correlation = corExp(form = ~ lon + lat, nugget = T), data = final.dat.c)
m3 <- gls(fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz, correlation = corGaus(form = ~ lon + lat, nugget = T), data = final.dat.c)
m4 <- gls(fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz, correlation = corSpher(form = ~ lon + lat, nugget = T), data = final.dat.c)
m5 <- gls(fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz, correlation = corRatio(form = ~ lon + lat, nugget = T), data = final.dat.c)

mod.lab <- c("noCor","Exp","Gaus","Spher","Ratio")
AIC.vec <- c(AICc(m1), AICc(m2), AICc(m3), AICc(m4), AIC(m5))
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
100*(1 - results.sort[5,5]/(results.sort[1,5])) # pseudo R2 - McFadden
100*(1 - results.sort[5,6]/(results.sort[1,6])) # pseudo R2 - Cox & Snell
100*(1 - results.sort[5,7]/(results.sort[1,7])) # pseudo R2 - Craig & Uhler


vario4 <- Variogram(m2, form = ~ lon + lat, resType = "pearson")
plot(vario4, smooth = TRUE)
vario4.nr <- Variogram(m2, form = ~ lon + lat, resType = "normalized")
plot(vario4.nr, smooth = TRUE)

# run with Spherical spatial autocorrelation
# model set
m1 <- "fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz"
m2 <- "fert ~ homevis + Feduc"
m3 <- "fert ~ homevis + CathMusl"
m4 <- "fert ~ homevis + infmort"
m5 <- "fert ~ homevis + anyContr"
m6 <- "fert ~ homevis + quality"
m7 <- "fert ~ homevis + housesiz"
m8 <- "fert ~ Feduc + infmort"
m9 <- "fert ~ Feduc + anyContr"
m10 <- "fert ~ Feduc + quality"
m11 <- "fert ~ Feduc + housesiz"
m12 <- "fert ~ CathMusl + infmort"
m13 <- "fert ~ CathMusl + anyContr"
m14 <- "fert ~ CathMusl + quality"
m15 <- "fert ~ CathMusl + housesiz"
m16 <- "fert ~ infmort + anyContr"
m17 <- "fert ~ infmort + quality"
m18 <- "fert ~ infmort + housesiz"
m19 <- "fert ~ homevis"
m20 <- "fert ~ Feduc"
m21 <- "fert ~ CathMusl"
m22 <- "fert ~ infmort"
m23 <- "fert ~ anyContr"
m24 <- "fert ~ quality"
m25 <- "fert ~ housesiz"
m26 <- "fert ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
SaveCount <- BIC.vec <- AICc.vec <- LL.vec <- k.vec <- psR2mcf <- psR2cs <- psR2cu <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- gls(as.formula(mod.vec[i]), correlation = corGaus(form = ~ lon + lat, nugget = T), method="ML", data = final.dat.c)
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



## calculate difference between RSS & SS of data dispersion (another relative goodness-of-fit measure)
mdat <- mean(final.dat.c$fert)
SCR <- sum((final.dat.c$fert - mdat)**2) # sum of the residuals
nmod <- n.mod
Vargls <- matrix(0, nmod, 1)
i <- 1
Cand.models <- list( )
i

#Vargls = calculated as the difference between the sum of square of the residuals of the generalised least-square model 
#         and the sum of squares of the dispersion of the data around the observed mean divided by the sum of squares of 
#         the dispersion of the data around the observed mean.

#Mod 1
Cand.models[[i]] <- gls.1 <- gls(fert ~ Feduc + housesiz, dat=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.1)
SCR2 <- sum(gls.1$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 2
Cand.models[[i]] <- gls.2 <- gls(fert ~ Feduc + anyContr, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.2)
SCR2 <- sum(gls.2$residuals)**2
Vargls[i,1]<-100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 3
Cand.models[[i]] <- gls.3 <- gls(fert ~ infmort + housesiz, data=final.dat.c, correlation=corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.3)
SCR2 <- sum(gls.3$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 4
Cand.models[[i]] <- gls.4 <- gls(fert ~ homevis + Feduc + CathMusl + infmort + anyContr + quality + housesiz, data=final.dat.c, correlation=corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.4)
SCR2 <- sum(gls.4$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 5
Cand.models[[i]] <- gls.5 <- gls(fert ~ Feduc + infmort, data=final.dat.c, correlation=corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.5)
SCR2 <- sum(gls.5$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 6
Cand.models[[i]] <- gls.6 <- gls(fert ~ homevis + housesiz, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
#output[i,1]<-AICc(gls.1, return.K = FALSE, second.ord = TRUE, nobs = NULL)
sdat <- predict(gls.6)
SCR2 <- sum(gls.6$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 7
Cand.models[[i]] <- gls.7 <- gls(fert ~ Feduc, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.7)
SCR2 <- sum(gls.7$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 8
Cand.models[[i]] <- gls.8 <- gls(fert ~ homevis + Feduc, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.8)
SCR2 <- sum(gls.8$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 9
Cand.models[[i]] <- gls.9 <- gls(fert ~ homevis + CathMusl, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.9)
SCR2 <- sum(gls.9$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 10
Cand.models[[i]] <- gls.10 <- gls(fert ~ homevis + infmort, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.10)
SCR2 <- sum(gls.10$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 11
Cand.models[[i]] <- gls.11 <- gls(fert ~ homevis + anyContr, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.11)
SCR2 <- sum(gls.11$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 12
Cand.models[[i]] <- gls.12 <- gls(fert ~ homevis + quality, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.12)
SCR2 <- sum(gls.12$residuals)**2
Vargls[i,1]<-100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 13
Cand.models[[i]] <- gls.13 <- gls(fert ~ Feduc + quality, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.13)
SCR2 <- sum(gls.13$residuals)**2
Vargls[i,1]<-100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 14
Cand.models[[i]] <- gls.14 <- gls(fert ~ CathMusl + infmort, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.14)
SCR2 <- sum(gls.14$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 15
Cand.models[[i]] <- gls.15 <- gls(fert ~ CathMusl + anyContr, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.15)
SCR2 <- sum(gls.15$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 16
Cand.models[[i]] <- gls.16 <- gls(fert ~ CathMusl + quality, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.16)
SCR2 <- sum(gls.16$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 17
Cand.models[[i]] <- gls.17 <- gls(fert ~ CathMusl + housesiz, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.17)
SCR2 <- sum(gls.17$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 18
Cand.models[[i]] <- gls.18 <- gls(fert ~ infmort + anyContr, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <-predict(gls.18)
SCR2 <- sum(gls.18$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 19
Cand.models[[i]] <- gls.19 <- gls(fert ~ infmort + quality, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.19)
SCR2 <- sum(gls.19$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 20
Cand.models[[i]] <- gls.20 <- gls(fert ~ homevis, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.20)
SCR2 <- sum(gls.20$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 21
Cand.models[[i]] <- gls.21 <- gls(fert ~ CathMusl, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.21)
SCR2 <- sum(gls.21$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 22
Cand.models[[i]] <- gls.22 <- gls(fert ~ infmort, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.22)
SCR2 <- sum(gls.22$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 23
Cand.models[[i]] <- gls.23 <- gls(fert ~ anyContr, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.23)
SCR2 <- sum(gls.23$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 24
Cand.models[[i]] <- gls.24 <- gls(fert ~ quality, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.24)
SCR2 <- sum(gls.24$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 25
Cand.models[[i]] <- gls.25 <- gls(fert ~ housesiz, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.25)
SCR2 <- sum(gls.25$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

#Mod 26
Cand.models[[i]] <- gls.26 <- gls(fert ~ 1, data=final.dat.c, correlation= corGaus(form = ~ lon + lat, nugget = T), method="ML",verbose=T)
sdat <- predict(gls.26)
SCR2 <- sum(gls.26$residuals)**2
Vargls[i,1] <- 100*(SCR-SCR2)/SCR
i <- i + 1
i

stat <- as.data.frame(aictab(Cand.models, modnames = NULL, second.ord = TRUE, nobs = NULL,sort = F)) #stats on GLS
stat$Vargls <- as.vector(Vargls)
stat.sort <- stat[order(stat[,6],decreasing=T),1:8]
stat.sort
sumtable2 <- sumtable
sumtable2$Vargls <- as.vector(Vargls)
sumtable2$VarglsRnk <- 27 - rank(sumtable2$Vargls)
summary.table2 <- sumtable2[order(sumtable2[,9],decreasing=T),] # order by wBIC
summary.table2

