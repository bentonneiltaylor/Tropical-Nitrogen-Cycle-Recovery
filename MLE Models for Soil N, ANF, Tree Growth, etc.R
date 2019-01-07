#### Maximum Likelihood Analyses for Successional Dynamics of N fixation project #######
##These analyses deal with all data other than SNF data, which is in a separate file
setwd("C:/Users/Benton/Desktop/Work Files/Research/Costa Rica/Successional Nitrogen & Fixation Dynamics/Data Analysis and Code Files")
source("C:\\RCode\\Homegrown Functions.R")
library(bbmle)
library(multcomp)
library(lme4)
library(MASS)
library(Hmisc)
library(VGAM)

#Read in the data 
sn<-read.csv("Soil N data_Individual Cores_all data.csv")
asym<-read.csv("ANF data_Individual Cores_all data.csv")
tr.sbplt<-read.csv("Tree Growth_Subplot Summary.csv")
tr<-read.csv("Tree-based N fixation and Soil N Data.csv")
nin<-read.csv("Nitrogen inputs by Age Table.csv")
plt.alldat<-read.csv("Concatenated Plot Summary Data.csv")

#Use this code to test for changes in successional trends with different assumed ages for LEPP
#sn$age<-ifelse(sn$plot=="LEPP", yes=200, no=sn$age)
#asym$stand.age<-ifelse(asym$Plot=="LEPP", yes=200, no=asym$stand.age)
#tr$age<-ifelse(tr$plot=="LEPP", yes=200, no=tr$age)

#Writing some functions to calculate and compare AICc's
AICc <- function(fitobj,dat){
  K <- length(coef(fitobj))+1
  n <- length(dat)#make this Nrow of data
  AICc <- AIC(fitobj) + 2*K*(K+1)/(n-K-1)
}

deltaAICs <- function(AICvec){
  best <- 0
  for(i in 1:length(AICvec)){
    print(AICvec[i]-min(AICvec))
    if(AICvec[i]==min(AICvec)) best <- i
  }
  best
}

##########################################################################################
##########################################################################################
##### DOES SOIL N VARY WITH FOREST AGE??? ################################################
##########################################################################################
##########################################################################################
### Soil Nitrogen
### Creating several MLE models for how Total N varies by stand age
nit.mle<-sn[!is.na(sn$n.tot),]
nit.mle$logN<-log(nit.mle$n.tot)

#No variation############################################################################
nitmean<-function(nit_u,sigma,nit){
  nit.pred<-nit_u
  -sum(dnorm(nit, nit.pred, exp(sigma), log=T))
}
nitmean_fit<-mle2(nitmean, start=list(nit_u=3, sigma=.1), data=list(nit=nit.mle$logN))
summary(nitmean_fit)
#Linear variation############################################################################
nitage <- function(nit_u,m, c, sigma, nit, age) {
  nit.pred = age * m + c
  -sum(dnorm(nit, nit.pred, exp(sigma), log = TRUE))
}
nitage_fit<-mle2(nitage, start=list(m=.001,c=.01,sigma=.1), data=list(nit=nit.mle$logN, age=nit.mle$age))
summary(nitage_fit)
#quadratic variation############################################################################
nitquadage <- function(nit_u,b1, b2, c, sigma, nit, age) {
  nit.pred = age * b1 + b2*age^2 + c
  -sum(dnorm(nit, nit.pred, exp(sigma), log = TRUE))
}
nitquadage_fit<-mle2(nitquadage, start=list(b1=.001,b2=.05,c=3,sigma=.1), data=list(nit=nit.mle$logN, age=nit.mle$age))
summary(nitquadage_fit)
#Gaussian version 2############################################################################
lnormNLL_nit_quadage <- function(c,a,nit,age,u_nit,sig,sigma){
  nit.pred<-log(c+(a-c)*exp((-(age-u_nit)^2)/(2*sig^2)))
  - sum(dlnorm(nit,mean=log(nit.pred),sd=exp(sigma),log=TRUE))
}
fit_nit_quadage <- mle2(lnormNLL_nit_quadage,
                        start=list(c=18,a=50,u_nit=28,sig=8,sigma=.5),
                        data=list(nit=nit.mle$logN,age=nit.mle$age))
print(summary(fit_nit_quadage))

# Compare MLE models############################################################################
#################################################################################################
deltaAICs(c(AICc(nitmean_fit,nit.mle),AICc(nitage_fit,nit.mle), AICc(nitquadage_fit,nit.mle),AICc(fit_nit_quadage,nit.mle)))
### So the best fit model is one with a slight exponential decay shape. 
#This is not really significantly better than a quadratic fit with a negative quad term
#################################################################################################





##########################################################################################
##########################################################################################
##### DOES THE PROPORTION OF N AS AMMONIUM VARY WITH FOREST AGE??? #######################
##########################################################################################
##########################################################################################
sn$prop.no.tot<-with(sn, no/n.tot)
propnit.mle<-sn[!is.na(sn$prop.no.tot),]

#No variation ############################################################################
propnitmean<-function(propnit_u,sigma,propnit){
  propnit.pred<-propnit_u
  -sum(dnorm(propnit, propnit.pred, exp(sigma), log=T))
}
propnitmean_fit<-mle2(propnitmean, start=list(propnit_u=.5, sigma=.05), data=list(propnit=propnit.mle$prop.no.tot))
summary(propnitmean_fit)
#Linear variation ############################################################################
propnitage <- function(propnit_u,m, c, sigma, propnit, age) {
  propnit.pred = age * m + c
  -sum(dnorm(propnit, propnit.pred, exp(sigma), log = TRUE))
}
propnitage_fit<-mle2(propnitage, start=list(m=.001,c=.01,sigma=.1), data=list(propnit=propnit.mle$prop.no.tot, age=propnit.mle$age))
summary(propnitage_fit)
#quadratic variation ############################################################################
propnitquadage <- function(propnit_u,b1, b2, c, sigma, propnit, age) {
  propnit.pred = age * b1 + b2*age^2 + c
  -sum(dnorm(propnit, propnit.pred, exp(sigma), log = TRUE))
}
propnitquadage_fit<-mle2(propnitquadage, start=list(b1=.001,b2=.05,c=.01,sigma=.1), data=list(propnit=propnit.mle$prop.no.tot, age=propnit.mle$age))
summary(propnitquadage_fit)
#inverse sigmoidal variation ############################################################################
propnitsigage <- function(propnit_u,b1, b2, b3,sigma, propnit, age) {
  propnit.pred = b1/(1+exp(b2*age-b3))
  -sum(dnorm(propnit, propnit.pred, exp(sigma), log = TRUE))
}
propnitsigage_fit<-mle2(propnitsigage, start=list(b1=.9,b2=.05,b3=35,sigma=.1), data=list(propnit=propnit.mle$prop.no.tot, age=propnit.mle$age))
summary(propnitsigage_fit)
#Gaussian version 2  ############################################################################
lnormNLL_propnit_quadage2 <- function(c,a,propnit,age,u_nit,sig,sigma){
  propnit.pred<-log(c+(a-c)*exp((-(age-u_nit)^2)/(2*sig^2)))
  - sum(dlnorm(propnit,mean=log(propnit.pred),sd=exp(sigma),log=TRUE))
}
fit_propnit_quadage2 <- mle2(lnormNLL_propnit_quadage2,
                            start=list(c=.1,a=.2,u_nit=28,sig=.5,sigma=.05),
                            data=list(propnit=log(1-propnit.mle$prop.no.tot),age=propnit.mle$age))
print(summary(fit_propnit_quadage2))

##########################################################################################################
deltaAICs(c(AICc(propnitmean_fit,propnit.mle),AICc(propnitage_fit,propnit.mle), AICc(propnitquadage_fit,propnit.mle),AICc(propnitsigage_fit,propnit.mle)))
#It looks like the contribution of ammonium to total N is high in all plots, but declines through
#succession in a quadratic nature. 
##########################################################################################################






##########################################################################################
##########################################################################################
##### DOES THE SOIL AMMONIUM VARY WITH FOREST AGE??? #####################################
##########################################################################################
##########################################################################################
nh.mle<-tr[!is.na(tr$nh),]
nh.mle$logNH<-log(nh.mle$nh)

#No variation #############################################################################
nhmean<-function(nh_u,sigma,nh){
  nh.pred<-nh_u
  -sum(dnorm(nh, nh.pred, exp(sigma), log=T))
}
nhmean_fit<-mle2(nhmean, start=list(nh_u=3, sigma=.1), data=list(nh=nh.mle$logNH))
summary(nhmean_fit)
#Linear variation #############################################################################
nhage <- function(nh_u,m, c, sigma, nh, age) {
  nh.pred = age * m + c
  -sum(dnorm(nh, nh.pred, exp(sigma), log = TRUE))
}
nhage_fit<-mle2(nhage, start=list(m=.001,c=.01,sigma=.1), data=list(nh=nh.mle$logNH, age=nh.mle$age))
summary(nhage_fit)
#quadratic variation #############################################################################
nhquadage <- function(nh_u,b1, b2, c, sigma, nh, age) {
  nh.pred = age * b1 + b2*age^2 + c
  -sum(dnorm(nh, nh.pred, exp(sigma), log = TRUE))
}
nhquadage_fit<-mle2(nhquadage, start=list(b1=.001,b2=.05,c=.01,sigma=.1), data=list(nh=nh.mle$logNH, age=nh.mle$age))
summary(nhquadage_fit)
#Gaussian version 2 #############################################################################
lnormNLL_nh_gaus <- function(c,a,nh,age,u_nit,sig,sigma){
  nh.pred<-log(c+(a-c)*exp((-(age-u_nit)^2)/(2*sig^2)))
  - sum(dlnorm(nh,mean=log(nh.pred),sd=exp(sigma),log=TRUE))
}
fit_nh_gaus <- mle2(lnormNLL_nh_gaus,
                    start=list(c=18,a=50,u_nit=28,sig=8,sigma=.5),
                    data=list(nh=nh.mle$logNH, age=nh.mle$age))
print(summary(fit_nh_gaus))

#################################################################################################
deltaAICs(c(AICc(nhmean_fit,nh.mle),AICc(nhage_fit,nh.mle), AICc(nhquadage_fit,nh.mle),AICc(fit_nh_gaus,nh.mle)))
#Linear fit with a negative slope correlating nh4 availability with age is the best fit.
##################################################################################################







##########################################################################################
##########################################################################################
##### DOES THE SOIL NITRATE VARY WITH FOREST AGE??? ######################################
##########################################################################################
##########################################################################################
no.mle<-sn[!is.na(sn$no),]
no.mle$logNO<-log(no.mle$no+1)
#No variation
nomean<-function(no_u,sigma,no){
  no.pred<-no_u
  -sum(dnorm(no, no.pred, exp(sigma), log=T))
}
nomean_fit<-mle2(nomean, start=list(no_u=3, sigma=.1), data=list(no=no.mle$logNO))
summary(nomean_fit)#remember this is ***log+1*** so the mean needs to be back-transformed
#Linear variation #########################################################################
noage <- function(no_u,m, c, sigma, no, age) {
  no.pred = age * m + c
  -sum(dnorm(no, no.pred, exp(sigma), log = TRUE))
}
noage_fit<-mle2(noage, start=list(m=.001,c=.01,sigma=.1), data=list(no=no.mle$logNO, age=no.mle$age))
summary(noage_fit)
#quadratic variation #########################################################################
noquadage <- function(no_u,b1, b2, c, sigma, no, age) {
  no.pred = age * b1 + b2*age^2 + c
  -sum(dnorm(no, no.pred, exp(sigma), log = TRUE))
}
noquadage_fit<-mle2(noquadage, start=list(b1=.001,b2=.05,c=.01,sigma=.1), data=list(no=no.mle$logNO, age=no.mle$age))
summary(noquadage_fit)
#Gaussian version 2 #########################################################################
#lnormNLL_no_gaus <- function(c,a,no,age,u_nit,sig,sigma){
#  no.pred<-log(c+(a-c)*exp((-(age-u_nit)^2)/(2*sig^2)))
#  - sum(dlnorm(no,mean=log(no.pred),sd=exp(sigma),log=TRUE))
#}
#fit_no_gaus <- mle2(lnormNLL_no_gaus,
#                    start=list(c=2,a=2.0,u_nit=28,sig=28,sigma=.5),
#                    data=list(no=no.mle$logNO, age=no.mle$age))
#print(summary(fit_no_gaus))
#############################################################################################
deltaAICs(c(AICc(nomean_fit,no.mle),AICc(noage_fit,no.mle), AICc(noquadage_fit,no.mle)))
#The quadratic fit with a slight positive slope and a negative quadratic term (hump-shape) is the best fit.
#############################################################################################







##########################################################################################
##########################################################################################
##### DOES ANF VARY WITH FOREST AGE??? ###################################################
##########################################################################################
##########################################################################################
asym$loganf<-log(asym$Nfixd.ha)

#No variation ###########################################################################
anfmean<-function(anf_u,sigma,anf){
  anf.pred<-anf_u
  -sum(dnorm(anf, anf.pred, exp(sigma), log=T))
}
anfmean_fit<-mle2(anfmean, start=list(anf_u=1, sigma=.1), data=list(anf=asym$loganf))
summary(anfmean_fit)
#Linear variation ###########################################################################
anfage <- function(anf_u,m, c, sigma, anf, age) {
  anf.pred = age * m + c
  -sum(dnorm(anf, anf.pred, exp(sigma), log = TRUE))
}
anfage_fit<-mle2(anfage, start=list(m=.001,c=.01,sigma=.1), data=list(anf=asym$loganf, age=asym$stand.age))
summary(anfage_fit)
#quadratic variation ###########################################################################
anfquadage <- function(anf_u,b1, b2, c, sigma, anf, age) {
  anf.pred = age * b1 + b2*age^2 + c
  -sum(dnorm(anf, anf.pred, exp(sigma), log = TRUE))
}
anfquadage_fit<-mle2(anfquadage, start=list(b1=.001,b2=.05,c=.01,sigma=.1), data=list(anf=asym$loganf, age=asym$stand.age))
summary(anfquadage_fit)
#Gaussian version 2 ###########################################################################
lnormNLL_anf_quadage <- function(c,a,age,anf_u,sig,sigma){
  anf.pred<-exp(c+(a-c)*exp((-(age-anf_u)^2)/(2*sig^2)))
  - sum(dlnorm(anf,mean=anf.pred,sd=exp(sigma),log=TRUE))
}
fit_anf_quadage <- mle2(lnormNLL_anf_quadage,
                        start=list(c=2,a=10,anf_u=35,sig=8,sigma=5),
                        data=list(anf=asym$Nfixd.ha, age=asym$stand.age))
print(summary(fit_anf_quadage))
###########################################################################################################
deltaAICs(c(AICc(anfmean_fit,asym),AICc(anfage_fit,asym), AICc(anfquadage_fit,asym),AICc(fit_anf_quadage,asym)))
### There's an issue here somewhere. The Gaussian is the best fit when looking at the plots
### but the deltaAIC doesn't say that. Need to figure this out.
###########################################################################################################









##########################################################################################
##########################################################################################
##### DO ANF INCUBATION RATES VARY WITH FOREST AGE??? ####################################
##########################################################################################
##########################################################################################

#No variation ############################################################################
rate.anfmean<-function(rate.anf_u,sigma,rate.anf){
  rate.anf.pred<-rate.anf_u
  -sum(dnorm(rate.anf, rate.anf.pred, exp(sigma), log=T))
}
rate.anfmean_fit<-mle2(rate.anfmean, start=list(rate.anf_u=.5, sigma=.1), data=list(rate.anf=asym$XNfixd))
summary(rate.anfmean_fit)
#Linear variation ############################################################################
rate.anfage <- function(rate.anf_u,m, c, sigma, rate.anf, age) {
  rate.anf.pred = age * m + c
  -sum(dnorm(rate.anf, rate.anf.pred, exp(sigma), log = TRUE))
}
rate.anfage_fit<-mle2(rate.anfage, start=list(m=.001,c=.01,sigma=.1), data=list(rate.anf=asym$XNfixd, age=asym$stand.age))
summary(rate.anfage_fit)
#quadratic variation ############################################################################
rate.anfquadage <- function(rate.anf_u,b1, b2, c, sigma, rate.anf, age) {
  rate.anf.pred = age * b1 + b2*age^2 + c
  -sum(dnorm(rate.anf, rate.anf.pred, exp(sigma), log = TRUE))
}
rate.anfquadage_fit<-mle2(rate.anfquadage, start=list(b1=.001,b2=.05,c=.01,sigma=.1), data=list(rate.anf=asym$XNfixd, age=asym$stand.age))
summary(rate.anfquadage_fit)
###################################################################################################
deltaAICs(c(AICc(rate.anfmean_fit,asym),AICc(rate.anfage_fit,asym), AICc(rate.anfquadage_fit,asym)))
#The actual rate of ANF per unit leaf-litter mass increases in a hump-shaped fashion throug succession
#similar to the way area-based ANF rates vary through succession
###################################################################################################







##########################################################################################
##########################################################################################
##### DOES LITTER DENSITY VARY WITH FOREST AGE??? ########################################
##########################################################################################
##########################################################################################

#No variation ############################################################################
ms.anfmean<-function(ms.anf_u,sigma,ms.anf){
  ms.anf.pred<-ms.anf_u
  -sum(dnorm(ms.anf, ms.anf.pred, exp(sigma), log=T))
}
ms.anfmean_fit<-mle2(ms.anfmean, start=list(ms.anf_u=2, sigma=.5), data=list(ms.anf=asym$tot.mass))
summary(ms.anfmean_fit)
#Linear variation ############################################################################
ms.anfage <- function(ms.anf_u,m, c, sigma, ms.anf, age) {
  ms.anf.pred = age * m + c
  -sum(dnorm(ms.anf, ms.anf.pred, exp(sigma), log = TRUE))
}
ms.anfage_fit<-mle2(ms.anfage, start=list(m=.001,c=.01,sigma=.1), data=list(ms.anf=asym$tot.mass, age=asym$stand.age))
summary(ms.anfage_fit)
#quadratic variation ############################################################################
ms.anfquadage <- function(ms.anf_u,b1, b2, c, sigma, ms.anf, age) {
  ms.anf.pred = age * b1 + b2*age^2 + c
  -sum(dnorm(ms.anf, ms.anf.pred, exp(sigma), log = TRUE))
}
ms.anfquadage_fit<-mle2(ms.anfquadage, start=list(b1=.001,b2=.05,c=.01,sigma=.1), data=list(ms.anf=asym$tot.mass, age=asym$stand.age))
summary(ms.anfquadage_fit)
###################################################################################################
deltaAICs(c(AICc(ms.anfmean_fit,asym),AICc(ms.anfage_fit,asym), AICc(ms.anfquadage_fit,asym)))
#Leaf litter mass per unit ground area does not change significantly through succession in these plots
###################################################################################################










##########################################################################################
##########################################################################################
##### DOES TREE GROWTH VARY WITH FOREST AGE??? ########################################
##########################################################################################
##########################################################################################
#We'll look at the subplot scale (to account for within-plot variation)

#No variation ###########################################################################
grwmean<-function(grw_u,sigma,grw){
  grw.pred<-grw_u
  -sum(dnorm(grw, grw.pred, exp(sigma), log=T))
}
grwmean_fit<-mle2(grwmean, start=list(grw_u=100, sigma=100), data=list(grw=tr.sbplt$BAI/100))
summary(grwmean_fit)
#Linear variation ###########################################################################
grwage <- function(grw_u,m, c, sigma, grw, age) {
  grw.pred = age * m + c
  -sum(dnorm(grw, grw.pred, exp(sigma), log = TRUE))
}
grwage_fit<-mle2(grwage, start=list(m=-10,c=100,sigma=10), data=list(grw=(tr.sbplt$BAI/100),age=tr.sbplt$stand.age))
summary(grwage_fit)
#quadratic variation ###########################################################################
grwage.quad <- function(grw_u,b1, b2, c, sigma, grw,age) {
  grw.pred = age * b1 + b2*age^2 + c
  -sum(dnorm(grw, grw.pred, exp(sigma), log = TRUE))
}
grwage.quad_fit<-mle2(grwage.quad, start=list(b1=-10,b2=-1,c=100,sigma=10), data=list(grw=(tr.sbplt$BAI/100),age=tr.sbplt$stand.age))
summary(grwage.quad_fit)
#sigmoidal decay ###########################################################################
grwage.sig <- function(grw_u,b1, b2, b3,sigma, grw,age) {
  grw.pred = b1/(1+exp(age - b2))+b3
  -sum(dnorm(grw, grw.pred, exp(sigma), log = TRUE))
}
grwage.sig_fit<-mle2(grwage.sig, start=list(b1=1.15,b2=40,b3=1,sigma=10), data=list(grw=(tr.sbplt$BAI/100),age=tr.sbplt$stand.age))
summary(grwage.sig_fit)
#################################################################################################
deltaAICs(c(AICc(grwmean_fit,tr.sbplt),AICc(grwage_fit,tr.sbplt), AICc(grwage.quad_fit,tr.sbplt),AICc(grwage.sig_fit,tr.sbplt)))
#The linear decay is the best model based on the delta AIC's (but the sigmoidal fit sure looks better)
#################################################################################################






##########################################################################################
##########################################################################################
##### DOES THE RELATIVE CONTRIBUTION OF ANF VS SNF VARY WITH FOREST AGE??? ###############
##########################################################################################
##########################################################################################

nin[,c(1,2,5,10:12)]
plt.alldat[,c(1,8,15,22,41)]







##########################################################################################
##########################################################################################
###############################################################################################
#### Does SNF change with relative abundance of N fixers? #####################################
###############################################################################################
##########################################################################################
##########################################################################################

#We'll look at the correlation of SNF with N fixer abundance first at the core level
with(tr, summary(lm(logsnf~fixNCI)))#fixer NCI around each core doesn't predict snf rates (P = .1798)

### Creating several MLE models for how SNF varies by N fixer abundance
#No variation
trabun<-tr[!is.na(tr$fixNCI),]

snfmean<-function(snf_u,sigma,snf){
  snf.pred<-snf_u
  -sum(dnorm(snf, snf.pred, exp(sigma), log=T))
}
snfmean_fit<-mle2(snfmean, start=list(snf_u=1, sigma=.1), data=list(snf=trabun$logsnf))
summary(snfmean_fit)
#Linear variation
snfabun <- function(snf_u,m, c, sigma, snf, abun) {
  snf.pred = abun * m + c
  -sum(dnorm(snf, snf.pred, exp(sigma), log = TRUE))
}
snfabun_fit<-mle2(snfabun, start=list(m=.00008,c=.02,sigma=1), data=list(snf=trabun$logsnf, abun=trabun$fixNCI))
summary(snfabun_fit)
#quadratic variation
snfquadabun <- function(snf_u,b1, b2, c, sigma, snf, abun) {
  snf.pred = abun * b1 + b2*abun^2 + c
  -sum(dnorm(snf, snf.pred, exp(sigma), log = TRUE))
}
snfquadabun_fit<-mle2(snfquadabun, start=list(b1=.001,b2=.05,c=.01,sigma=.1), data=list(snf=trabun$logsnf, abun=trabun$fixNCI))
summary(snfquadabun_fit)

deltaAICs(c(AIC(snfmean_fit),AIC(snfabun_fit), AIC(snfquadabun_fit)))

### Creating several MLE models for how SNF varies by N fixer abundance
#No variation
snfmean<-function(snf_u,sigma,snf){
  snf.pred<-snf_u
  -sum(dnorm(snf, snf.pred, exp(sigma), log=T))
}
snfmean_fit<-mle2(snfmean, start=list(snf_u=1, sigma=.1), data=list(snf=tr$logsnf))
summary(snfmean_fit)
#Linear variation
snfabun.tr <- function(snf_u,m, c, sigma, snf, abun.tr) {
  snf.pred = abun.tr * m + c
  -sum(dnorm(snf, snf.pred, exp(sigma), log = TRUE))
}
snfabun.tr_fit<-mle2(snfabun.tr, start=list(m=.00008,c=.02,sigma=1), data=list(snf=tr$logsnf, abun.tr=tr$tr.fixNCI))
summary(snfabun.tr_fit)
#quadratic variation
snfquadabun.tr <- function(snf_u,b1, b2, c, sigma, snf, abun.tr) {
  snf.pred = abun.tr * b1 + b2*abun.tr^2 + c
  -sum(dnorm(snf, snf.pred, exp(sigma), log = TRUE))
}
snfquadabun.tr_fit<-mle2(snfquadabun.tr, start=list(b1=.001,b2=.05,c=.01,sigma=.1), data=list(snf=tr$logsnf, abun.tr=tr$tr.fixNCI))
summary(snfquadabun.tr_fit)

deltaAICs(c(AIC(snfmean_fit),AIC(snfabun.tr_fit), AIC(snfquadabun.tr_fit)))

#Now lets see if snf is correlated with the fixer NCI of the tree the core was taken around
with(tr, summary(lm(logsnf~tr.fixNCI)))#fixer NCI around the focal tree doesn't predict snf rates (P=.711)
#And we'll check to see if snf is correlated with the latest growth data of the focal tree
with(tr, summary(lm(logsnf~tr.Growth)))#Nope, no correlation between snf and focal tree growth (P=.631)

#Does plot-level BAI relate to SNF?
with(plt.alldat, summary(lm(BAI~snf.kg.ha))) #No, P=.3386
#Does plot-level N-fixer BAI relate to SNF?
with(plt.alldat, summary(lm(fixerBAI~snf.kg.ha))) #No, P=.3386

#Finally, let's see if plot-level SNF is correlated with plot fixer basal area
with(plt.alldat, summary(lm(snf.kg.ha~fixerBA)))#The basal area of fixers in a plot is significantly NEGATIVELY correlated with snf


#Checking to see if the nodulation itself (hasnods) or the amount of nodules for those nodulated individuals varies with fixer abundance
#testing for nodulation as a binary variable
nodulation.core<-glmer(hasnods~fixNCI+(1|Sample.Tree)+(1|plot), data=tr, family=binomial(link="logit"),nAGQ=1)
summary(nodulation.core)#No, fixer NCI around the core is not a significant predictor of whether a core will have nodules
nodulation.tr<-glmer(hasnods~tr.fixNCI+(1|plot), data=tr, family=binomial(link="logit"),nAGQ=10)
summary(nodulation.tr)#No, fixer crowding around the focal tree is not a significant predictor of whether a core will have nodules

#Now making a dataframe of only cores with nodules to test for effects of fixer abundance on the amount of nodules
nodulated<-tr[tr$nod.bio>0,]
with(nodulated, summary(lm(logsnf~fixNCI)))#For cores with nodules, crowding of fixers around the core isn't a predictor of snf
with(nodulated, summary(lm(logsnf~tr.fixNCI)))#For cores with nodules, crowding of fixers around the focal tree isn't a predictor of snf

#Does SNF vary with the density of N-fixer stems in a plot or the ratio of fixers to non-fixers?
with(plt.alldat, summary(lm(snf.g.m2~fix.den)))#Not significant (P=.345)
with(plt.alldat, summary(lm(snf.g.m2~fnf.ratio)))#Not significant (P=.227)

###############################################################################
###############################################################################
###################### END ####################################################
###############################################################################
###############################################################################