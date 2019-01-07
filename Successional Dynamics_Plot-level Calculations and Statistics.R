### Analyses for Successional Dynamics of Nitrogen in La Selva Bosques Plots ###
#Ben Taylor
#Code Started 11/20/2017

setwd("C:/Users/Benton/Desktop/Work Files/Research/Costa Rica/Successional Nitrogen & Fixation Dynamics/Data Analysis and Code Files")
source("C:\\RCode\\Homegrown Functions.R")
library(bbmle)
library(multcomp)
library(lme4)
library(MASS)
library(Hmisc)
library(VGAM)
#################################################################################
############# Tree Growth #######################################################
#################################################################################
tr.plt<-read.csv("Tree Growth_Plot Summary.csv")
tr.sbplt<-read.csv("Tree Growth_Subplot Summary.csv")

#################################################################################
############# Asymbiotic Nitrogen Fixation ######################################
#################################################################################
asym<-read.csv("Bosques Asymbiotic Fixation Data_2016.csv")
asym$stand.age<-c(rep(19,16),rep(100,8),rep(37,8),rep(29,8))
#area of sample in m2
area<-(3.14*(.04)^2)
#Jar Volume replaced (proportion), Enrichment of normal atmosphere, enrichment of added gas
jvr<-.5
X15Nair<-.003663
X15Ngas<-.98
#Expected Enrichment of Jar Headspace during Incubation
enrich<-(X15Ngas*jvr)+(X15Nair*(1-jvr))
#Now the Percent of N atoms fixed in each sample calculated by 
#taking the difference between the measured X15N and X15N of fixation and dividing by expected enrichment
asym$XNfixd<-with(asym, ((X15Nat/100)-.003663)/(enrich-.003663))
#Calculating the total amount of N in each sample
asym$tot.N<-with(asym, tot.mass*(X.N/100))
#calculating the amount of N fixed in each sample (gN fixed per jar area per incubation period)
asym$Nfixd<-with(asym, tot.N*XNfixd)
#Now the amount of N fixed in g per m2 per year
asym$Nfixd.m2<-with(asym, (Nfixd/area)*365)
#Now in kgN/ha/yr
asym$Nfixd.ha<-with(asym, (Nfixd.m2*10))
#write.csv(asym, file="ANF data_Individual Cores_all data.csv", row.names=F)

asym.plt<-data.frame("Plot"=unique(asym$Plot),
                     "stand.age"=c(19,19,100,37,29),
                     "Nfixd.m2"=with(asym, tapply(Nfixd.m2, Plot, mean)),
                     "Nfixd.m2.se"=with(asym, tapply(Nfixd.m2, Plot, sefun)),
                     "Nfixd.ha"=with(asym, tapply(Nfixd.ha, Plot, mean)),
                     "Nfixd.ha.se"=with(asym, tapply(Nfixd.ha, Plot, sefun)),
                     "Lit.mass"=with(asym, tapply(tot.mass, Plot, mean)),
                     "Lit.mass.se"=with(asym, tapply(tot.mass, Plot, sefun)),
                     "fix.rt"=with(asym, tapply(XNfixd, Plot, mean)),
                     "fix.rt.se"=with(asym, tapply(XNfixd, Plot, sefun)))
asym.plt.fudge<-asym.plt
asym.plt.fudge$stand.age<-as.factor(ifelse(asym.plt$Plot=="JE", yes=20, no=as.character(asym.plt$stand.age)))
#write.csv(asym.plt, file="Asymbiotic Fixation_Plot Summary Data.csv",row.names=F)
#write.csv(asym.plt.fudge, file="Asymbiotic Fixation_Plot Summary_Year Fudge.csv", row.names=F)

#################################################################################
############# Nodule Biomass ####################################################
#################################################################################
### Read in the dataframes, one for tree-based and one for plot-based sampling
adlt<-read.csv("Tree-based Nodule Data_2017.csv")
adlt<-adlt[adlt$tree.id!="3-063-000006",]
grw.mrg<-read.csv("Sample Tree Growth Data_2014.csv")
adlt<-merge(adlt, grw.mrg, by="tree.id", all.x=T)
plt<-read.csv("Bosques Nodule Biomass_Data_2016.csv")
plt.age<-unique(plt[,c("Plot","Stand.Age")]) #creating a list of plots and their ages to add back into the combined "nod" dataframe later
colnames(plt.age)<-c("plot","age")
plt.nci<-read.csv("Plot-based nodule NCI_to merge.csv")
plt.nci$mrg.col<-with(plt.nci, paste(Locality_Code,"_",X,"-",Y))
plt$mrg.col<-with(plt, paste(Plot,"_",X,"-",Y))
plt.nci<-plt.nci[,c(10,7,8)]
plt<-merge(plt[,-10],plt.nci,by="mrg.col",all.x=T)
### Because fixer NCI doesn't affect Nodule Biomass, now we'll combine both sampling-scheme datasets
nod<-data.frame("plot"=as.factor(c(as.character(plt$Plot),as.character(adlt$plot))),
                "core"=as.factor(c(as.character(plt$Ind_ID),as.character(adlt$sample.id))),
                "X"=c(plt$X,adlt$core.x),
                "Y"=c(plt$Y,adlt$core.y),
                "nod.num"=c(plt$Nodule.Number,adlt$nodule.num),
                "nod.bio"=c(plt$Nodule.Biomass,adlt$nodule.mass..g.),
                "NCI"=c(plt$NCI,adlt$NCI),
                "fixNCI"=c(plt$fixNCI,adlt$fixNCI),
                "tr.DAP"=c(rep(NA,nrow(plt)),adlt$tr.DAP),
                "tr.height"=c(rep(NA,nrow(plt)),adlt$tr.height),
                "tr.NCI"=c(rep(NA,nrow(plt)),adlt$tr.NCI),
                "tr.fixNCI"=c(rep(NA,nrow(plt)),adlt$tr.fixNCI),
                "tr.nh"=c(rep(NA,nrow(plt)),adlt$tr.nh),
                "tr.RH"=c(rep(NA,nrow(plt)),adlt$tr.RH))
nod<-merge(nod, plt.age,by="plot",all.x=T,all.y=T)
nod$nod.size.mg<-with(nod, (nod.bio/nod.num)*1000) #calculates the average biomass (in mg) for each nodule 
### Calculate nodule biomass per m2 (nodule biomass for each core divided by core area in m2)
core.ar<-(pi*(4^2))*.0001  #This is the total surface area (in m2) cored around each tree (sum of 3 8-cm diameter cores)
core.ar.tot<-core.ar*80 #Gives the total area sampled per plot (core area x number of cores per plot (80))
nod$nod.bio.m2<-nod$nod.bio/core.ar

noninfl.gmn<-function(vec) exp(mean(log(vec[vec>0])))
gmn<-function(vec) (1-(length(vec[vec==0])/length(vec)))*exp(mean(log(vec[vec>0])))
gse.hi<-function(vec) (1-(length(vec[vec==0])/length(vec)))*(exp(mean(log(vec[vec>0]))+sefun(log(vec[vec>0]))))
gse.lo<-function(vec) (1-(length(vec[vec==0])/length(vec)))*(exp(mean(log(vec[vec>0]))-sefun(log(vec[vec>0]))))

### Making a dataframe for nodulation in each plot
nod.plt<-data.frame("plot"=sort(unique(nod$plot)),#Each plot name
                    "age"=plt.age$age,#Plot stand age
                    "nod.bio.gmn"=with(nod, tapply(nod.bio.m2, plot, gmn)),
                    "nod.bio.gse.hi"=with(nod, tapply(nod.bio.m2, plot, gse.hi)),
                    "nod.bio.gse.lo"=with(nod, tapply(nod.bio.m2, plot, gse.lo)),
                    "nod.bio"=with(nod, tapply(nod.bio.m2, plot, mean)),#mean nodule biomass (g/m2) for each plot
                    "nod.bio.se"=with(nod, tapply(nod.bio.m2, plot, sefun)),#standard error for nodule biomass (g/m2) for each plot
                    "nod.num"=with(nod, tapply(nod.num, plot, sum)),#total number of individual nodules found in each plot
                    "nod.size"=with(nod, tapply(nod.size.mg, plot, mean, na.rm=T)),#mean mass of each individual nodule (in mg). This comes from the mean nodule size in each core then averaged for each plot
                    "nod.size.se"=with(nod, tapply(nod.size.mg, plot, sefun)))#standard error for nodule size. This is the standard error of core averages for each plot
nod.plt.fudge<-nod.plt
nod.plt.fudge$age<-ifelse(nod.plt.fudge$plot=="JE", yes=20, no=nod.plt.fudge$age)
#write.csv(nod.plt, file="Nodule Biomass_Plot Summary Data.csv",row.names=F)
#write.csv(nod.plt.fudge, file="Nodule Biomass_Plot Summary_Year Fudge.csv", row.names=F)

#################################################################################
############# Symbiotic Nitrogen Fixation #######################################
#################################################################################
# First we need to calculate fixation rates for the nodules we did 15N incubations on
# We'll get these fixation rates in units of gN/g nodule/yr so we can multiply it 
# by g nodules in each sample that we took from inside the plots

# N fixation in incubation samples
sym<-read.csv("Bosques Nodule Incubation Data_2016.csv")
#Syringe Volume replaced (proportion), Enrichment of normal atmosphere, enrichment of added gas
svr<-.2
X15Nair<-.003663
X15Ngas<-.98
#Expected Enrichment of syringe Headspace during Incubation
syr.enrich<-(X15Ngas*svr)+(X15Nair*(1-svr))
#Now the Percent of N atoms fixed in each sample calculated by 
#taking the difference between the measured X15N and X15N of fixation and dividing by expected enrichment
sym$XNfixd<-with(sym, ((X15Nat/100)-.003663)/(syr.enrich-.003663))
#Calculating the total amount of N in each sample
sym$tot.N<-with(sym, nod.mass*(X.N/100))
#calculating the amount of N fixed in each sample (gN fixed per sample per incubation period)
sym$Nfixd.nod<-with(sym, tot.N*XNfixd)
#Now we'll convert N fixation to units of (gN fixed per gram of nodule per year)
time.conv<-17520 #the number of "half hours" in a year
sym$Nfixd<-with(sym, (Nfixd.nod/nod.mass)*time.conv)
#Putting together a plot summary data frame to apply to the nodule biomass data above
nod.inc<-data.frame("plot"=sort(unique(sym$Plot)),
                    "Nfix.rate"=with(sym, tapply( Nfixd, Plot, mean)),#This is in units of g N fixed per g nodule per year
                    "Nfix.rate.sd"=with(sym, tapply(Nfixd, Plot, sd)),#This is in units of g N fixed per g nodule per year
                    "Nfix.rate.se"=with(sym, tapply(Nfixd, Plot, sefun)))#This is in units of g N fixed per g nodule per year
#Merging the nodule biomass and fixation rate data
nod<-merge(nod, nod.inc, by="plot", all.x=T)
nod$Nfixd.core<-nod$nod.bio*nod$Nfix.rate #N fixation in gN/core/yr
nod$Nfixd.g.m2<-nod$nod.bio.m2*nod$Nfix.rate #N fixation estimate in gN/m2/yr
nod$Nfixd.kg.ha<-(nod$Nfixd.g.m2/1000)*10000 #N fixation estimate in kgN/ha/yr
nod$logsnf<-log(nod$Nfixd.kg.ha+1)
#write.csv(nod, file="Plot- and tree-based SNF data combined.csv", row.names=F)

#This creates a function to calculate the arithmetic mean of zero-inflated lognormal distributed data.
#It takes a vector of data that are not log-transformed (i.e. raw SNF data) including all the zeros
zinflamn<-function(vec) (1-(length(vec[vec==0])/length(vec)))*exp(mean(log(vec[vec>0])) + (sd(log(vec[vec>0]))*sd(log(vec[vec>0])))/2)

### Making a dataframe for Symbiotic Fixation in each plot
sym.plt<-data.frame("plot"=sort(unique(nod$plot)),#Each plot name
                    "age"=plt.age$age,#Plot stand age
                    #"snf.g.m2.zln"=with(nod, tapply(Nfixd.g.m2, plot, zinflgmn)),#mean SNF (g/m2) for each plot
                    "snf.g.m2.gmn"=with(nod, tapply(Nfixd.g.m2, plot, gmn)),
                    "snf.g.m2.gse.hi"=with(nod, tapply(Nfixd.g.m2, plot, gse.hi)),
                    "snf.g.m2.gse.lo"=with(nod, tapply(Nfixd.g.m2, plot, gse.lo)),
                    "snf.g.m2"=with(nod, tapply(Nfixd.g.m2, plot, mean)),
                    "snf.g.m2.sd"=with(nod, tapply(Nfixd.g.m2, plot, sd)),
                    "snf.g.m2.se"=with(nod, tapply(Nfixd.g.m2, plot, sefun)),#standard error for SNF (g/m2) for each plot
                    "snf.kg.ha.zln"=with(nod, tapply(Nfixd.kg.ha, plot, zinflamn)),#0-infl lognormal mean SNF (kg/ha) for each plot
                    "snf.kg.ha.gmn"=with(nod, tapply(Nfixd.kg.ha, plot, gmn)),#0-infl geometric mean SNF for each plot
                    "snf.kg.ha.gse.hi"=with(nod, tapply(Nfixd.kg.ha, plot, gse.hi)),
                    "snf.kg.ha.gse.lo"=with(nod, tapply(Nfixd.kg.ha, plot, gse.lo)),
                    "snf.kg.ha"=with(nod, tapply(Nfixd.kg.ha, plot, mean)),#normal arithmetic mean
                    "snf.kg.ha.sd"=with(nod, tapply(Nfixd.kg.ha, plot, sd)),#normal arithmetic sd
                    "snf.kg.ha.se"=with(nod, tapply(Nfixd.kg.ha, plot, sefun)),#normal arithmetic se
                    "snf.rate"=with(nod, tapply(Nfix.rate, plot, mean)),#mean snf rate per nodule biomass
                    "snf.rate.se"=with(nod, tapply(Nfix.rate.se, plot, mean)))#standard error SNF (kg/ha) for each plot
sym.plt.fudge<-sym.plt
sym.plt.fudge$age<-ifelse(sym.plt.fudge$plot=="JE", yes=20, no=sym.plt.fudge$age)
#write.csv(sym.plt, file="Symbiotic N Fixation_Plot Summary Data.csv",row.names=F)
#write.csv(sym.plt.fudge, file="Symbiotic N Fixation_Plot Summary_Year Fudge.csv", row.names=F)

########## Creating a Dataframe of N inputs by Age #############################
nin<-data.frame("stand.age"=c(19,29,37,100),
                "anf"=with(asym, tapply(Nfixd.ha, stand.age, mean, na.rm=T)),
                "anf.se"=with(asym, tapply(Nfixd.ha, stand.age, sefun)),
                "snf"=with(nod, tapply(Nfixd.kg.ha, age, mean, na.rm=T)),
                "snf.se"=with(nod, tapply(Nfixd.kg.ha,age,sefun)),
                "snf.gmn"=with(nod, tapply(Nfixd.kg.ha, age, gmn)),
                "snf.se.hi"=with(nod, tapply(Nfixd.kg.ha, age, gse.hi)),
                "snf.se.lo"=with(nod, tapply(Nfixd.kg.ha, age, gse.lo)))
nin$snf.anf.ratio<-with(nin, snf/anf)
nin$anf.snf.ratio<-with(nin, anf/snf)
nin$prop.snf<-with(nin, snf/(snf+anf))
nin$snfg.anf.ratio<-with(nin, snf.gmn/anf)
nin$anf.snfg.ratio<-with(nin, anf/snf.gmn)
nin$prop.snfg<-with(nin, snf.gmn/(snf.gmn+anf))
#write.csv(nin, file="Nitrogen inputs by Age Table.csv", row.names=F)
#################################################################################
############# Soil Nitrogen Availability ########################################
#################################################################################
sn<-read.csv("Bosques Soil N Data_All.csv")
#Calculating Inorganic N concentrations in ug/g of soil 
sn$nh<-with(sn, (nh.con*.03)/soil.dry.mass)#uses the volume of KCl and soil dry weight to give NH4 in units of ug per g of soil
sn$no<-with(sn, (no.con*.03)/soil.dry.mass)#uses the volume of KCl and soil dry weight to give NO3 in units of ug per g of soil
sn$n.tot<-sn$nh+sn$no
sn$age<-ifelse(sn$plot%in%c("BEJ","JE"), yes=19, no=100)
sn$age<-ifelse(sn$plot=="LEPS", yes=37, no=sn$age)
sn$age<-ifelse(sn$plot=="LSUR", yes=29, no=sn$age)
sn$samp.schm<-c(rep("T",300),rep("S",150))
sn$prop.nh.tot<-with(sn, nh/n.tot)
sn$prop.no.tot<-with(sn, no/n.tot)
#write.csv(sn, file="Soil N data_Individual Cores_all data.csv", row.names=F)
#Converting concentrations (ug/g) to gN m-2 using bulk density

#Creating a dataframe of plot means for making figures
nit.plt<-data.frame("plot"=sort(unique(sn$plot)),
                    "age"=c(19,19,100,37,29),
                    "no"=with(sn, tapply(no, plot, mean, na.rm=T)),
                    "no.se"=with(sn, tapply(no, plot, sefun)),
                    "nh"=with(sn, tapply(nh, plot, mean,na.rm=T)),
                    "nh.se"=with(sn, tapply(nh, plot, sefun)),
                    "n.tot"=with(sn, tapply(n.tot, plot, mean,na.rm=T)),
                    "n.tot.se"=with(sn, tapply(n.tot, plot, sefun)),
                    "prop.no"=with(sn, tapply(prop.no.tot,plot,mean, na.rm=T)),
                    "prop.no.se"=with(sn, tapply(prop.no.tot,plot,sefun)))
nit.plt.fudge<-nit.plt
nit.plt.fudge$age<-ifelse(nit.plt.fudge$plot=="JE", yes=20, no=nit.plt.fudge$age)
#write.csv(nit.plt, file="Soil Nitrogen_Plot Summary Data.csv",row.names=F)
#write.csv(nit.plt.fudge, file="Soil Nitrogen_Plot Summary_Year Fudge.csv", row.names=F)

#################################################################################################
##### Combining Tree-based Soil N and SNF data ##################################################
#################################################################################################
tr<-read.csv("Tree-based Soil N Data_2017.csv")
tr$nh<-with(tr, (nh.con*.03)/dry.mass)#uses the volume of KCl and soil dry weight to give NH4 in units of ug per g of soil
tr$no<-with(tr, (no.con*.03)/dry.mass)#uses the volume of KCl and soil dry weight to give NO3 in units of ug per g of soil
tr$n.tot<-tr$nh+tr$no
nod$core.id.mrg<-as.factor(paste(nod$plot,"_T_",nod$core))
nod.mrg<-nod[,c("core.id.mrg","nod.num","nod.bio","NCI","age","nod.size.mg","nod.bio.m2","Nfix.rate","Nfix.rate.se","Nfixd.g.m2","Nfixd.kg.ha")]
tr<-merge(tr,nod.mrg,by="core.id.mrg",all.x=T,all.y=F)
adlt$core.id.mrg<-as.factor(paste(adlt$plot,"_T_",adlt$sample.id))
adlt.mrg<-adlt[,c("core.id.mrg","fixNCI","tr.DAP","tr.height","tr.NCI","tr.fixNCI","tr.nh","tr.RH", "subplot.ID","tr.Growth")]
tr<-merge(tr,adlt.mrg,by="core.id.mrg",all.x=T,all.y=T)
tr$hasnods<-as.factor(ifelse(tr$nod.bio>0, yes=1, no=0))
tr$logsnf<-log(tr$Nfixd.kg.ha+1)
#write.csv(tr, "Tree-based N Fixation and Soil N Data.csv", row.names=F)
#######################################################################################
##### Putting all our information together at the plot level ##########################
#######################################################################################
stemden<-read.csv("Plot-level stem densities.csv")
plt.alldat<-cbind(tr.plt,stemden[,c(2:4)],sym.plt[,c(10:17)],asym.plt[,c(3:10)],nit.plt[,c(3:8)],nod.plt[,c(3:10)])
plt.alldat$BAI.se<-with(tr.sbplt, tapply(BAI*100, plot, sefun))
plt.alldat$fixerBA.se<-with(tr.sbplt, tapply(fixBA*100, plot, sefun))
plt.alldat$fixerBAI.se<-with(tr.sbplt, tapply(fixBAI*100, plot, sefun))
plt.alldat$Ninputs.gmn.tot<-with(plt.alldat, snf.kg.ha.gmn+Nfixd.ha)
plt.alldat$Ninputs.gmn.sehi<-with(plt.alldat, Ninputs.gmn.tot+(snf.kg.ha.gse.hi-snf.kg.ha.gmn)+Nfixd.ha.se)
plt.alldat$Ninputs.gmn.selo<-with(plt.alldat, Ninputs.gmn.tot-(snf.kg.ha.gmn-snf.kg.ha.gse.lo)-Nfixd.ha.se)
plt.alldat$Ninputs.amn.tot<-with(plt.alldat, snf.kg.ha+Nfixd.ha)#total N inputs using normal arithmetic mean
plt.alldat$Ninputs.amn.totse<-with(plt.alldat, snf.kg.ha.se+Nfixd.ha.se)#SE of total N inputs using normal arithmetic mean
plt.alldat$Win.p<-(with(nod[nod$nod.bio>0,], tapply(nod.bio,plot,length))/80)#calculates p from Winbourne et al. 2018 (the probability of finding nodules in a given core)
plt.alldat$Win.v<-with(nod[nod$nod.bio>0,], (tapply(nod.bio,plot,sd))/(tapply(nod.bio,plot,mean)))#calculates v from Winbourne et al. 2018 (the variance in nodulation for cores that have nodules)
#write.csv(plt.alldat, file="Concatenated Plot Summary Data.csv",row.names=F)
##############################################################################################################
##### Calculating Bootstrapped Standard Errors for total N inputs at the plot level ##########################
##############################################################################################################
distdat<-data.frame("plot"=c(rep("BEJ",1000),rep("JE",1000),rep("LEPP",1000),rep("LEPS",1000),rep("LSUR",1000)))
m<-c(2.570269,1.896448,2.824337,8.652544,3.388757)
s<-c(1.423358,1.139941,1.351248,4.189098,2.097591)
loc<-log(m^2/sqrt((s^2+m^2)))
sh<-sqrt(log(1+(s^2/m^2)))
distdat$anf<-rlnorm(n=5000,meanlog = c(rep(loc[1],1000),rep(loc[2],1000),rep(loc[3],1000),rep(loc[4],1000),rep(loc[5],1000)),
                    sdlog=c(rep(sh[1],1000),rep(sh[2],1000),rep(sh[3],1000),rep(sh[4],1000),rep(sh[5],1000)))
sym.nz<-nod[nod$Nfixd.kg.ha>0,]
nm<-with(sym.nz, tapply(Nfixd.kg.ha,plot,mean))
ns<-with(sym.nz, tapply(Nfixd.kg.ha,plot,sd))
nloc<-log(nm^2/sqrt((ns^2+nm^2)))
nsh<-sqrt(log(1+(ns^2/nm^2)))
distdat$snf<-c(rep(0,800),rep(995,200),rep(0,850),rep(996,150),rep(0,888),rep(997,112),rep(0,925),rep(998,75),rep(0,950),rep(999,50))
distdat$snf<-ifelse(distdat$snf==995,yes=rlnorm(n=200,meanlog=rep(nloc[1],200),sdlog=rep(nsh[1],200)),no=distdat$snf)
distdat$snf<-ifelse(distdat$snf==996,yes=rlnorm(n=150,meanlog=rep(nloc[1],150),sdlog=rep(nsh[1],150)),no=distdat$snf)
distdat$snf<-ifelse(distdat$snf==997,yes=rlnorm(n=112,meanlog=rep(nloc[1],112),sdlog=rep(nsh[1],112)),no=distdat$snf)
distdat$snf<-ifelse(distdat$snf==998,yes=rlnorm(n=75,meanlog=rep(nloc[1],75),sdlog=rep(nsh[1],75)),no=distdat$snf)
distdat$snf<-ifelse(distdat$snf==999,yes=rlnorm(n=50,meanlog=rep(nloc[1],50),sdlog=rep(nsh[1],50)),no=distdat$snf)
distdat$nin<-distdat$anf+distdat$snf
distdat$lognin<-log(distdat$nin)

#Now we have all the dataframes set up. Each variable of interest has it's own dataframe
#Tree growth=tr.sbplt
#asym=ANF
#nod=nodule biomass (also includes the SNF data)
#sym=SNF 
#sn=Soil Nitrogen
#There is also a "tr" dataframe of combined variables that were sampled in a compatible way at the tree level.
#Finally, there is a "plt.alldat" dataframe with plot averages for all variables
plot.correlations<-rcorr(as.matrix(plt.alldat[,-1]))
#write.csv(plot.correlations$P, file="Plot-level correlation table.csv",row.names=T)


##################################################################################
##################################################################################
#### Creating a Table 1 of Plot-Level Metrics for Manuscript #####################
##################################################################################
##################################################################################

tab1<-data.frame("plot"=plt.alldat$plot,
                 "age"=c(19,19,"OG",37,29),
                 "BA"=plt.alldat$BA/10000,
                 "fixerBA"=plt.alldat$fixerBA/10000,
                 "soilN"=plt.alldat$n.tot,
                 "soilN.se"=plt.alldat$n.tot.se,
                 "anf"=plt.alldat$Nfixd.ha,
                 "anf.se"=plt.alldat$Nfixd.ha.se,
                 "snf.mn"=sym.plt$snf.kg.ha, #non lognormal arithmetic mean (don't report this)
                 "snf.se"=sym.plt$snf.kg.ha.se,
                 "snf.gmn"=sym.plt$snf.kg.ha.gmn,#0-inflated geometric mean
                 "snf.gse.hi"=sym.plt$snf.kg.ha.gse.hi,#0-inflated geometric mean plus se
                 "snf.gse.lo"=sym.plt$snf.kg.ha.gse.lo,#0-inflated geometric mean minus se
                 "p"=plt.alldat$Win.p,
                 "v"=plt.alldat$Win.v)
tab1<-tab1[c(1,2,5,4,3),]
#write.csv(tab1,file="Table 1_Plot summary stats.csv",row.names=F)
