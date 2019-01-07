setwd("C:/Users/Benton/Desktop/Work Files/Research/Costa Rica/Successional Nitrogen & Fixation Dynamics/Data Analysis and Code Files")
library(ggplot2);library(RColorBrewer);library(scales);library(grid);library(stringr)
source("C:\\RCode\\Homegrown Functions.R")
source("C:\\RCode\\fte_ggplot_theme.R")

#Reading in the dataframes
#Most of these are summary data put together in the "Successional Dynamics_Calculations and Statistics" R file
#The "Year Fudge" part of these filenames indicates that plot JE's age has been fudged to 20 for plotting purposes even though it's actually 19
tree<-read.csv("Tree Growth_Plot Summary_Year Fudge.csv")
asym<-read.csv("Asymbiotic Fixation_Plot Summary_Year Fudge.csv")
nod<-read.csv("Nodule Biomass_Plot Summary_Year Fudge.csv")
sym<-read.csv("Symbiotic N Fixation_Plot Summary_Year Fudge.csv")
sn<-read.csv("Soil Nitrogen_Plot Summary_Year Fudge.csv")
nod.tr<-read.csv("Tree-based N Fixation and Soil N Data.csv")
pltsm<-read.csv("Concatenated Plot Summary Data.csv")
ninage<-read.csv("Nitrogen inputs by Age Table.csv")
trsdf<-read.csv("Tree Growth_Subplot Summary_Year Fudge.csv")
trsnf<-read.csv("Tree-based N Fixation and Soil N Data.csv")

t1<-data.frame("plot"=sort(unique(trsdf$plot)),
               "stand.age"=c(19,20,100,37,29),
               "BA"=with(trsdf, tapply((BA*(100/10000)), plot, mean)),
               "BA.se"=with(trsdf, tapply((BA*(100/10000)), plot, sefun)),
               "fixerBA"=with(trsdf, tapply((fixBA*(100/10000)), plot, mean)),
               "fixerBA.se"=with(trsdf, tapply((fixBA*(100/10000)), plot, sefun)),
               "nonBA"=with(trsdf, tapply((nonBA*(100/10000)), plot, mean)),
               "nonBA.se"=with(trsdf, tapply((nonBA*(100/10000)), plot, sefun)),
               "BAI"=with(trsdf, tapply((BAI*(100/10000)), plot, mean)),
               "BAI.se"=with(trsdf, tapply((BAI*(100/10000)), plot, sefun)),
               "fixerBAI"=with(trsdf, tapply((fixBAI*(100/10000)), plot, mean)),
               "fixerBAI.se"=with(trsdf, tapply((fixBAI*(100/10000)), plot, sefun)),
               "nonBAI"=with(trsdf, tapply((nonBAI*(100/10000)), plot, mean)),
               "nonBAI.se"=with(trsdf, tapply((nonBAI*(100/10000)), plot, sefun)),
               "Fixer.Type"="Total")

t2<-data.frame("plot"=sort(unique(trsdf$plot)),
               "stand.age"=c(19,19.2,45,37,29),
               "BA"=with(trsdf, tapply((BA*(100/10000)), plot, mean)),
               "BA.se"=with(trsdf, tapply((BA*(100/10000)), plot, sefun)),
               "fixerBA"=with(trsdf, tapply((fixBA*(100/10000)), plot, mean)),
               "fixerBA.se"=with(trsdf, tapply((fixBA*(100/10000)), plot, sefun)),
               "nonBA"=with(trsdf, tapply((nonBA*(100/10000)), plot, mean)),
               "nonBA.se"=with(trsdf, tapply((nonBA*(100/10000)), plot, sefun)),
               "BAI"=with(trsdf, tapply((BAI*(100/10000)), plot, mean)),
               "BAI.se"=with(trsdf, tapply((BAI*(100/10000)), plot, sefun)),
               "fixerBAI"=with(trsdf, tapply((fixBAI*(100/10000)), plot, mean)),
               "fixerBAI.se"=with(trsdf, tapply((fixBAI*(100/10000)), plot, sefun)),
               "nonBAI"=with(trsdf, tapply((nonBAI*(100/10000)), plot, mean)),
               "nonBAI.se"=with(trsdf, tapply((nonBAI*(100/10000)), plot, sefun)),
               "Fixer.Type"="Total")

#### Figure 1 Version 2 ##############################################
#### Bar graphs of growth, N inputs, and N availability ##############
######################################################################
sn2<-sn
sn2$age<-c(19,19.2,45,37,29)
snx1<-seq(19,40,.1)
c<-15.38823
a<-58.04385
sig<-3.87257
u_nit<-26.38389
sny1<-log(c+(a-c)*exp((-(snx1-u_nit)^2)/(2*sig^2)))
sncrv1<-data.frame("x"=snx1,"y"=exp(sny1))
snx2<-seq(97,102,.1)
sny2<-log(c+(a-c)*exp((-(snx2-u_nit)^2)/(2*sig^2)))
sncrv2<-data.frame("x"=(snx2-54),"y"=exp(sny2))
#ylb<-strwrap(x=expression("Soil Inorganic N"~(mg~N~kg^{-1}~soil)),width=20,simplify = F)
sn.2.fig<-ggplot(sn2, aes(x=age, y=n.tot))+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(n.tot-n.tot.se), ymax=(n.tot+n.tot.se)),width=.3)+
  xlab("Forest Age (yr)")+
  ylab(expression("Inorganic N"~(mg~N~kg^{-1}~soil)))+
  #ylab(ylb)+
  scale_x_continuous(breaks=c(20,25,30,35,40,45),labels=c("20","25","30","35","40","OG"),)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),axis.title.y=element_text(size=15))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=17,xend=17,y=0,yend=60),colour="black")+
  geom_segment(aes(x=17,xend=41.5,y=0,yend=0),colour="black")+
  geom_segment(aes(x=43,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=41.25,xend=41.75,y=-3.33,yend=3.33),colour="black")+
  geom_segment(aes(x=42.75,xend=43.25,y=-3.33,yend=3.33),colour="black")+
  geom_line(aes(x=x,y=y), data=sncrv1,colour="black")+
  geom_line(aes(x=x,y=y), data=sncrv2,colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.25))


asym2<-asym
asym2$stand.age<-c(19,19.2,45,37,29)
anfx1<-seq(19,40,.1)
anfy1<-exp(.758+(2.104-.758)*exp((-(anfx1-36.096)^2)/(2*4.286^2)))
anfcrv1<-data.frame("x"=anfx1,"y"=anfy1)
anfx2<-seq(97,102,.1)
anfy2<-exp(.758+(2.104-.758)*exp((-(anfx2-36.096)^2)/(2*4.286^2)))
anfcrv2<-data.frame("x"=(anfx2-54),"y"=anfy2)
anf.2.fig<-ggplot(asym2, aes(x=stand.age, y=Nfixd.ha))+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(Nfixd.ha-Nfixd.ha.se), ymax=(Nfixd.ha+Nfixd.ha.se)),width=.3)+
  xlab("Forest Age (yr)")+
  ylab(expression("ANF"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_x_continuous(breaks=c(20,25,30,35,40,45),labels=c("20","25","30","35","40","OG"),)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=17,xend=17,y=0,yend=10),colour="black")+
  geom_segment(aes(x=17,xend=41.5,y=0,yend=0),colour="black")+
  geom_segment(aes(x=43,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=41.25,xend=41.75,y=-.666,yend=.666),colour="black")+
  geom_segment(aes(x=42.75,xend=43.25,y=-.666,yend=.666),colour="black")+
  geom_line(aes(x=x,y=y), data=anfcrv1,colour="black")+
  geom_line(aes(x=x,y=y), data=anfcrv2,colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.257))

sym2<-sym
sym2$age<-c(19,19.2,45,37,29)
snfx1<-seq(19,40,.1)
snfy1<-(1-.882493)*exp((3.640921-.01976*snfx1))
snfcrv1<-data.frame("x"=snfx1,"y"=snfy1)
snfx2<-seq(97,102,.1)
snfy2<-(1-.882493)*exp((3.640921-.01976*snfx2))
snfcrv2<-data.frame("x"=(snfx2-54),"y"=snfy2)
snf.2.fig<-ggplot(sym2, aes(x=age, y=snf.kg.ha.gmn))+
  geom_point(colour="grey50", size=5, shape=17)+
  geom_errorbar(aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi)),width=.3, colour="grey50")+
  geom_point(data=sym2, mapping=aes(x=age-.35, y=snf.kg.ha), size=5, colour="black")+
  geom_errorbar(data=sym2, mapping=aes(x=age-.35,ymin=(snf.kg.ha-snf.kg.ha.se),ymax=(snf.kg.ha+snf.kg.ha.se)),width=.3)+
  xlab("Forest Age (yr)")+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_x_continuous(breaks=c(20,25,30,35,40,45),labels=c("20","25","30","35","40","OG"),)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=17,xend=17,y=0,yend=20),colour="black")+
  geom_segment(aes(x=17,xend=41.5,y=0,yend=0),colour="black")+
  geom_segment(aes(x=43,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=41.25,xend=41.75,y=-1.3,yend=1.3),colour="black")+
  geom_segment(aes(x=42.75,xend=43.25,y=-1.3,yend=1.3),colour="black")+
  geom_line(aes(x=x,y=y), data=snfcrv1,colour="black")+
  geom_line(aes(x=x,y=y), data=snfcrv2,colour="black")+
  ggtitle("c)")+
  theme(plot.title=element_text(hjust=-.38))

t2.2<-t2
t2.2$age<-c(19,19.2,45,37,29)
t2x1<-seq(19,40,.1)
t2y1<-(t2x1*-.0088917)+1.3867114
t2crv1<-data.frame("x"=t2x1,"y"=t2y1)
t2x2<-seq(97,102,.1)
t2y2<-(t2x2*-.0088917)+1.3867114
t2crv2<-data.frame("x"=t2x2-54,"y"=t2y2)
tbai.2.fig<-ggplot(t2, aes(x=stand.age, y=BAI))+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(BAI-BAI.se), ymax=(BAI+BAI.se)),width=.3)+
  xlab("Forest Age (yr)")+
  ylab(expression("BAI"~(m^{2}~ha^{-1}~yr^{-1})))+
  scale_x_continuous(breaks=c(20,25,30,35,40,45),labels=c("20","25","30","35","40","OG"),)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=17,xend=17,y=0,yend=1.5),colour="black")+
  geom_segment(aes(x=17,xend=41.5,y=0,yend=0),colour="black")+
  geom_segment(aes(x=43,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=41.25,xend=41.75,y=-.1,yend=.1),colour="black")+
  geom_segment(aes(x=42.75,xend=43.25,y=-.1,yend=.1),colour="black")+
  geom_line(aes(x=x,y=y), data=t2crv1,colour="black")+
  geom_line(aes(x=x,y=y), data=t2crv2,colour="black")+
  ggtitle("d)")+
  theme(plot.title=element_text(hjust=-.255))

#png(filename = "Fig 1_Plot-Level Successional N Dynamics_V2.png", width=9, height=8, units="in", res=300)
multiplot(sn.2.fig,snf.2.fig,anf.2.fig,tbai.2.fig, cols=2)
#dev.off()

png(filename = "Fig 1.V2_Plot-Level Successional N Dynamics_V2.png", width=9, height=8, units="in", res=300)
multiplot(sn.2.fig,snf.2.fig,anf.2.fig,tbai.2.fig, cols=2)
dev.off()
#### Figure 2 ########################################################
#### Total N inputs and Growth #######################################
######################################################################
pltsm$Ninputs.gmn.tot2<-pltsm$Ninputs.gmn.tot^2
pltsm$Ninputs.amn.tot2<-pltsm$Ninputs.amn.tot^2
quad<-with(pltsm, lm(BAI/10000~Ninputs.gmn.tot+Ninputs.gmn.tot2))
amn.lm<-with(pltsm, lm(BAI/10000~Ninputs.amn.tot))
#x<-seq(3,12,.1)
quadfun<-function(x){coef(quad)[1]+(coef(quad)[2]*x)+(coef(quad)[3]*x^2)}
amn.lmfun<-function(x){coef(amn.lm)[1]+(coef(amn.lm)[2]*x)}
grw.ninputs_V1<-ggplot(pltsm, aes(x=Ninputs.gmn.tot, y=(BAI/10000)))+
  stat_function(fun=quadfun,xlim = c(3.2,11.5))+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=((BAI-BAI.se)/10000), ymax=((BAI+BAI.se)/10000),width=.3))+
  geom_errorbarh(aes(xmin=(Ninputs.gmn.selo),xmax=(Ninputs.gmn.sehi),height=.06))+
  xlab(expression("N Fixation"~(kg~N~ha^{-1}~yr^{-1})))+
  ylab(expression("BAI"~(m^{2}~ha^{-1}~yr^{-1})))+
  scale_x_continuous(breaks=c(4,6,8,10,12))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=2.5,xend=2.5,y=.2,yend=1.5),colour="black")+
  geom_segment(aes(x=2.5,xend=15,y=.2,yend=.2),colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.15))

grw.ninputs_V2<-ggplot(pltsm, aes(x=Ninputs.amn.tot, y=(BAI/10000)))+
  #stat_function(fun=amn.lmfun,xlim = c(3.2,15.5))+
  geom_smooth(method="lm", colour="black",se=F)+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=((BAI-BAI.se)/10000), ymax=((BAI+BAI.se)/10000),width=.3))+
  geom_errorbarh(aes(xmin=(Ninputs.amn.tot-Ninputs.amn.totse),xmax=(Ninputs.amn.tot+Ninputs.amn.totse),height=.06))+
  xlab(expression("N Fixation"~(kg~N~ha^{-1}~yr^{-1})))+
  ylab(expression("BAI"~(m^{2}~ha^{-1}~yr^{-1})))+
  scale_x_continuous(breaks=c(4,8,12,16,20))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=2.5,xend=2.5,y=.2,yend=1.5),colour="black")+
  geom_segment(aes(x=2.5,xend=20,y=.2,yend=.2),colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.15))

pltanf.bai<-ggplot(pltsm, aes(y=(BAI/10000), x=Nfixd.ha))+
  geom_point(colour="black", size=5)+
  geom_errorbarh(aes(xmin=(Nfixd.ha-Nfixd.ha.se), xmax=(Nfixd.ha+Nfixd.ha.se),height=.05))+
  geom_errorbar(aes(ymin=((BAI-BAI.se)/10000), ymax=((BAI+BAI.se)/10000),width=.3))+
  xlab(expression("ANF"~(kg~N~ha^{-1}~yr^{-1})))+
  ylab(expression("BAI"~(m^{2}~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=.35,yend=1.6),colour="black")+
  geom_segment(aes(x=0,xend=10,y=.35,yend=.35),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.15))

pltsnf.bai<-ggplot(pltsm, aes(y=(BAI/10000), x=snf.kg.ha.gmn))+
  geom_point(colour="black", size=5)+
  geom_errorbarh(aes(xmin=(snf.kg.ha.gse.lo), xmax=(snf.kg.ha.gse.hi),height=.05))+
  geom_errorbar(aes(ymin=((BAI-BAI.se)/10000), ymax=((BAI+BAI.se)/10000),width=.3))+
  xlab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  ylab(expression("BAI"~(m^{2}~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=.35,yend=1.6),colour="black")+
  geom_segment(aes(x=0,xend=10,y=.35,yend=.35),colour="black")+
  ggtitle("c)")+
  theme(plot.title=element_text(hjust=-.15))

pltsnf.bai.2<-ggplot(pltsm, aes(y=(BAI/10000), x=snf.kg.ha))+
  geom_point(colour="black", size=5)+
  geom_errorbarh(aes(xmin=(snf.kg.ha-snf.kg.ha.se), xmax=(snf.kg.ha+snf.kg.ha.se),height=.05))+
  geom_errorbar(aes(ymin=((BAI-BAI.se)/10000), ymax=((BAI+BAI.se)/10000),width=.3))+
  xlab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  ylab(expression("BAI"~(m^{2}~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=.35,yend=1.6),colour="black")+
  geom_segment(aes(x=0,xend=20,y=.35,yend=.35),colour="black")+
  ggtitle("c)")+
  theme(plot.title=element_text(hjust=-.15))

png(filename = "Fig 2_N inputs and Tree Growth.png", width=6, height=11, units="in", res=300)
multiplot(grw.ninputs_V1, pltanf.bai, pltsnf.bai, cols=1)
dev.off()

png(filename = "Fig 2.V2_N inputs and Tree Growth.png", width=6, height=11, units="in", res=300)
multiplot(grw.ninputs_V2, pltanf.bai, pltsnf.bai.2, cols=1)
dev.off()
#### Figure 3 ########################################################
#### Relative importance of ANF vs. SNF ##############################
######################################################################
ninage2<-data.frame("age"=rep(ninage$stand.age,2),
                    "n.input"=c(ninage$anf,ninage$snf.gmn),
                    "n.input.se.hi"=c((ninage$anf.se), (ninage$snf.se.hi-ninage$snf.gmn)),
                    "n.input.se.lo"=c((ninage$anf.se), (ninage$snf.gmn-ninage$snf.se.lo)),
                    "input.type"=c(rep("ANF",4),rep("SNF",4)))
ninage2$bartop<-c((ninage2[c(1:4),2]+ninage2[c(5:8),2]),ninage2[c(5:8),2])

ninage3<-data.frame("age"=rep(ninage$stand.age,2),
                    "n.input"=c(ninage$anf,ninage$snf),
                    "n.input.se"=c((ninage$anf.se), (ninage$snf.se)),
                    #"n.input.se.lo"=c((ninage$anf.se), (ninage$snf.gmn-ninage$snf.se.lo)),
                    "input.type"=c(rep("ANF",4),rep("SNF",4)))
ninage3$bartop<-c((ninage3[c(1:4),2]+ninage3[c(5:8),2]),ninage3[c(5:8),2])

input.plt<-ggplot(data=ninage2,aes(x=as.factor(age),y=n.input,fill=input.type))+
  geom_bar(stat="identity", colour="black")+
  geom_errorbar(aes(ymax=(bartop+n.input.se.hi),ymin=(bartop-n.input.se.lo)),width=0.5,position = position_dodge(width=0))+
  scale_fill_manual(values=c("gray88","gray35"))+
  scale_y_continuous(breaks=c(2.5,5,7.5,10))+
  xlab("Forest Age (yr)")+
  ylab(expression("N Fixation"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_x_discrete(breaks=c(19,29,37,100),labels=c("19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="top",legend.title=element_blank())+
  geom_segment(aes(x=0,xend=0,y=0,yend=12),colour="black")+
  geom_segment(aes(x=0,xend=5,y=0,yend=0),colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.175))

input.plt2<-ggplot(data=ninage3,aes(x=as.factor(age),y=n.input,fill=input.type))+
  geom_bar(stat="identity", colour="black")+
  geom_errorbar(aes(ymax=(bartop+n.input.se),ymin=(bartop-n.input.se)),width=0.5,position = position_dodge(width=.1))+
  scale_fill_manual(values=c("gray88","gray35"))+
  #scale_y_continuous(breaks=c(2.5,5,7.5,10))+
  xlab("Forest Age (yr)")+
  ylab(expression("N Fixation"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_x_discrete(breaks=c(19,29,37,100),labels=c("19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="top",legend.title=element_blank())+
  geom_segment(aes(x=0,xend=0,y=0,yend=15),colour="black")+
  geom_segment(aes(x=0,xend=5,y=0,yend=0),colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.175))

propinput.plt<-ggplot(data=ninage,aes(x=as.factor(stand.age),y=(prop.snfg*100)))+
  geom_bar(stat="identity", colour="black")+
  scale_fill_manual("gray35")+
  xlab("Forest Age (yr)")+
  ylab("SNF Contribution to N Fixation (%)")+
  scale_x_discrete(breaks=c(19,29,37,100),labels=c("19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),axis.title.y=element_text(size=14))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="top",legend.title=element_blank())+
  geom_segment(aes(x=0,xend=5,y=50,yend=50),colour="black", linetype="dashed")+
  geom_segment(aes(x=0,xend=0,y=0,yend=70),colour="black")+
  geom_segment(aes(x=0,xend=5,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.09))

propinput.plt2<-ggplot(data=ninage,aes(x=as.factor(stand.age),y=(prop.snf*100)))+
  geom_bar(stat="identity", colour="black")+
  scale_fill_manual("gray35")+
  xlab("Forest Age (yr)")+
  ylab("SNF Contribution to N Fixation (%)")+
  scale_x_discrete(breaks=c(19,29,37,100),labels=c("19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),axis.title.y=element_text(size=14))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="top",legend.title=element_blank())+
  geom_segment(aes(x=0,xend=5,y=50,yend=50),colour="black", linetype="dashed")+
  geom_segment(aes(x=0,xend=0,y=0,yend=90),colour="black")+
  geom_segment(aes(x=0,xend=5,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.09))

png(filename = "Fig 3_Components of N Inputs.png", width=6, height=8, units="in", res=300)
multiplot(input.plt,propinput.plt, cols=1)
dev.off()

png(filename = "Fig 3.V2_Components of N Inputs.png", width=6, height=8, units="in", res=300)
multiplot(input.plt2,propinput.plt2, cols=1)
dev.off()
#### Figure 4 ########################################################
#### N-fixer abundance and SNF #######################################
######################################################################

with(trsnf, summary(lm(logsnf~log(tr.fixNCI+1))))
with(trsnf, lm(Nfixd.kg.ha~tr.fixNCI))
fixncisnf<-ggplot(trsnf, aes(x=log(tr.fixNCI+1), y=logsnf))+
  #geom_smooth(method="lm",colour="black",linetype=2,se=F)+
  geom_point(colour="black", size=3.5)+
  xlab("Crowding from N Fixers")+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_x_continuous(breaks=c(0.000999,2.398,3.932,6.2166,8.517,9.903),labels=c(".001","10","50","500","5,000","20,000"))+
  theme(axis.text.x=element_text(angle=30))+
  scale_y_continuous(breaks=c(0.000999,2.398,3.932,6.2166),labels=c("0","10","50","500"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=6),colour="black")+
  geom_segment(aes(x=0,xend=11,y=0,yend=0),colour="black")+
  geom_abline(mapping=NULL, data=NULL, slope=1, linetype="dotted")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.2))+
  annotate("text", label="1:1", x=5.5, y=6.2, size=5, color="black")

trsnf.grwth<-ggplot(nod.tr, aes(x=(tr.Growth), y=logsnf))+
  geom_point(colour="black", size=3)+
  xlab(expression("Tree Growth"~(cm^{2}~yr^{-1})))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_y_continuous(breaks=c(0.000999,.6931,2.398,3.932,5.9939),labels=c("0","1","10","50","400"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=6),colour="black")+
  geom_segment(aes(x=0,xend=.2,y=0,yend=0),colour="black")+
  ggtitle("c)")+
  theme(plot.title=element_text(hjust=-.15))

pltsm.lm<-with(pltsm, lm(snf.kg.ha.gmn~fixerBA))
fixsnf<-ggplot(pltsm, aes(x=(fixerBA/10000), y=snf.kg.ha.gmn))+
  #geom_smooth(method="lm",colour="black", linetype="dashed",se=F)+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi)),width=.3)+
  geom_errorbarh(aes(xmin=((fixerBA/10000)-(fixerBA.se/10000)), xmax=((fixerBA/10000)+(fixerBA.se/10000))),height=.3)+
  xlab(expression("Nitrogen Fixer Basal Area"~(m^2~ha^{-1})))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=4,xend=4,y=0,yend=10),colour="black")+
  geom_segment(aes(x=4,xend=12,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.15))

pltsm.lm2<-with(pltsm, lm(snf.kg.ha~fixerBA))
fixsnf2<-ggplot(pltsm, aes(x=(fixerBA/10000), y=snf.kg.ha))+
  geom_smooth(method="lm",colour="black", se=F)+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(snf.kg.ha-snf.kg.ha.se), ymax=(snf.kg.ha+snf.kg.ha.se)),width=.3)+
  geom_errorbarh(aes(xmin=((fixerBA/10000)-(fixerBA.se/10000)), xmax=((fixerBA/10000)+(fixerBA.se/10000))),height=.3)+
  xlab(expression("Nitrogen Fixer Basal Area"~(m^2~ha^{-1})))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=4,xend=4,y=0,yend=20),colour="black")+
  geom_segment(aes(x=4,xend=12,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.15))

fixsnf3<-ggplot(pltsm, aes(x=(fixerBA/10000), y=snf.kg.ha))+
  geom_smooth(method="lm",colour="grey50", se=F)+
  geom_point(colour="grey50", size=5, shape=17)+
  geom_errorbar(aes(ymin=(snf.kg.ha-snf.kg.ha.se), ymax=(snf.kg.ha+snf.kg.ha.se)),width=.3)+
  geom_errorbarh(aes(xmin=((fixerBA/10000)-(fixerBA.se/10000)), xmax=((fixerBA/10000)+(fixerBA.se/10000))),height=.3)+
  geom_point(data=pltsm, mapping=aes(x=(fixerBA/10000),y=snf.kg.ha.gmn),colour="black",size=5)+
  geom_errorbar(data=pltsm, aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi)),width=.3)+
  geom_errorbarh(data=pltsm, aes(xmin=((fixerBA/10000)-(fixerBA.se/10000)), xmax=((fixerBA/10000)+(fixerBA.se/10000))),height=.3)+
  xlab(expression("Nitrogen Fixer Basal Area"~(m^2~ha^{-1})))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=4,xend=4,y=0,yend=20),colour="black")+
  geom_segment(aes(x=4,xend=12,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.15))

fixgrwsnf<-ggplot(pltsm, aes(x=(fixerBAI/10000), y=snf.kg.ha.gmn))+
  #geom_smooth(method="lm",colour="black", linetype="dashed",se=F)+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi)),width=.025)+
  geom_errorbarh(aes(xmin=((fixerBAI/10000)-(fixerBAI.se/10000)), xmax=((fixerBAI/10000)+(fixerBAI.se/10000))),height=.3)+
  xlab(expression("Nitrogen Fixer BAI"~(m^2~ha^{-1}~yr^{-1})))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0.1,xend=0.1,y=0,yend=10),colour="black")+
  geom_segment(aes(x=0.1,xend=.6,y=0,yend=0),colour="black")+
  ggtitle("d)")+
  theme(plot.title=element_text(hjust=-.2))

fixgrwsnf2<-ggplot(pltsm, aes(x=(fixerBAI/10000), y=snf.kg.ha))+
  #geom_smooth(method="lm",colour="black", linetype="dashed",se=F)+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(snf.kg.ha-snf.kg.ha.se), ymax=(snf.kg.ha+snf.kg.ha.se)),width=.015)+
  geom_errorbarh(aes(xmin=((fixerBAI/10000)-(fixerBAI.se/10000)), xmax=((fixerBAI/10000)+(fixerBAI.se/10000))),height=.8)+
  xlab(expression("Nitrogen Fixer BAI"~(m^2~ha^{-1}~yr^{-1})))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0.1,xend=0.1,y=0,yend=20),colour="black")+
  geom_segment(aes(x=0.1,xend=.6,y=0,yend=0),colour="black")+
  ggtitle("d)")+
  theme(plot.title=element_text(hjust=-.2))

fixgrwsnf3<-ggplot(pltsm, aes(x=(fixerBAI/10000), y=snf.kg.ha))+
  #geom_smooth(method="lm",colour="black", linetype="dashed",se=F)+
  geom_point(size=5, shape=17, colour="grey50")+
  geom_errorbar(aes(ymin=(snf.kg.ha-snf.kg.ha.se), ymax=(snf.kg.ha+snf.kg.ha.se)),width=.015,colour="grey50")+
  geom_errorbarh(aes(xmin=((fixerBAI/10000)-(fixerBAI.se/10000)), xmax=((fixerBAI/10000)+(fixerBAI.se/10000))),height=.8,colour="grey50")+
  geom_point(data=pltsm, mapping=aes((x=fixerBAI/10000), y=snf.kg.ha.gmn), size=5)+
  geom_errorbar(data=pltsm, aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi)),width=.015)+
  geom_errorbarh(data=pltsm, aes(xmin=((fixerBAI/10000)-(fixerBAI.se/10000)), xmax=((fixerBAI/10000)+(fixerBAI.se/10000))),height=.8)+
  xlab(expression("Nitrogen Fixer BAI"~(m^2~ha^{-1}~yr^{-1})))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0.1,xend=0.1,y=0,yend=20),colour="black")+
  geom_segment(aes(x=0.1,xend=.6,y=0,yend=0),colour="black")+
  ggtitle("d)")+
  theme(plot.title=element_text(hjust=-.2))

png(filename = "Fig 4_SNF vs. Fixer Abundance and Growth.png", width=11, height=10, units="in", res=300)
multiplot(fixncisnf,trsnf.grwth,fixsnf,fixgrwsnf, cols=2)
dev.off()

png(filename = "Fig 4.V2_SNF vs. Fixer Abundance and Growth.png", width=11, height=10, units="in", res=300)
multiplot(fixncisnf,trsnf.grwth,fixsnf2,fixgrwsnf2, cols=2)
dev.off()

png(filename = "Fig 4.V3_SNF vs. Fixer Abundance and Growth.png", width=11, height=10, units="in", res=300)
multiplot(fixncisnf,trsnf.grwth,fixsnf3,fixgrwsnf3, cols=2)
dev.off()
#### Supplemental Figure 1 ##########################################
#### Example Plot Map (BEJ) ##########################################
######################################################################
bejtrs<-read.csv("Bosques Tree Data_All Trees File_11-16.csv")
bejtrs<-bejtrs[bejtrs$Locality_Code=="BEJ"&bejtrs$Year==2014,]
bejtrs$fixer<-as.factor(bejtrs$fixer)
bejnit<-read.csv("Tree-based N Fixation and Soil N Data.csv")
bejnit<-bejnit[bejnit$plot=="BEJ",]
bejplotnit<-read.csv("Bosques Soil N Data_All.csv")
bejplotnit<-bejplotnit[c(301:330),]
bejasym<-read.csv("Bosques Asymbiotic Fixation Data_2016.csv")
bejasym<-bejasym[bejasym$Plot=="BEJ",]
bejpltsnf<-read.csv("Bosques Nodule Biomass_Data_2016.csv")
bejpltsnf<-bejpltsnf[bejpltsnf$Plot=="BEJ",]

#png(filename = "Fig S1_Example Plot Map.png", width=8, height=14, units="in", res=300)
par(mar=c(5,5,2,2))
plot(1, type="n", axes=F, xlab="Short Edge of Plot (m)", ylab="Long edge of Plot (m)", xlim=c(0,80), ylim=c(0,200), bty="n",cex.lab=2.25, cex.axis=2.25)
axis(1, at=c(0,10,20,30,40,50), labels=c("0","10","20","30","40","50"))
axis(2)
points(bejtrs$X, bejtrs$Y, cex=(bejtrs$logDAP-1), bg=c(addTrans("blue",255),addTrans("red",255))[bejtrs$fixer], 
     xlab="Short Edge of Plot (m)", ylab="Long edge of Plot (m)", pch=21, bty="n")
rect(0,0,50,200, col=addTrans("white",210))
abline(v = c(0,10,20,30,40,50),lty=c(1,2,2,2,2,1))
abline(h = c(seq(0,200,10)),lty=c(1,rep(2,19),1))
rect(-5,-10,0,210, col="white", border=NA)
rect(50.05,-10,90,210, col="white", border=NA)
rect(-1,200.1,51,210, col="white", border=NA)
rect(-1,-10,51,-.1, col="white", border=NA)
points(bejnit$core.x, bejnit$core.y, pch=22, cex=2.5, bg="forestgreen")
points(bejplotnit$core.x, bejplotnit$core.y, pch=23,cex=2.5, bg="yellow")
points(bejpltsnf$X, bejpltsnf$Y, pch=24,cex=2.5,bg="orange")
points(bejasym$X, bejasym$Y, pch=25,cex=2.75, bg="magenta")
legend(50,150,
       c("Non-fixing Tree","Fixing Tree","Tree SNF & Soil N","Plot Soil N","Plot SNF","ANF"),
       pch=c(21,21,22,23,24,25),
       pt.bg=c("blue","red","forestgreen","yellow","orange","magenta"),bty="n",cex=1.5)
#dev.off()

#### Supplemental Figure 2 ##########################################
#### Breakdown of Soil N #############################################
######################################################################
snage1<-data.frame("age"=rep(sn$age,2),
                    "n"=c(sn$nh,sn$no),
                   "n.se"=c(sn$nh.se,sn$no.se),
                    "n.type"=c(rep("NH4",5),rep("NO3",5)))
snage2<-data.frame("age"=rep(sn$age,2),
                   "n"=c(sn$no,sn$nh),
                   "n.se"=c(sn$no.se,sn$nh.se),
                   "n.type"=c(rep("Nitrate",5),rep("Ammonium",5)))
snage2$bartop<-c(snage2[c(1:5),2],(snage2[c(1:5),2]+snage2[c(6:10),2]))
sntype.plt<-ggplot(data=snage2,aes(x=as.factor(age),y=n,fill=n.type))+
  geom_bar(stat="identity", colour="black")+
  geom_errorbar(aes(ymax=(snage2$bartop+snage2$n.se),ymin=(snage2$bartop-snage2$n.se)),width=0.5,position = position_dodge(width=.1))+
  scale_fill_manual(values=c("gray","gray40"))+
  xlab("Forest Age (yr)")+
  ylab(expression("Soil N"~(mg~N~kg^{-1}~soil)))+
  scale_x_discrete(breaks=c(19,20,29,37,100),labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), 
        legend.position="top",legend.title=element_blank())+
  geom_segment(aes(x=0,xend=0,y=0,yend=50),colour="black")+
  geom_segment(aes(x=0,xend=5,y=0,yend=0),colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.15))

prop.no<-sn
prop.no$age<-c(19,19.2,45,37,29)
prop.no$prop.no<-prop.no$prop.no*100
prop.no$prop.no.se<-prop.no$prop.no.se*100
prop.no.plt<-ggplot(data=prop.no,aes(x=as.factor(age),y=prop.no))+
  geom_bar(stat="identity", colour="black")+
  geom_errorbar(aes(ymin=(prop.no-prop.no.se), ymax=(prop.no+prop.no.se)),width=.3)+
  xlab("Forest Age (yr)")+
  ylab("Soil N Nitrate (%)")+
  scale_fill_manual(values=c("gray40"))+
  scale_x_discrete(breaks=c(19,19.2,29,37,45),labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="top",legend.title=element_blank())+
  geom_segment(aes(x=0,xend=0,y=0,yend=12),colour="black")+
  geom_segment(aes(x=0,xend=5,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.15))

#png(filename = "Fig S2_Components of soil N.png", width=5, height=9, units="in", res=300)
multiplot(sntype.plt,prop.no.plt, cols=1)
#dev.off()

#### Supplemental Figure 3 ##########################################
#### Breakdown of ANF ###############################################
######################################################################
area<-(3.14*(.04)^2)
asym$Lit.mass.m2<-asym$Lit.mass/area
asym$Lit.mass.m2.se<-asym$Lit.mass.se/area
lit.fig<-ggplot(asym, aes(x=as.factor(stand.age), y=Lit.mass.m2))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(Lit.mass.m2-Lit.mass.m2.se), ymax=(Lit.mass.m2+Lit.mass.m2.se)),width=.1)+
  xlab("Forest Age (yr)")+
  ylab(expression("Litter Mass"~(g~m^{-2})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  theme(plot.margin=unit(c(.2,.5,.2,1),"cm"))+
  geom_segment(aes(x=0,xend=0,y=0,yend=600),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  #annotate("text", x=c(1,2,3,4,5), y=42, label=c("b","ab","a","b","ab"), col=1,size=6)+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.42))

anf.rt.fig<-ggplot(asym, aes(x=as.factor(stand.age), y=(fix.rt*1000)))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=((fix.rt*1000)-(fix.rt.se*1000)), ymax=((fix.rt*1000)+(fix.rt.se*1000))),width=.1)+
  xlab("Forest Age (yr)")+
  ylab(expression("ANF Activity"~(g~N~kg^{-1}~yr^{-1})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  theme(plot.margin=unit(c(.2,.5,.2,1),"cm"))+
  geom_segment(aes(x=0,xend=0,y=0,yend=.25),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  #annotate("text", x=c(1,2,3,4,5), y=42, label=c("b","ab","a","b","ab"), col=1,size=6)+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.42))

#png(filename = "Fig S3_ANF Breakdown Figure.png", width=5, height=8, units="in", res=300)
par(mar=c(5,5.5,2,2))
multiplot(lit.fig, anf.rt.fig, cols=1)
#dev.off()

#### Supplemental Figure 4 ##########################################
#### Breakdown of SNF  ##############################################
######################################################################

nod.fig<-ggplot(nod, aes(x=as.factor(age), y=nod.bio.gmn))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(nod.bio.gse.lo), ymax=(nod.bio.gse.hi)),width=.1)+
  xlab("Forest Age (yr)")+
  ylab(expression("Nodulation"~(g~m^{-2})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=2.5),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  #annotate("text", x=c(1,2,3,4,5), y=8, label=c("a","a","a","a","a"), col=1,size=6)+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.1))

snf.rt.fig<-ggplot(sym, aes(x=as.factor(age), y=snf.rate))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(snf.rate-snf.rate.se), ymax=(snf.rate+snf.rate.se)),width=.1)+
  xlab("Forest Age (yr)")+
  ylab(expression("SNF Activity"~(g~N~g^{-1}~yr^{-1})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=.6),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  #annotate("text", x=c(1,2,3,4,5), y=42, label=c("b","ab","a","b","ab"), col=1,size=6)+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.1))

#png(filename = "Fig S4_SNF Breakdown Figure.png", width=5, height=8, units="in", res=300)
multiplot(nod.fig, snf.rt.fig, cols=1)
#dev.off()

#### Supplemental Figure 5 ##########################################
#### Breakdown of Nodulation  ########################################
######################################################################
area<-(3.14*(.04)^2)
core.ar.tot<-area*80 #Gives the total area sampled per plot (core area x number of cores per plot (80))
nod.num.fig<-ggplot(nod, aes(x=as.factor(age), y=(nod.num/core.ar.tot)))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  #geom_errorbar(aes(ymin=(nod.bio-nod.bio.se), ymax=(nod.bio+nod.bio.se)),width=.1)+
  xlab("Forest Age (yr)")+
  ylab(expression("Nodule Abundance"~(nodules~m^{-2})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=200),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  #annotate("text", x=c(1,2,3,4,5), y=8, label=c("a","a","a","a","a"), col=1,size=6)+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.1))

nod.size.fig<-ggplot(nod, aes(x=as.factor(age), y=nod.size))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(nod.size-nod.size.se), ymax=(nod.size+nod.size.se)),width=.1)+
  xlab("Forest Age (yr)")+
  ylab(expression("Nodule size"~(mg~nodule^{-1})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=110),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  #annotate("text", x=c(1,2,3,4,5), y=8, label=c("a","a","a","a","a"), col=1,size=6)+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.1))

#png(filename = "Fig S5_Nodulation Breakdown Figure.png", width=6, height=9, units="in", res=300)
multiplot(nod.num.fig, nod.size.fig, cols=1)
#dev.off()

#### Supplemental Figure 6 ########################################################
#### Relationship between SNF and soil N ########################################
######################################################################

crsnf.no<-ggplot(nod.tr, aes(x=(no), y=logsnf))+
  #geom_smooth(method="lm",colour="black")+
  geom_point(colour="black", size=3)+
  #geom_errorbar(aes(ymin=(snf.kg.ha-snf.kg.ha.se), ymax=(snf.kg.ha+snf.kg.ha.se)),width=.3)+
  xlab(expression("Soil NO"[3]^"-"*~(mg~kg^{-1}~"of NO"[3]^"-"*"-N")))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_y_continuous(breaks=c(0.000999,2.398,3.932,5.9939),labels=c(".001","10","50","400"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=6),colour="black")+
  geom_segment(aes(x=0,xend=12,y=0,yend=0),colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.15))

pltsnf.no<-ggplot(pltsm, aes(x=no, y=snf.kg.ha.gmn))+
  #geom_smooth(method="lm",colour="black")+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi),width=.1))+
  geom_errorbarh(aes(xmin=(no-no.se),xmax=(no+no.se),width=.3))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  xlab(expression("Soil NO"[3]^"-"*~(mg~kg^{-1}~"of NO"[3]^"-"*"-N")))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=.5,xend=.5,y=0,yend=14),colour="black")+
  geom_segment(aes(x=.5,xend=3,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.15))

nvec<-seq(5,92,.1)
svec<-(1-.902)*exp(2.113+.0361*nvec)
crv1<-data.frame("x"=nvec,"y"=log(svec))
crsnf.nh<-ggplot(nod.tr, aes(x=(nh), y=logsnf))+
  #geom_smooth(method="lm",colour="black")+
  geom_point(colour="black", size=3)+
  #geom_errorbar(aes(ymin=(snf.kg.ha-snf.kg.ha.se), ymax=(snf.kg.ha+snf.kg.ha.se)),width=.3)+
  xlab(expression("Soil NH"[4]^"+"*~(mg~kg^{-1}~"of NH"[4]^"+"*"-N")))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_y_continuous(breaks=c(0.000999,2.398,3.932,5.9939),labels=c(".001","10","50","400"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_line(aes(x=x,y=y), data=crv1,colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=6),colour="black")+
  geom_segment(aes(x=0,xend=100,y=0,yend=0),colour="black")+
  ggtitle("c)")+
  theme(plot.title=element_text(hjust=-.15))

pltsnf.nh<-ggplot(pltsm, aes(x=nh, y=snf.kg.ha.gmn))+
  #geom_smooth(method="lm",colour="black")+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi),width=.51))+
  geom_errorbarh(aes(xmin=(nh-nh.se),xmax=(nh+nh.se),width=.3))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  xlab(expression("Soil NH"[4]^"+"*~(mg~kg^{-1}~"of NH"[4]^"+"*"-N")))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=20,xend=20,y=0,yend=14),colour="black")+
  geom_segment(aes(x=20,xend=50,y=0,yend=0),colour="black")+
  ggtitle("d)")+
  theme(plot.title=element_text(hjust=-.15))

nvec<-seq(5,92,.1)
svec<-(1-.902)*exp(2.048+.0366*nvec)
crv1<-data.frame("x"=nvec,"y"=log(svec))
crsnf.ntot<-ggplot(nod.tr, aes(x=(n.tot), y=logsnf))+
  #geom_smooth(method="lm",colour="black")+
  geom_point(colour="black", size=3)+
  #geom_errorbar(aes(ymin=(snf.kg.ha-snf.kg.ha.se), ymax=(snf.kg.ha+snf.kg.ha.se)),width=.3)+
  xlab(expression("Inorganic N (mg"~kg^{-1}~"of NO"[3]^"-"*~"& NH"[4]^"+"*"-N)"))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  scale_y_continuous(breaks=c(0.000999,2.398,3.932,5.9939),labels=c(".001","10","50","400"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_line(aes(x=x,y=y), data=crv1,colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=6),colour="black")+
  geom_segment(aes(x=0,xend=100,y=0,yend=0),colour="black")+
  ggtitle("e)")+
  theme(plot.title=element_text(hjust=-.15))

pltsnf.ntot<-ggplot(pltsm, aes(x=n.tot, y=snf.kg.ha.gmn))+
  #geom_smooth(method="lm",colour="black")+
  geom_point(colour="black", size=5)+
  geom_errorbar(aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi),width=.51))+
  geom_errorbarh(aes(xmin=(n.tot-n.tot.se),xmax=(n.tot+n.tot.se),width=.3))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  xlab(expression("Inorganic N (mg"~kg^{-1}~"of NO"[3]^"-"*~"& NH"[4]^"+"*"-N)"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=20,xend=20,y=0,yend=14),colour="black")+
  geom_segment(aes(x=20,xend=50,y=0,yend=0),colour="black")+
  ggtitle("f)")+
  theme(plot.title=element_text(hjust=-.15))

#png(filename = "Fig S5_SNF vs Soil N.png", width=12, height=11, units="in", res=300)
multiplot(crsnf.no,crsnf.nh,crsnf.ntot,pltsnf.no,pltsnf.nh,pltsnf.ntot, cols=2)
#dev.off()

#### Supplemental Figure 7 ########################################################
#### Breakdown of Tree Growth ########################################
######################################################################
t3<-data.frame("plot"=sort(unique(trsdf$plot)),
               "stand.age"=c(19,20,100,37,29),
               "BA"=with(trsdf, tapply((fixBA*(100/10000)), plot, mean)),
               "BA.se"=with(trsdf, tapply((fixBA*(100/10000)), plot, sefun)),
               "BAI"=with(trsdf, tapply((fixBAI*(100/10000)), plot, mean)),
               "BAI.se"=with(trsdf, tapply((fixBAI*(100/10000)), plot, sefun)),
               "Fixer.Type"="Fixer")
t4<-data.frame("plot"=sort(unique(trsdf$plot)),
               "stand.age"=c(19,20,100,37,29),
               "BA"=with(trsdf, tapply((nonBA*(100/10000)), plot, mean)),
               "BA.se"=with(trsdf, tapply((nonBA*(100/10000)), plot, sefun)),
               "BAI"=with(trsdf, tapply((nonBAI*(100/10000)), plot, mean)),
               "BAI.se"=with(trsdf, tapply((nonBAI*(100/10000)), plot, sefun)),
               "Fixer.Type"="Non-Fixer")
t5<-t1[,c(1,2,3,4,9,10,15)]
trbadf<-rbind(t3,t4)
trbadf$bartop.ba<-c((trbadf[c(1:5),3]+trbadf[c(6:10),3]),trbadf[c(6:10),3])
ft.ba.fig<-ggplot(trbadf, aes(x=as.factor(stand.age), y=BA, fill=Fixer.Type))+
  geom_bar(stat="identity", colour="black")+
  geom_errorbar(aes(ymin=trbadf$bartop.ba-trbadf$BA.se, ymax=trbadf$bartop.ba+trbadf$BA.se),width=.1)+
  scale_fill_manual(values=c("gray","gray40"))+
  theme(text=element_text(size=20),axis.text.x=element_text(size=12, colour="black"),
        axis.text.y=element_text(size=15, colour="black"))+
  xlab("Forest Age (yr)")+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  ylab(expression("BA"~(m^{2}~ha^{-1})))+
  #scale_fill_manual(values=c("red","blue","darkgray"), labels=c("Fixer","Non-Fixer","Total"), name="Fixer Type")+
  theme(panel.background=element_rect(fill="white", color="white"),
        legend.position="top", legend.title = element_blank())+
  geom_segment(aes(x=0,xend=0,y=0,yend=35),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.1))

trbadf$bartop.bai<-c((trbadf[c(1:5),5]+trbadf[c(6:10),5]),trbadf[c(6:10),5])
ft.bai.fig<-ggplot(trbadf, aes(x=as.factor(stand.age), y=BAI, fill=Fixer.Type))+
  geom_bar(stat="identity", colour="black")+
  geom_errorbar(aes(ymin=trbadf$bartop.bai-trbadf$BAI.se, ymax=trbadf$bartop.bai+trbadf$BAI.se),width=.1)+
  scale_fill_manual(values=c("gray","gray40"))+
  theme(text=element_text(size=20),axis.text.x=element_text(size=12, colour="black"),
        axis.text.y=element_text(size=15, colour="black"))+
  xlab("Forest Age (yr)")+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  ylab(expression("BAI"~(m^{2}~ha^{-1}~yr^{-1})))+
  #scale_fill_manual(values=c("red","blue","darkgray"), labels=c("Fixer","Non-Fixer","Total"), name="Fixer Type")+
  theme(panel.background=element_rect(fill="white", color="white"),
        legend.position="none", legend.title = element_blank())+
  geom_segment(aes(x=0,xend=0,y=0,yend=1.5),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.1))

#png(filename = "Tree BA and BAI by Fixer Type.png", width=5, height=8, units="in", res=300)
multiplot(ft.ba.fig,ft.bai.fig, cols=1)
#dev.off()

#### Supplemental Figure 8 ########################################################
#### Relationship between SNF and Tree Growth ########################################
######################################################################
#trsnf.grwth<-ggplot(nod.tr, aes(x=(tr.Growth), y=logsnf))+
#  geom_point(colour="black", size=3)+
#  xlab(expression("Tree Growth"~(cm^{2}~yr^{-1})))+
#  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
#  scale_y_continuous(breaks=c(0.000999,2.398,3.932,5.9939),labels=c(".001","10","50","400"))+
#  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
#        axis.text.y=element_text(size=20, colour="black"))+
#  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
#  geom_segment(aes(x=0,xend=0,y=0,yend=6),colour="black")+
#  geom_segment(aes(x=0,xend=.2,y=0,yend=0),colour="black")+
#  ggtitle("b)")+
#  theme(plot.title=element_text(hjust=-.15))
fixncisnf.ln<-ggplot(trsnf, aes(x=tr.fixNCI, y=Nfixd.kg.ha))+
  geom_point(colour="black", size=3.5)+
  xlab("Crowding from N Fixers")+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=15, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=400),colour="black")+
  geom_segment(aes(x=0,xend=25000,y=0,yend=0),colour="black")+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.225))

trsnf.lin.grwth<-ggplot(nod.tr, aes(x=(tr.Growth), y=Nfixd.kg.ha))+
  geom_point(colour="black", size=3)+
  xlab(expression("Tree Growth"~(cm^{2}~yr^{-1})))+
  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=400),colour="black")+
  geom_segment(aes(x=0,xend=.2,y=0,yend=0),colour="black")+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.225))  
  
#trsnf.comb<-trsnf.grwth+annotation_custom(ggplotGrob(trsnf.lin.grwth),xmin=.1,xmax=.21,ymin=2.1,ymax=6.2)

#pltsnf.bai<-ggplot(pltsm, aes(x=(BAI/10000), y=snf.kg.ha.gmn))+
#  #geom_smooth(method="lm",colour="black")+
#  geom_point(colour="black", size=5)+
#  geom_errorbar(aes(ymin=(snf.kg.ha.gse.lo), ymax=(snf.kg.ha.gse.hi),width=.05))+
#  #geom_errorbarh(aes(xmin=(n.tot-n.tot.se),xmax=(n.tot+n.tot.se),width=.3))+
#  ylab(expression("SNF"~(kg~N~ha^{-1}~yr^{-1})))+
#  xlab(expression("BAI"~(m^{2}~ha^{-1}~yr^{-1})))+
#  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
#        axis.text.y=element_text(size=20, colour="black"))+
#  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
#  geom_segment(aes(x=.4,xend=.4,y=0,yend=12),colour="black")+
#  geom_segment(aes(x=.4,xend=1.6,y=0,yend=0),colour="black")+
#  ggtitle("c)")+
#  theme(plot.title=element_text(hjust=-.15))


#png(filename = "SNF vs. Tree Growth and Crowding_linear axes.png", width=6, height=8, units="in", res=300)
multiplot(fixncisnf.ln, trsnf.lin.grwth, cols=1)
#dev.off()

#### Unused version of Main Figure 1 #################################
#### Bar graphs of growth, N inputs, and N availability ##############
######################################################################

tba.fig<-ggplot(t1, aes(x=as.factor(stand.age), y=BA))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(BA-BA.se), ymax=(BA+BA.se)),width=.1)+
  xlab("Forest Age")+
  ylab(expression("BA"~(m^{-2}~ha^{-1})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=40),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  annotate("text", x=c(1,2,3,4,5), y=c(10,10,10,10,10), label=c("ab","a","b","ab","b"), col=1,size=6)+
  ggtitle("a)")+
  theme(plot.title=element_text(hjust=-.1))

tbai.fig<-ggplot(t1, aes(x=as.factor(stand.age), y=BAI))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(BAI-BAI.se), ymax=(BAI+BAI.se)),width=.1)+
  xlab("Forest Age")+
  ylab(expression("BAI"~(m^{-2}~ha^{-1})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=1.5),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  annotate("text", x=c(1,2,3,4,5), y=.3, label=c("a","a","a","a","b"), col=1,size=6)+
  ggtitle("b)")+
  theme(plot.title=element_text(hjust=-.1))

anf.fig<-ggplot(asym, aes(x=as.factor(stand.age), y=Nfixd.m2))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(Nfixd.m2-Nfixd.m2.se), ymax=(Nfixd.m2+Nfixd.m2.se)),width=.1)+
  xlab("Forest Age")+
  ylab(expression("ANF"~(gN~m^{-2}~yr^{-1})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=1),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  annotate("text", x=c(1,2,3,4,5), y=.6, label=c("b","ab","a","b","ab"), col=1,size=6)+
  ggtitle("e)")+
  theme(plot.title=element_text(hjust=-.1))

nod.fig<-ggplot(nod, aes(x=as.factor(age), y=nod.bio))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(nod.bio-nod.bio.se), ymax=(nod.bio+nod.bio.se)),width=.1)+
  xlab("Forest Age")+
  ylab(expression("Nodulation"~(g~m^{-2})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=8),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  annotate("text", x=c(1,2,3,4,5), y=8, label=c("a","a","a","a","a"), col=1,size=6)+
  ggtitle("c)")+
  theme(plot.title=element_text(hjust=-.1))

snf.fig<-ggplot(sym, aes(x=as.factor(age), y=snf.kg.ha))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(snf.kg.ha-snf.kg.ha.se), ymax=(snf.kg.ha+snf.kg.ha.se)),width=.1)+
  xlab("Forest Age")+
  ylab(expression("SNF"~(kgN~ha^{-1}~yr^{-1})))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=20),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  annotate("text", x=c(1,2,3,4,5), y=18, label=c("a","a","a","a","a"), col=1,size=6)+
  ggtitle("d)")+
  theme(plot.title=element_text(hjust=-.1))

sn.fig<-ggplot(sn, aes(x=as.factor(age), y=n.tot))+
  geom_bar(stat="identity", width=.5, fill="grey", colour="black")+
  geom_errorbar(aes(ymin=(n.tot-n.tot.se), ymax=(n.tot+n.tot.se)),width=.1)+
  xlab("Forest Age")+
  ylab(expression("Soil N"~(gN~kg^{-1}~soil)))+
  scale_x_discrete(labels=c("19","19","29","37","OG"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=50),colour="black")+
  geom_segment(aes(x=0,xend=6,y=0,yend=0),colour="black")+
  annotate("text", x=c(1,2,3,4,5), y=c(10,10,10,10,10), label=c("a","d","c","c","b"), col=1,size=6)+
  ggtitle("f)")+
  theme(plot.title=element_text(hjust=-.1))

png(filename = "Plot-Level Successional N Dynamics.png", width=8, height=9, units="in", res=300)
multiplot(tba.fig,nod.fig,anf.fig,tbai.fig,snf.fig,sn.fig, cols=2)
dev.off()