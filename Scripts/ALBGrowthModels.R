rm(list = ls())
#Growth modeling
library(openxlsx)
library(ggplot2)
library(FSA)
library(nlme)
library(mgcv)
library(car)
library(AICcmodavg)
library(lubridate)
library(minpack.lm)
library(nlstools)
library(patchwork)
library(tidyverse)
library(viridis)

#load in dataset with both ages and lengths
alb <- read.xlsx("Data/alb.xlsx")

#name functions from FSA
l1 <- makeGrowthFun(type="logistic")
g1 <- makeGrowthFun(type="Gompertz",pname="Ricker1")
vb1 <- makeGrowthFun(type = "von Bertalanffy", pname="Typical")

#To generate columns for decimal age and Sept age
alb <- alb %>%
  mutate(
    yday=yday(make_date(year=Collection.Year,month=Collection.Month,day=Collection.Day)),
    yday_may1 = as.numeric((yday-121) %% 365) + 1,
    Dec = (yday_may1-1)/365,
    AgeSep = Final.ID1+(123/365),
    AgeAdj = ifelse(is.na(Dec),AgeSep,Dec + Final.ID1),
  )

#####Compare Whole age vs Adj Age vs Sept Age for VBGF####
#VBGF
#For whole ages
data1=alb
sl=data1$FL
age=data1$Final.ID1
svvb <- list(Linf=130,K=0.3,t0=-2)
fitvbwh <- nlsLM(sl~vb1(age,Linf,K,t0),data=data1,start=svvb,trace=TRUE)
summary(fitvbwh)
bootvbwh <- nlsBoot(fitvbwh)
cbind(coef(fitvbwh),confint(bootvbwh))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbwh <- Boot(fitvbwh,f=predict2vb)
preds1vbwh <- data.frame(ages,
                         predict(fitvbwh,data.frame(age=ages)),
                         confint(f.boot2vbwh))
names(preds1vbwh) <- c("age","fit","LCI","UCI")
headtail(preds1vbwh)
preds1vbwh$Model <- "Whole Age"

#For age adjusted from May 1
age=data1$AgeAdj
fitvbdec <- nlsLM(sl~vb1(age,Linf,K,t0),data=data1,start=svvb,trace=TRUE)
summary(fitvbdec)
bootvbdec <- nlsBoot(fitvbdec)
cbind(coef(fitvbdec),confint(bootvbdec))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbdec <- Boot(fitvbdec,f=predict2vb)
preds1vbdec <- data.frame(ages,
                          predict(fitvbdec,data.frame(age=ages)),
                          confint(f.boot2vbdec))
names(preds1vbdec) <- c("age","fit","LCI","UCI")
headtail(preds1vbdec)
preds1vbdec$Model <- "Decimal Age"


#For AgeSep (all fish on Sep 1st)
age=data1$AgeSep
fitvbsep <- nlsLM(sl~vb1(age,Linf,K,t0),data=data1,start=svvb,trace=TRUE)
summary(fitvbsep)
bootvbsep <- nlsBoot(fitvbsep)
cbind(coef(fitvbsep),confint(bootvbsep))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbsep <- Boot(fitvbsep,f=predict2vb)
preds1vbsep <- data.frame(ages,
                          predict(fitvbsep,data.frame(age=ages)),
                          confint(f.boot2vbsep))
names(preds1vbsep) <- c("age","fit","LCI","UCI")
headtail(preds1vbsep)
preds1vbsep$Model <- "Sept 1st"

#run sexes combined model above to get fitvb
aictab(list(fitvbwh,fitvbdec,fitvbsep),c("Whole Age","Decimal Age", "Sept 1st"))
BIC(fitvbwh,fitvbdec,fitvbsep)

predsdec <- rbind(preds1vbwh,preds1vbdec,preds1vbsep)
predsdec$Model <- factor(predsdec$Model, levels=c("Whole Age","Decimal Age","Sept 1st"))
col.pal3 <- c("#331D9E","#719AC7","#E67C82")
#plot of Whole Age, Adj age, and Sept age, not in publication
ggplot() + 
  geom_point(data=alb,aes(y=FL,x=AgeAdj),size=2,alpha=0.1)+
  geom_line(data=predsdec,aes(y=fit,x=age,color=Model),linewidth=1)+ #line for aged fish only
  #geom_ribbon(data=predsdec,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.5)+
  scale_color_manual(values=col.pal3)+
  scale_fill_manual(values=col.pal3)+
  labs(x="Age (years)", y = "Fork Length (cm)")+
  scale_x_continuous(limits=c(0,16),breaks=seq(0,15,1))+
  scale_y_continuous(limits=c(20,140),breaks=seq(20,140,10))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.2),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

####Compare Whole age vs Adj Age vs Sept Age for Gompertz####
#Not changing the names of the fit models due to effort
#GOMPERTZ
#Whole ages
data1=alb
sl=data1$FL
age=data1$Final.ID1
svvb <- list(Linf=130,gi=0.5,ti=0.2)
fitvbwh <- nlsLM(sl~g1(age,Linf,gi,ti),data=data1,start=svvb,trace=TRUE)
summary(fitvbwh)
bootvbwh <- nlsBoot(fitvbwh)
cbind(coef(fitvbwh),confint(bootvbwh))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbwh <- Boot(fitvbwh,f=predict2vb)
preds1vbwh <- data.frame(ages,
                         predict(fitvbwh,data.frame(age=ages)),
                         confint(f.boot2vbwh))
names(preds1vbwh) <- c("age","fit","LCI","UCI")
headtail(preds1vbwh)
preds1vbwh$Model <- "Whole Age"

#Age adjusted at May 1 (decimal)
age=data1$AgeAdj
fitvbdec <- nlsLM(sl~g1(age,Linf,gi,ti),data=data1,start=svvb,trace=TRUE)
summary(fitvbdec)
bootvbdec <- nlsBoot(fitvbdec)
cbind(coef(fitvbdec),confint(bootvbdec))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbdec <- Boot(fitvbdec,f=predict2vb)
preds1vbdec <- data.frame(ages,
                          predict(fitvbdec,data.frame(age=ages)),
                          confint(f.boot2vbdec))
names(preds1vbdec) <- c("age","fit","LCI","UCI")
headtail(preds1vbdec)
preds1vbdec$Model <- "Decimal Age"

#For AgeSep (all fish on Sep 1st)
age=data1$AgeSep
fitvbsep <- nlsLM(sl~g1(age,Linf,gi,ti),data=data1,start=svvb,trace=TRUE)
summary(fitvbsep)
bootvbsep <- nlsBoot(fitvbsep)
cbind(coef(fitvbsep),confint(bootvbsep))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbsep <- Boot(fitvbsep,f=predict2vb)
preds1vbsep <- data.frame(ages,
                          predict(fitvbsep,data.frame(age=ages)),
                          confint(f.boot2vbsep))
names(preds1vbsep) <- c("age","fit","LCI","UCI")
headtail(preds1vbsep)
preds1vbsep$Model <- "Sept 1st"

#run sexes combined model above to get fitvb
aictab(list(fitvbwh,fitvbdec,fitvbsep),c("Whole Age","Decimal Age", "Sept 1st"))
BIC(fitvbwh,fitvbdec,fitvbsep)


####Compare Whole age vs Adj Age vs Sept Age for Logistic####
#Not changing the names of the fit models due to effort
#LOGISTIC
#Whole ages
data1=alb
sl=data1$FL
age=data1$Final.ID1
svvb <- list(Linf=130,gninf=0.5,ti=2)
fitvbwh <- nlsLM(sl~l1(age,Linf,gninf,ti),data=data1,start=svvb,trace=TRUE)
summary(fitvbwh)
bootvbwh <- nlsBoot(fitvbwh)
cbind(coef(fitvbwh),confint(bootvbwh))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbwh <- Boot(fitvbwh,f=predict2vb)
preds1vbwh <- data.frame(ages,
                         predict(fitvbwh,data.frame(age=ages)),
                         confint(f.boot2vbwh))
names(preds1vbwh) <- c("age","fit","LCI","UCI")
headtail(preds1vbwh)
preds1vbwh$Model <- "Whole Age"

#Decimal ages adjusted to May 1
age=data1$AgeAdj
fitvbdec <- nlsLM(sl~l1(age,Linf,gninf,ti),data=data1,start=svvb,trace=TRUE)
summary(fitvbdec)
bootvbdec <- nlsBoot(fitvbdec)
cbind(coef(fitvbdec),confint(bootvbdec))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbdec <- Boot(fitvbdec,f=predict2vb)
preds1vbdec <- data.frame(ages,
                          predict(fitvbdec,data.frame(age=ages)),
                          confint(f.boot2vbdec))
names(preds1vbdec) <- c("age","fit","LCI","UCI")
headtail(preds1vbdec)
preds1vbdec$Model <- "Decimal Age"

#For AgeSep (all fish on Sep 1st)
age=data1$AgeSep
fitvbsep <- nlsLM(sl~l1(age,Linf,gninf,ti),data=data1,start=svvb,trace=TRUE)
summary(fitvbsep)
bootvbsep <- nlsBoot(fitvbsep)
cbind(coef(fitvbsep),confint(bootvbsep))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbsep <- Boot(fitvbsep,f=predict2vb)
preds1vbsep <- data.frame(ages,
                          predict(fitvbsep,data.frame(age=ages)),
                          confint(f.boot2vbsep))
names(preds1vbsep) <- c("age","fit","LCI","UCI")
headtail(preds1vbsep)
preds1vbsep$Model <- "Sept 1st"

#run sexes combined model above to get fitvb
aictab(list(fitvbwh,fitvbdec,fitvbsep),c("Whole Age","Decimal Age", "Sept 1st"))
BIC(fitvbwh,fitvbdec,fitvbsep)


####Compare Whole age vs Adj Age vs Sept Age for Richards####
#Not changing the names of the fit models due to effort
#RICHARDS
#Whole Ages
data1=alb
sl=data1$FL
age=data1$Final.ID1
svvb <- list(Linf=130, k=0.5, a=1, b=15)
fitvbwh <- nlsLM(sl~Linf*(1-a*exp(-k*age))^b,data=data1,start=svvb,trace=TRUE)
summary(fitvbwh)
bootvbwh <- nlsBoot(fitvbwh)
cbind(coef(fitvbwh),confint(bootvbwh))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbwh <- Boot(fitvbwh,f=predict2vb)
preds1vbwh <- data.frame(ages,
                         predict(fitvbwh,data.frame(age=ages)),
                         confint(f.boot2vbwh))
names(preds1vbwh) <- c("age","fit","LCI","UCI")
headtail(preds1vbwh)
preds1vbwh$Model <- "Whole Age"

#Decimal ages from May 1
age=data1$AgeAdj
fitvbdec <- nlsLM(sl~Linf*(1-a*exp(-k*age))^b,data=data1,start=svvb,trace=TRUE)
summary(fitvbdec)
bootvbdec <- nlsBoot(fitvbdec)
cbind(coef(fitvbdec),confint(bootvbdec))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbdec <- Boot(fitvbdec,f=predict2vb)
preds1vbdec <- data.frame(ages,
                          predict(fitvbdec,data.frame(age=ages)),
                          confint(f.boot2vbdec))
names(preds1vbdec) <- c("age","fit","LCI","UCI")
headtail(preds1vbdec)
preds1vbdec$Model <- "Decimal Age"

#For AgeSep (all fish on Sep 1st)
age=data1$AgeSep
fitvbsep <- nlsLM(sl~Linf*(1-a*exp(-k*age))^b,data=data1,start=svvb,trace=TRUE)
summary(fitvbsep)
bootvbsep <- nlsBoot(fitvbsep)
cbind(coef(fitvbsep),confint(bootvbsep))
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- seq(1,15,by=0.2)
f.boot2vbsep <- Boot(fitvbsep,f=predict2vb)
preds1vbsep <- data.frame(ages,
                          predict(fitvbsep,data.frame(age=ages)),
                          confint(f.boot2vbsep))
names(preds1vbsep) <- c("age","fit","LCI","UCI")
headtail(preds1vbsep)
preds1vbsep$Model <- "Sept 1st"


#run sexes combined model above to get fitvb
aictab(list(fitvbwh,fitvbdec,fitvbsep),c("Whole Age","Decimal Age", "Sept 1st"))
BIC(fitvbwh,fitvbdec,fitvbsep)

####Running Growth models on Sept 1st, a little redundant, but already set up for graphing and residuals####
#Fitting growth models
#name data
data1=alb
sl=data1$FL
#age=data1$Final.ID1
age=data1$AgeSep

#gompertz
svG1 <- list(Linf=130,gi=0.5,ti=0.2)
fitG1 <- nlsLM(sl~g1(age,Linf,gi,ti),data=data1,start=svG1,trace=TRUE)
gvals=coef(fitG1)
summary(fitG1)
bootG1 <- nlsBoot(fitG1)
cbind(coef(fitG1),confint(bootG1))
#Get residuals for plotting, plot at bottom
resG1 <- data.frame(
  fitted=fitted(fitG1),
  residuals=resid(fitG1)/summary(fitG1)$sigma
)

# predicting for plotting
predict2G1 <- function(x) predict(x,data.frame(age=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(1.3,15.3,by=0.2)
f.boot2G1 <- Boot(fitG1,f=predict2G1)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1G1 <- data.frame(ages,
                       predict(fitG1,data.frame(age=ages)),
                       confint(f.boot2G1))
names(preds1G1) <- c("age","fit","LCI","UCI")
headtail(preds1G1)
preds1G1$Model <- "Gompertz"


#logistic
svL1 <- list(Linf=130,gninf=0.5,ti=2)
fitL1 <- nlsLM(sl~l1(age,Linf,gninf,ti),data=data1,start=svL1,trace=TRUE)
lvals=coef(fitL1)
summary(fitL1)
bootL1 <- nlsBoot(fitL1)
cbind(coef(fitL1),confint(bootL1))
#Get residuals for plotting, plot at bottom
resL1 <- data.frame(
  fitted=fitted(fitL1),
  residuals=resid(fitL1)/summary(fitL1)$sigma
)

# predicting another way for plotting
predict2L1 <- function(x) predict(x,data.frame(age=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(1.3,15.3,by=0.2)
f.boot2L1 <- Boot(fitL1,f=predict2L1)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1L1 <- data.frame(ages,
                       predict(fitL1,data.frame(age=ages)),
                       confint(f.boot2L1))
names(preds1L1) <- c("age","fit","LCI","UCI")
headtail(preds1L1)
preds1L1$Model <- "Logistic"


#von bertalanffy 3-parameter
svvb <- list(Linf=130,K=0.3,t0=-2)
fitvb <- nlsLM(sl~vb1(age,Linf,K,t0),data=data1,start=svvb,trace=TRUE)
vbvals=coef(fitvb)
summary(fitvb)
bootvb <- nlsBoot(fitvb)
cbind(coef(fitvb),confint(bootvb))
#Get residuals for plotting, plot at bottom
resvb1 <- data.frame(
  fitted=fitted(fitvb),
  residuals=resid(fitvb)/summary(fitvb)$sigma
)

# predicting another way for plotting
predict2vb <- function(x) predict(x,data.frame(age=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(1.3,15.3,by=0.2)
f.boot2vb <- Boot(fitvb,f=predict2vb)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb <- data.frame(ages,
                       predict(fitvb,data.frame(age=ages)),
                       confint(f.boot2vb))
names(preds1vb) <- c("age","fit","LCI","UCI")
headtail(preds1vb)
preds1vb$Model <- "von Bertalanffy"


#Richards
svr1 <- list(Linf=130, k=0.5, a=1, b=15)
fitR1 <- nlsLM(sl~Linf*(1-a*exp(-k*age))^b,data=data1,start=svr1,trace=TRUE)
summary(fitR1)
bootR1 <- nlsBoot(fitR1)
cbind(coef(fitR1),confint(bootR1))
rvals=coef(fitR1)
#Get residuals for plotting, plot at bottom
resR1 <- data.frame(
  fitted=fitted(fitR1),
  residuals=resid(fitR1)/summary(fitR1)$sigma
)


# predicting another way for plotting
predict2R1 <- function(x) predict(x,data.frame(age=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(1.3,15.3,by=0.2)
f.boot2R1 <- Boot(fitR1,f=predict2R1)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1R1 <- data.frame(ages,
                       predict(fitR1,data.frame(age=ages)),
                       confint(f.boot2R1))
names(preds1R1) <- c("age","fit","LCI","UCI")
headtail(preds1R1)
preds1R1$Model <- "Richards"


#create AIC table
aictab(list(fitvb,fitL1,fitG1,fitR1),c("von Bertalanffy","logistic","Gompertz","Richards"))

#obtain BIC values
BIC(fitL1,fitG1,fitvb,fitR1)

#view parameter ests
vbvals
gvals
lvals
rvals

# plotting with ggplot
preds3.all <- rbind(preds1vb,preds1G1,preds1L1,preds1R1)
preds3.all$Model <- factor(preds3.all$Model, levels=c("von Bertalanffy","Gompertz","Logistic","Richards"))

col.pal2 <- c("#332288","#117733","#DDCC77","#661100")

#plot of VB, Gompertz, Logistic, AND Richards, not a figure in the paper
ggplot() + 
  geom_point(data=alb,aes(y=FL,x=AgeSep),size=2,alpha=0.1)+
  geom_line(data=preds3.all,aes(y=fit,x=age,color=Model),linewidth=1,linetype=2)+ #line for aged fish only
  geom_ribbon(data=preds3.all,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.5)+
  scale_color_manual(values=col.pal2)+
  scale_fill_manual(values=col.pal2)+
  labs(x="Age (years)", y = "Fork Length (cm)")+
  scale_x_continuous(limits=c(0,16),breaks=seq(0,15,1))+
  scale_y_continuous(limits=c(20,140),breaks=seq(20,140,10))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.2),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

ggsave("Output/growthmodelsandrichards.jpeg", height=5, width=6)

#Plotting residuals####
sG1 <- ggplot(resG1, aes(fitted,residuals))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dashed")+
  labs(x="Fitted", y="Standardized Residuals")+
  theme_classic(base_size=14)
hG1 <- ggplot(resG1,aes(residuals))+
  geom_histogram(fill="grey",color="black")+
  labs(x="Standardized Residuals",y="Count")+
  theme_classic(base_size=14)

svb1 <- ggplot(resvb1, aes(fitted,residuals))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dashed")+
  labs(x="Fitted", y="Standardized Residuals")+
  theme_classic(base_size=14)
hvb1 <- ggplot(resvb1,aes(residuals))+
  geom_histogram(fill="grey",color="black")+
  labs(x="Standardized Residuals",y="Count")+
  theme_classic(base_size=14)

sL1 <- ggplot(resL1, aes(fitted,residuals))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dashed")+
  labs(x="Fitted", y="Standardized Residuals")+
  theme_classic(base_size=14)
hL1 <- ggplot(resL1,aes(residuals))+
  geom_histogram(fill="grey",color="black")+
  labs(x="Standardized Residuals",y="Count")+
  theme_classic(base_size=14)

sR1 <- ggplot(resR1, aes(fitted,residuals))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dashed")+
  labs(x="Fitted", y="Standardized Residuals")+
  theme_classic(base_size=14)
hR1 <- ggplot(resR1,aes(residuals))+
  geom_histogram(fill="grey",color="black")+
  labs(x="Standardized Residuals",y="Count")+
  theme_classic(base_size=14)

design <-'
AB
CD
EF
GH'
wrap_plots(A=svb1,B=hvb1,C=sR1,D=hR1,E=sG1,F=hG1,G=sL1,H=hL1,design=design) + plot_annotation(tag_levels=list(c("A","","B","","C","","D","")))
ggsave("Output/growthmodelresidualerrors.jpeg", height=20, width=14)

