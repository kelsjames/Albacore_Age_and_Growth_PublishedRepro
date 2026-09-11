#Compare male vs female for best model (prob VBGF)
rm(list = ls())
#Growth modeling
library(openxlsx)
library(ggplot2)
library(FSA)
library(FSAdata)
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
library(lmtest)

alb <- read.xlsx("Data/alb.xlsx")
alb$AgeSep <- alb$Final.ID1+(123/365)
alb$Sex <- as.factor(alb$Sex)
datafm <- subset(alb, Sex == "Female" | Sex == "Male")

#name data
data1=datafm
sl=data1$FL
age=data1$AgeSep

#name functions from FSA
vb1 <- makeGrowthFun(type = "von Bertalanffy", pname="Typical")

#starting values
svvb <- list(Linf=120,K=0.15,t0=-2)

####von bertalanffy 3-parameter sexes combined - unknowns excluded####
fitvb <- nlsLM(sl~vb1(age,Linf,K,t0),data=data1,start=svvb,trace=TRUE)
vbvals=coef(fitvb)
summary(fitvb)
bootvb <- nlsBoot(fitvb)
cbind(coef(fitvb),confint(bootvb))

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
preds1vb$Model <- "Sexes Combined"


####von bertalanffy 3-parameter sexes separate####
svsex=lapply(svvb,rep,2) #replication of starting values for other models becuase model by Sex needs 6 starting parameters, not just 3

fitfvm = nls(sl~vb1(age,Linf[Sex],K[Sex],t0[Sex]),data=data1,start=svsex,trace=TRUE)
vbfvmvals=coef(fitfvm)
summary(fitfvm)
bootfvm <- nlsBoot(fitfvm)
cbind(coef(fitfvm),confint(bootfvm))

# predicting for plotting, have to predict M & F separately
#Females
fem <- subset(alb, Sex == "Female")
data3=fem
sl3=data3$FL
age3=data3$AgeSep

fitf = nlsLM(sl3~vb1(age3,Linf,K,t0),data=data3,start=svvb,trace=TRUE)
summary(fitf)
bootf <- nlsBoot(fitf)
cbind(coef(fitf),confint(bootf))
predict2f <- function(x) predict(x,data.frame(age3=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals for females
ages <- seq(1.3,13.3,by=0.2)
f.boot2f <- Boot(fitf,f=predict2f)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1f <- data.frame(ages,
                       predict(fitf,data.frame(age3=ages)),
                       confint(f.boot2f))
names(preds1f) <- c("age","fit","LCI","UCI")
headtail(preds1f)
preds1f$Model <- "Females"

#Males
mal <- subset(alb, Sex == "Male")
data4=mal
sl4=data4$FL
age4=data4$AgeSep

fitm = nlsLM(sl4~vb1(age4,Linf,K,t0),data=data4,start=svvb,trace=TRUE)
summary(fitm)
bootm <- nlsBoot(fitm)
cbind(coef(fitm),confint(bootm))
predict2m <- function(x) predict(x,data.frame(age4=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals for females
ages <- seq(1.3,15.3,by=0.2)
f.boot2m <- Boot(fitm,f=predict2m)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1m <- data.frame(ages,
                        predict(fitm,data.frame(age4=ages)),
                        confint(f.boot2m))
names(preds1m) <- c("age","fit","LCI","UCI")
headtail(preds1m)
preds1m$Model <- "Males"

#Get Reg VB on all data (including unknowns) for supplemental figure
data5=alb
sl5=data5$FL
age5=data5$AgeSep

fitvbr = nlsLM(sl5~vb1(age5,Linf,K,t0),data=data5,start=svvb,trace=TRUE)
summary(fitvbr)
bootvbr <- nlsBoot(fitvbr)
cbind(coef(fitvbr),confint(bootvbr))
predict2vbr <- function(x) predict(x,data.frame(age5=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals for females
ages <- seq(1.3,15.3,by=0.2)
f.boot2vbr <- Boot(fitvbr,f=predict2vbr)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vbr <- data.frame(ages,
                      predict(fitvbr,data.frame(age5=ages)),
                      confint(f.boot2vbr))
names(preds1vbr) <- c("age","fit","LCI","UCI")
headtail(preds1vbr)
preds1vbr$Model <- "Sexes Combined with Unknowns"

####Compare####
#create AICc table (to get AIC add secord.ord = FALSE, default is TRUE, which is AICc)
aictab(list(fitvb,fitfvm),c("Sexes Combined","Sexes Separate"))

#obtain BIC values
BIC(fitvb,fitfvm)

#Likelihood ratio test
lrtest(fitvb,fitfvm) #The p-value is not significant, so we should use the nested model (fitvb)
#rather than the complex model (fitfvm). The complex model does not out-perform the nested model.
#In other words the LogLik of fvm is higher than the LogLik of vb (-949 vs -951) so the fvm fits better,
#but not better enough at the expense of extra degrees of freedom (7 vs 4). Hence the p > 0.05
#lrt(fitvb,fitfvm) - since comparing only one model (fvm) to the nested model (vb), lrtest above
#is fine to use rather than the lrt() from FSA. lrt() needs an input of a common model that I don't have here
#Extra sum-of-squares test - assumes the underlying distribution is normal
#extraSS(fitvb,fitfvm)

#view parameter ests
vbvals
vbfvmvals

####Plot####
# plotting with ggplot for Figure 11
preds1.all <- rbind(preds1vb,preds1f,preds1m)
preds1.all$Model <- factor(preds1.all$Model,levels=c("Sexes Combined","Females","Males"))

col.pal <- c("red","blue","gray30","red","blue")
col.pal1 <- c("gray30", "red", "blue")
datafm <- datafm %>%
  mutate(Model = case_when(Sex == "Female" ~ "Females",
                           Sex == "Male" ~ "Males"))

ggplot() + 
  geom_point(data=datafm,aes(y=FL,x=AgeSep, color = Sex),size=2,alpha=0.1, show.legend = F)+
  geom_line(data=preds1.all,aes(y=fit,x=age, color=Model,linetype=Model,group=Model),linewidth=1, show.legend = F)+ 
  geom_ribbon(data=preds1.all,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.2)+
  scale_color_manual(values=col.pal)+
  scale_fill_manual(values=col.pal1)+
  scale_linetype_manual(values=c("dashed","solid","solid"))+
  guides(fill=guide_legend(override.aes=list(alpha=1)))+
  labs(x="Age (years)", y = "Fork Length (cm)")+
  scale_x_continuous(limits=c(0.5,15.5),breaks=seq(0,15,1))+
  scale_y_continuous(limits=c(45,125),breaks=seq(20,140,10))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.3),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

ggsave("Output/SexSpecificVB.jpeg", height=6, width=8)


# plot sexes combined with all fish (including unknowns) and just M&F (excluding unknowns) with ggplot for supplemental figure
preds5.all <- rbind(preds1vb,preds1vbr)
preds5.all$Model <- factor(preds5.all$Model,levels=c("Sexes Combined","Sexes Combined with Unknowns"))

ggplot() + 
  geom_line(data=preds5.all,aes(y=fit,x=age, color=Model, linetype=Model, group=Model),linewidth=1, show.legend = F)+ 
  geom_ribbon(data=preds5.all,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.2)+
  scale_color_manual(values=c("#E880DA","black"))+ 
  scale_fill_manual(values=c("#E880DA","black"))+
  scale_linetype_manual(values=c("solid","dashed"))+
  guides(fill=guide_legend(override.aes=list(alpha=1)))+
  labs(x="Age (years)", y = "Fork Length (cm)")+
  scale_x_continuous(limits=c(0.5,15.5),breaks=seq(0,15,1))+
  scale_y_continuous(limits=c(45,125),breaks=seq(20,140,10))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.3),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

ggsave("Output/SexesCombinedVBUnk.jpeg", height=6, width=8)

