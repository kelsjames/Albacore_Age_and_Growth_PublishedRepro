rm(list = ls())
#Compare to previous studies, involves comparing models
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

#Bring in age and length data
alb <- read.xlsx("Data/alb.xlsx")
alb$AgeSep <- alb$Final.ID1+(123/365)

#bring in other age data
ages.other <- read.csv("Data/AlbacoreRpub.csv")

#combine all data, add AgeSep for Snodgrass
alb$Source <- "Snodgrass"
alb$Region <- ifelse(alb$Collection.Region == "CPO","Central","Eastern")
alb.trans <- alb %>% select(Source, ALB_ID, AgeSep, Final.ID1, FL, Collection.Year, 
                            Region, Sex, Collection.Month, Collection.Day)
colnames(alb.trans) <- c("Source", "ID", "Ageadj", "Age", "FL", "Year", "Region", "Sex", "Month", "Day")
ages.other.trans <- ages.other %>% select(Source, ID, Ageadj, Age, FL, Year, Region, Sex, Month, Day)
alb.all <- rbind(ages.other.trans, alb.trans)


##Set up growth models for previous studies
# setting up von Bertalanffy parameters for all models
vb1 <- vbFuns()
svvb <- list(Linf=250,K=0.3,t0=-2)

#Use below as an example and modify for comparison
#Can't do AIC or BIC since different datasets, but overlap of CI are a good indicator of 'different'
# Snodgrass data set (this study) ###### 
data1=alb
sl=data1$FL
age=data1$AgeSep

fitvb.org <- nlsLM(sl~vb1(age,Linf,K,t0),data=data1,start=svvb,trace=TRUE)
vbvals.org=coef(fitvb.org)
summary(fitvb.org)
bootvb.org <- nlsBoot(fitvb.org)
vbvals.org <- as.data.frame(t(cbind(coef(fitvb.org),confint(bootvb.org))))
vbvals.org$Dataset <- "Snodgrass"
rownames(vbvals.org)<-c("coeff","LCI","UCI")
vbvals.org$Parameter <- rownames(vbvals.org)

# predicting another way for plotting
predict2vb.org <- function(x) predict(x,data.frame(age=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(1.3,15.3,by=0.2)
f.boot2vb.org <- Boot(fitvb.org,f=predict2vb.org)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb.org <- data.frame(ages,
                           predict(fitvb.org,data.frame(age=ages)),
                           confint(f.boot2vb.org))
names(preds1vb.org) <- c("age","fit","LCI","UCI")
headtail(preds1vb.org)
preds1vb.org$Model <- "This Study"

# Wells data set ###### 
age.wells <- alb.all %>% filter(Source=="Wells")
data2=age.wells
sl2=data2$FL
age2=data2$Ageadj

fitvb.wells <- nlsLM(sl2~vb1(age2,Linf,K,t0),data=data2,start=svvb,trace=TRUE)
vbvals.wells=coef(fitvb.wells)
summary(fitvb.wells)
bootvb.wells <- nlsBoot(fitvb.wells)
vbvals.wells <- as.data.frame(t(cbind(coef(fitvb.wells),confint(bootvb.wells))))
vbvals.wells$Dataset <- "Wells"
rownames(vbvals.wells)<-c("coeff","LCI","UCI")
vbvals.wells$Parameter <- rownames(vbvals.wells)

# predicting another way for plotting
predict2vb.wells <- function(x) predict(x,data.frame(age2=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(1.3,15.3,by=0.2)
f.boot2vb.wells <- Boot(fitvb.wells,f=predict2vb.wells)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb.wells <- data.frame(ages,
                           predict(fitvb.wells,data.frame(age2=ages)),
                           confint(f.boot2vb.wells))
names(preds1vb.wells) <- c("age","fit","LCI","UCI")
headtail(preds1vb.wells)
preds1vb.wells$Model <- "Wells et al. 2013"

# Comparing output of this study to Wells for Figure 10######
preds1.snod.wells <- rbind(preds1vb.org,preds1vb.wells)
preds1.snod.wells$Model <- factor(preds1.snod.wells$Model,levels=c("This Study","Wells et al. 2013"))
col.pal7 <- c("black","#249846","black","#249846")


alb.snodwells <- subset(alb.all, Source != "Chen")
alb.snodwells$S1 <- ifelse(alb.snodwells$Source == "Wells", "Wells et al. 2013", "This Study")
alb.snodwells$S1 <- as.factor(alb.snodwells$S1)
ggplot() + 
  geom_point(data=alb.snodwells,aes(y=FL,x=Ageadj, color=S1), size=2,alpha=0.3)+
  geom_line(data=preds1.snod.wells,aes(y=fit,x=age, color=Model,linetype=Model,group=Model),linewidth=1,show.legend = FALSE)+ #line only for observed ages
  geom_ribbon(data=preds1.snod.wells,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.2)+
  scale_color_manual(values=col.pal7)+
  scale_fill_manual(values=col.pal7)+
  scale_linetype_manual(values=c("dashed","solid"))+
  guides(fill=guide_legend(override.aes=list(alpha=1)))+
  labs(x="Age (years)", y = "Fork length (cm)",color="Dataset",fill="Dataset")+
  scale_x_continuous(limits=c(0.5,15.5),breaks=seq(0,15,1))+
  scale_y_continuous(limits=c(45,125),breaks=seq(20,140,10))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.2),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

ggsave("Output/GrowthModelOwynWells.jpeg", height=6, width=8)


# Snodgrass F v M (this study)####
# predicting for plotting, have to predict M & F separately
#Females
fem <- subset(alb, Sex == "Female")
data4=fem
sl4=data4$FL
age4=data4$AgeSep

fitf = nlsLM(sl4~vb1(age4,Linf,K,t0),data=data4,start=svvb,trace=TRUE)
vbvals.f=coef(fitf)
summary(fitf)
bootf <- nlsBoot(fitf)
cbind(coef(fitf),confint(bootf))
predict2f <- function(x) predict(x,data.frame(age4=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals for females
ages <- seq(1.3,13.3,by=0.2)
f.boot2f <- Boot(fitf,f=predict2f)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1f <- data.frame(ages,
                      predict(fitf,data.frame(age4=ages)),
                      confint(f.boot2f))
names(preds1f) <- c("age","fit","LCI","UCI")
headtail(preds1f)
preds1f$Model <- "Females This Study"

#Males
mal <- subset(alb, Sex == "Male")
data5=mal
sl5=data5$FL
age5=data5$AgeSep

fitm = nlsLM(sl5~vb1(age5,Linf,K,t0),data=data5,start=svvb,trace=TRUE)
vbvals.m=coef(fitm)
summary(fitm)
bootm <- nlsBoot(fitm)
cbind(coef(fitm),confint(bootm))
predict2m <- function(x) predict(x,data.frame(age5=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals for females
ages <- seq(1.3,15.3,by=0.2)
f.boot2m <- Boot(fitm,f=predict2m)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1m <- data.frame(ages,
                      predict(fitm,data.frame(age5=ages)),
                      confint(f.boot2m))
names(preds1m) <- c("age","fit","LCI","UCI")
headtail(preds1m)
preds1m$Model <- "Males This Study"

# Chen and Wells F v M####
#predicting for plotting, have to predict M & F separately
#Females
femc <- subset(ages.other, Sex == "F")
data6=femc
sl6=data6$FL
age6=data6$Ageadj

fitfc = nlsLM(sl6~vb1(age6,Linf,K,t0),data=data6,start=svvb,trace=TRUE)
vbvals.fc=coef(fitfc)
summary(fitfc)
bootfc <- nlsBoot(fitfc)
cbind(coef(fitfc),confint(bootfc))
predict2fc <- function(x) predict(x,data.frame(age6=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals for females
ages <- seq(1.3,11.3,by=0.2)
f.boot2fc <- Boot(fitfc,f=predict2fc)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1fc <- data.frame(ages,
                      predict(fitfc,data.frame(age6=ages)),
                      confint(f.boot2fc))
names(preds1fc) <- c("age","fit","LCI","UCI")
headtail(preds1fc)
preds1fc$Model <- "Females James et al. 2020a"

#Males
malc <- subset(ages.other, Sex == "M")
data7=malc
sl7=data7$FL
age7=data7$Ageadj

fitmc = nlsLM(sl7~vb1(age7,Linf,K,t0),data=data7,start=svvb,trace=TRUE)
vbvals.mc=coef(fitmc)
summary(fitmc)
bootmc <- nlsBoot(fitmc)
cbind(coef(fitmc),confint(bootmc))
predict2mc <- function(x) predict(x,data.frame(age7=ages))
# constructing mean lengths at ages with bootstrapped confidence intervals for females
ages <- seq(1.3,14.3,by=0.2)
f.boot2mc <- Boot(fitmc,f=predict2mc)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1mc <- data.frame(ages,
                      predict(fitmc,data.frame(age7=ages)),
                      confint(f.boot2mc))
names(preds1mc) <- c("age","fit","LCI","UCI")
headtail(preds1mc)
preds1mc$Model <- "Males James et al. 2020a"

# Comparing output sexes separate For Figure 12 ######
# getting coefficients for each model
vbvals.all.fm <- rbind(vbvals.f,vbvals.m,vbvals.fc,vbvals.mc)
vbvals.all.fm <- vbvals.all.fm  %>% mutate_if(is.numeric, ~round(., 2))
#write.xlsx(vbvals.all.fm, "Output/GrowthModelCompAcrossStudies_SexSeparate.xlsx")

# plotting with ggplot
preds1.all <- rbind(preds1f,preds1m,preds1fc,preds1mc)
preds1.all$Model <- factor(preds1.all$Model,levels=c("Females This Study","Males This Study","Females James et al. 2020a","Males James et al. 2020a"))

#col.pal <- c("#CC6677","#117733","#332288","#DDCC77","#661100")
col.pal2 <- c("red","blue","#FFD700","#2FDCC0")
alb.fm <- subset(alb.all, Sex == "F" | Sex == "M" | Sex == "Female" | Sex == "Male")
alb.fm <- alb.fm %>%
  mutate(Model = case_when(Source %in% c("Chen","Wells") & Sex == "F" ~ "Females Chen & Wells",
                           Source %in% c("Chen","Wells") & Sex == "M" ~ "Males Chen & Wells",
                           Source == "Snodgrass" & Sex == "Female" ~ "Females Snodgrass",
                           Source == "Snodgrass" & Sex == "Male" ~ "Males Snodgrass"))
alb.fm$line <- ifelse(alb.fm$Source == "Snodgrass","A","B")

ggplot() + 
  #geom_point(data=alb.fm,aes(y=FL,x=Ageadj,color=Model), size=2,alpha=0.3)+
  geom_line(data=preds1.all,aes(y=fit,x=age, linetype=Model),linewidth=1)+ #line only for observed ages
  geom_line(data=preds1.all,aes(y=fit,x=age, color=Model,linetype=Model),linewidth=1)+ #line only for observed ages
  geom_ribbon(data=preds1.all,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.2,show.legend = FALSE)+
  scale_color_manual(values=col.pal2)+
  scale_fill_manual(values=col.pal2)+
  #scale_fill_manual(values=c("black","black","black","black"))+
  scale_linetype_manual(values=c("solid","dotdash","dashed","dotted"))+
  guides(fill=guide_legend(override.aes=list(alpha=1)))+
  labs(x="Age (years)", y = "Fork length (cm)",color="Dataset",fill="Dataset")+
  scale_x_continuous(limits=c(0.5,15),breaks=seq(0,15,1))+
  scale_y_continuous(limits=c(45,125),breaks=seq(20,140,10))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.3),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

ggsave("Output/GrowthModelOwynChen_SexesSeparate_wpoints.jpeg", height=6, width=8)
ggsave("Output/GrowthModelOwynChen_SexesSeparate_wopoints.jpeg", height=6, width=8)
ggsave("Output/GrowthModelOwynChen_SexesSeparate_wopoints_linetypes.jpeg", height=6, width=8)
ggsave("Output/GrowthModelOwynChen_SexesSeparate_wopoints_linetypesnocolor.jpeg", height=6, width=8)
ggsave("Output/GrowthModelOwynChen_SexesSeparate_wopoints_linetypescolor.jpeg", height=6, width=8)