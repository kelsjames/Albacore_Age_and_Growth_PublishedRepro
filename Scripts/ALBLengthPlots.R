#Length frequency plots
rm(list = ls())

####Data wrangling
library(openxlsx)
library(tidyverse)
library(ggplot2)

#Load data
alb <- read.xlsx("Data/alb.xlsx")

#Sample collection by year and region - Figure 2 aging dataset
ggplot(alb, aes(x=Collection.Year,fill = Collection.Region))+
  geom_histogram(stat="count",boundary=0, color = "black")+
  scale_fill_manual(values=c("sienna1","lightskyblue"))+
  scale_x_continuous(breaks=seq(2013,2025,by=1))+
  scale_y_continuous(limits=c(0,140),breaks=seq(0,140,20))+
  labs(x="Year",y="Count",fill="Region")+
  theme_classic(base_size = 14)

ggsave("Output/SampleCountbyYearRegion.jpeg", height=6, width=8)


#Length frequency by Region - Figure 5
ggplot(alb, aes(x=FL,fill=Collection.Region))+
  geom_histogram(binwidth=10, boundary = 0, color = "black")+
  scale_fill_manual(values=c("sienna1","lightskyblue"))+
  labs(y="Count", x = "Fork length (cm)", fill="Region")+
  theme_classic(base_size = 14)

ggsave("Output/LengthFreqbyRegion.jpeg", height=6, width=8)

#Age frequency by Region - Figure 6
ggplot(alb, aes(x=as.factor(Final.ID1),fill=Collection.Region))+
  geom_histogram(stat="count", color = "black")+
  scale_fill_manual(values=c("sienna1","lightskyblue"))+
  labs(y="Count", x = "Age", fill="Region")+
  theme_classic(base_size = 14)

ggsave("Output/agefrequency.jpeg", height=6, width=8)

#Length Frequency by Sex with UNKNOWNS Figure 7
ggplot(alb, aes(x=FL,fill=Sex))+
  geom_histogram(binwidth=10, boundary = 0, color = "black")+
  scale_fill_manual(values=c("red","blue","grey"), labels=c("Female","Male","Unknown"))+
  facet_wrap(vars(Collection.Region),nrow=2)+
  labs(y="Count", x = "Fork length (cm)", fill="Sex")+
  theme_classic(base_size = 14)

ggsave("Output/LengthFreqbySexUnk.jpeg", height=12, width=8)

#Age Frequency by Sex with Unkonwns - Figure 8
ggplot(alb, aes(as.factor(x=Final.ID1),fill=Sex))+
  geom_histogram(stat="count", color = "black")+
  scale_fill_manual(values=c("red","blue","grey"), labels=c("Female","Male","Unknown"))+
  facet_wrap(vars(Collection.Region),nrow=2)+
  labs(y="Count", x = "Age", fill="Sex")+
  theme_classic(base_size = 14)

ggsave("Output/agefrequencybySexUnk.jpeg", height=12, width=8)

####Mean length and age for whole dataset and by sex and region####
#Length by Region for text
summary(alb)
alb %>%
  group_by(Collection.Region) %>%
  summarise_at(vars(FL), list(Min = min, Mean = mean, Max = max, SD = sd))
table(alb$Collection.Region)

epo <- subset(alb, Collection.Region == "EPO")
SEepo <- sd(epo$FL)/sqrt(233)
cpo <- subset(alb, Collection.Region == "CPO")
SEcpo <- sd(cpo$FL)/sqrt(214)
SEall <- sd(alb$FL)/sqrt(447)

#Age by Region for text
alb %>%
  group_by(Collection.Region) %>%
  summarise_at(vars(Final.ID1), list(Min = min, Mean = mean, Max = max, SD = sd))
table(alb$Collection.Region)

SEageepo <- sd(epo$Final.ID1)/sqrt(233)
SEagecpo <- sd(cpo$Final.ID1)/sqrt(214)
SEageall <- sd(alb$Final.ID1)/sqrt(447)

alb %>%
  group_by(Sex) %>%
  summarise_at(vars(FL), list(Min = min, Mean = mean, Max = max, SD = sd))
table(alb$Sex)

alb %>%
  group_by(Sex) %>%
  summarise_at(vars(Final.ID1), list(Min = min, Mean = mean, Max = max, SD = sd))

albf <- subset(alb, Sex == "Female")
SEf <- sd(albf$Final.ID1)/sqrt(109)
albm <- subset(alb, Sex == "Male")
SEm <- sd(albm$Final.ID1)/sqrt(203)

####Length at age comparisons####
#Pairwise tests
table(alb$Sex,alb$Final.ID1,useNA="always")

fem <- subset(alb, Sex == "Female")
mal <- subset(alb, Sex == "Male")

#Go through the next 6 lines of code for each age: 1-13 (no females at 14 or 15)
f <- subset(fem, Final.ID1 == 1)
m <- subset(mal, Final.ID1 == 1)

shapiro.test(f$FL)
shapiro.test(m$FL)
t.test(f$FL, m$FL)
#If non parametric required
wilcox.test(f$FL,m$FL)

#For Table 8
options(pillar.sigfig=6)
meanf <- fem %>% 
  group_by(as.factor(Final.ID1)) %>%
  summarise(
    mean = mean(FL, na.rm = TRUE),
    sd   = sd(FL, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n)
  )
meanf

meanm <- mal %>% 
  group_by(as.factor(Final.ID1)) %>%
  summarise(
    mean = mean(FL, na.rm = TRUE),
    sd   = sd(FL, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n)
  )
meanm