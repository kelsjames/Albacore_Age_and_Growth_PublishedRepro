rm(list = ls())

#Alb Bias and Precision
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(FSA)

#load in data
alb.ages <- read.xlsx("Data/alb.xlsx")

#Bias
#consider final age of ID1 to be best
ab.1R21RF <- ageBias(AgeID1.R2~Final.ID1, data = alb.ages,
                   ref.lab = "Final ID1", nref.lab = "ID1 R2")

#Normal Age-bias plot
plotAB(ab.1R21RF)
plot(ab.1R21RF, col.CIsig = "black")
plot(ab.1R21RF, show.CI = TRUE, show.range=FALSE)
plot(ab.1R21RF, xvals = "mean")
#the summary gives all the stats, which I like, including bias tests
summary(ab.1R21RF)

#Compare Owyn 2 to Owyn 1
ab.1R11R2 <- ageBias(AgeID1.R1~AgeID1.R2, data = alb.ages,
                   ref.lab = "ID1 R2", nref.lab = "ID1 R1")
plotAB(ab.1R11R2)
plot(ab.1R11R2)
summary(ab.1R11R2)

#Compare ID2 to ID1, shrink the dataset to only those ID2 aged
small <- subset(alb.ages, AgeID2.R1 != "NA")
ab.2R21RF <- ageBias(AgeID2.R2~Final.ID1, data = small,
                   ref.lab = "Final ID1", nref.lab = "ID2 R2")
plotAB(ab.2R21RF)
plot(ab.2R21RF)
summary(ab.2R21RF)

#ID2 R2 to ID2 R2
ab.2R21R2 <- ageBias(AgeID2.R2~AgeID1.R2, data = small,
                   ref.lab = "ID1 R2", nref.lab = "ID2 R2")
plotAB(ab.2R21R2)
plot(ab.2R21R2)
summary(ab.2R21R2)

#ID2 R1 to ID1 R1
ab.2R11R1 <- ageBias(AgeID2.R1~AgeID1.R1, data = small,
                   ref.lab = "ID1 R1", nref.lab = "ID2 R1")
plotAB(ab.2R11R1)
plot(ab.2R11R1)
summary(ab.2R11R1)

#ID2 R1 vs R2
ab.2R12R2 <- ageBias(AgeID2.R1~AgeID2.R2, data = alb.ages,
                   ref.lab = "ID2 R2", nref.lab = "ID2 R1")
plotAB(ab.2R12R2)
plot(ab.2R12R2)
summary(ab.2R12R2)

#Precision
#ID1
ap.1 <- agePrecision(~AgeID1.R2+Final.ID1,data=alb.ages)
summary(ap.1, what = "difference")
summary(ap.1, what = "precision")

ap.1R21 <- agePrecision(~AgeID1.R1+AgeID1.R2,data=alb.ages)
summary(ap.1R21, what = "difference")
summary(ap.1R21, what = "precision")

#ID2
ap.2R21 <- agePrecision(~AgeID2.R2+AgeID2.R1,data=alb.ages)
summary(ap.2R21, what = "difference")
summary(ap.2R21, what = "precision")

#ID2 R2 to ID1 Final
ap.ID21F <- agePrecision(~AgeID2.R2+Final.ID1,data=alb.ages)
summary(ap.ID21F, what = "difference")
summary(ap.ID21F, what = "precision")

#ID2 R2 to ID1 R2
ap.ID21 <- agePrecision(~AgeID2.R2+AgeID1.R2,data=alb.ages)
summary(ap.ID21, what = "difference")
summary(ap.ID21, what = "precision")
#Can also compare more than 2 at a time (Ch 4 Age Comparisons, Ogle)

#ID2 R1 to ID1 R1
ap.ID2111 <- agePrecision(~AgeID2.R1+AgeID1.R1,data=alb.ages)
summary(ap.ID2111, what = "difference")
summary(ap.ID2111, what = "precision")
