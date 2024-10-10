#Libraries
library(ggplot2)
library(ggpubr)
library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(lubridate)

#Upload datasets

#153 net zooplankton
ZoopforR6_18_24 <- read_csv("ZoopforR6.18.24.csv")
View(ZoopforR6_18_24)

#64 net zooplankton
Zoop64countforR <- read_csv("Zoop64countforR.csv")
View(Zoop64countforR)

#zooplankton net calibration information to calculate density
CalibrationMeta <- read_csv("CalibrationMeta.csv")
CM <- (CalibrationMeta)
View(CalibrationMeta)

#water chemistry
WaterChemistry<- read_csv("CSMIHuron2_2.csv")

#Cleaning datasets

#melt datasets in order to calculate density
Zoopcount153long<-melt(ZoopforR6_18_24, value.name="Count", 
                       variable.name = "Species")
#Merging our metadata with longform species data
Zoopcount153Mer<-merge(Zoopcount153long,CM, by="SiteID")
#creating column density by taking count/volume columns
Zoopcount153Mer$Density153<-Zoopcount153Mer$Count/Zoopcount153Mer$Volume

Zoopcount64long<-melt(Zoop64countforR, value.name = "Count", variable.name = "Species")
Zoopcount64Mer<-merge(Zoopcount64long, CM, by="SiteID")
Zoopcount64Mer$Density64<-Zoopcount64Mer$Count/Zoopcount64Mer$Volume

#go back to wide using density
Zoop64denswide<-reshape(Zoopcount64Mer, 
                        idvar = c ("SiteID","Station","Collection Date","DFS","Month","Area"), 
                        timevar="Species",direction = "wide", 
                        drop=c("Count","MeterEnd","MeterMCoe","Volume"))

#go back to wide using density 153 net
Zoop153denswide<-reshape(Zoopcount153Mer, 
                        idvar = c ("SiteID","Station","Collection Date","DFS","Month","Area"), 
                        timevar="Species",direction = "wide", 
                        drop=c("Count","MeterEnd","MeterMCoe","Volume"))

#merge to make one large dataset with 64 net, 153 net, and water chemistry
HuronCSMIWide<-cbind(Zoop64denswide, Zoop153denswide, by=c("SiteID","Station","Collection Date","DFS","Month","Area"))

#Analysis

#Cor
zoop64cors<-cor(Zoop64denswide[,7:ncol(Zoop64denswide)])
