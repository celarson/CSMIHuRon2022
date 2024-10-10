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

#64 net zooplankton
Zoop64countforR <- read_csv("Zoop64countforR.csv")

#zooplankton net calibration information to calculate density
CalibrationMeta <- read_csv("CalibrationMeta.csv")
CM <- (CalibrationMeta)

#water chemistry
WaterChemistry<- read_csv("CSMIHuron2_2.csv")

#Cleaning datasets

#melt datasets in order to calculate density
Zoopcount153long<-melt(ZoopforR6_18_24, value.name="Count", 
                       variable.name = "Species")
#Merging our metadata with longform species data
Zoopcount153Mer<-merge(Zoopcount153long,CM, by="SiteID")

Zoopcount64long<-melt(Zoop64countforR, value.name = "Count", variable.name = "Species")
Zoopcount64Mer<-merge(Zoopcount64long, CM, by="SiteID")

#Create one long total zooplankton dataset
Zoopcount<-rbind(Zoopcount153Mer,Zoopcount64Mer)

#create density variables
#creating column density by taking count/volume columns
Zoopcount153Mer$Density153<-Zoopcount153Mer$Count/Zoopcount153Mer$Volume
Zoopcount64Mer$Density64<-Zoopcount64Mer$Count/Zoopcount64Mer$Volume
Zoopcount$Density<-Zoopcount$Count/Zoopcount$Volume

#go back to wide using density
Zoop64denswide<-reshape(Zoopcount64Mer, 
                        idvar = c ("SiteID","Site","Date","DFS","Month","Area"), 
                        timevar="Species",direction = "wide", 
                        drop=c("Count","MeterEnd","MeterMCoe","Volume"))

#go back to wide using density 153 net
Zoop153denswide<-reshape(Zoopcount153Mer, 
                        idvar = c ("SiteID","Site","Date","DFS","Month","Area"), 
                        timevar="Species",direction = "wide", 
                        drop=c("Count","MeterEnd","MeterMCoe","Volume"))

#find difference between two datasets
setdiff(paste(Zoop153denswide$Site,Zoop153denswide$Date),
              paste(Zoop64denswide$Site,Zoop64denswide$Date))
#Parry Sound 18 8/9/2022 only collected for 153

#merge to make one large dataset with 64 net, 153 net, and water chemistry
HuronCSMIWidezoop<-merge(Zoop64denswide, Zoop153denswide, 
                     by=c("Site","Date","DFS","Month","Area"),all=T)

setdiff(paste(HuronCSMIWidezoop$Site,HuronCSMIWidezoop$Date),
        paste((subset(WaterChemistry, Depth=="Epi"))$Site,(subset(WaterChemistry, Depth=="Epi"))$Date))

HuronCSMIWide<-merge(HuronCSMIWidezoop, subset(WaterChemistry, Depth=="Epi"), 
                         by=c("Site","Date","DFS","Area"))


#Analysis

#Cor - include all numerical variables
cors<-cor(HuronCSMIWide[,7:ncol(HuronCSMIWide)])

#Visualizations
#Make dataframe for faceted dreissena and swf density graph
