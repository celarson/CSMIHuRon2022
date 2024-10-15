#Libraries
library(ggplot2)
library(ggpubr)
library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(lubridate)
library(readr)

#Upload datasets

#153 net zooplankton
ZoopforR6_18_24 <- read_csv("ZoopforR6.18.24.csv")

#64 net zooplankton
Zoop64countforR <- read_csv("Zoop64countforR.csv")

#zooplankton net calibration information to calculate density
CalibrationMeta <- read_csv("CalibrationMeta.csv")

#water chemistry
WaterChemistry<- read_csv("CSMIHuron2_2.csv")

#Cleaning datasets

#melt datasets in order to calculate density
Zoopcount153long<-melt(ZoopforR6_18_24, value.name="Count", 
                       variable.name = "Species")
#Merging our metadata with longform species data
Zoopcount153Mer<-merge(Zoopcount153long,CalibrationMeta, by="SiteID")

Zoopcount64long<-melt(Zoop64countforR, value.name = "Count", variable.name = "Species")
Zoopcount64Mer<-merge(Zoopcount64long, CalibrationMeta, by="SiteID")

#Create one long total zooplankton dataset
Zoopcount<-rbind(Zoopcount153Mer,Zoopcount64Mer)

#create density variables
#creating column density by taking count/volume columns
Zoopcount153Mer$Density153<-Zoopcount153Mer$Count/Zoopcount153Mer$Volume
Zoopcount64Mer$Density64<-Zoopcount64Mer$Count/Zoopcount64Mer$Volume
Zoopcount$Density<-Zoopcount$Count/Zoopcount$Volume

#go back to wide using density
Zoop64denswide<-reshape(Zoopcount64Mer, 
                        idvar = c ("SiteID","Site","UnifiedDate","DFS","Month","Area"), 
                        timevar="Species",direction = "wide", 
                        drop=c("Count","MeterEnd","MeterMCoe","Volume"))

#go back to wide using density 153 net
Zoop153denswide<-reshape(Zoopcount153Mer, 
                        idvar = c ("SiteID","Site","UnifiedDate","DFS","Month","Area"), 
                        timevar="Species",direction = "wide", 
                        drop=c("Count","MeterEnd","MeterMCoe","Volume"))

#find difference between two datasets
setdiff(paste(Zoop153denswide$Site,Zoop153denswide$UnifiedDate),
              paste(Zoop64denswide$Site,Zoop64denswide$UnifiedDate))
#Parry Sound 18 8/9/2022 only collected for 153

#merge to make one large dataset with 64 net, 153 net, and water chemistry
HuronCSMIWidezoop<-merge(Zoop64denswide, Zoop153denswide, 
                     by=c("Site","UnifiedDate","DFS","Month","Area"),all=T)

setdiff(paste(HuronCSMIWidezoop$Site,HuronCSMIWidezoop$UnifiedDate),
        paste((subset(WaterChemistry, Depth=="Epi"))$Site,
              (subset(WaterChemistry, Depth=="Epi"))$UnifiedDate))

HuronCSMIWide<-merge(HuronCSMIWidezoop, subset(WaterChemistry, Depth=="Epi"), 
                         by=c("Site","UnifiedDate","DFS","Area"))


#Analysis

#Cor - include all numerical variables
cors<-cor(HuronCSMIWide[,7:ncol(HuronCSMIWide)])

#Visualizations
#Make dataframe for faceted dreissena and swf density graph

#######Time-series line graphs
#Uploading Genus key to cut down on how many panels - binned species into genus
GenusKey <- read_csv("GenusKey.csv")

#subset data
BythoPrey2<-Zoopcount %>% filter(Species %in% c("Leptodiaptomusminutus", "Leptodiaptomussicilis", "Daphniaspp", "Daphniaparvula", "Daphniaretrocurva", "Leptodiaptomusashlandi", "Daphnialongiremis", "Leptodorakindti", "Daphniagaleatamendotae", "Diaphanosomaspp", "Bythotrepheslongimanus"))
#aggregate
BythoPreyagg<-aggregate(`Density` ~ Species + Area + Month, data = BythoPrey2, FUN = mean)

#Merge bythopreyagg and genus key
BythoPrey22<-merge(BythoPrey2,GenusKey, by="Species")

#aggregate again
BythoPreyMEAN<-aggregate(`Density` ~ Genus + Area + Month, data = BythoPrey22, FUN = mean)
BythoPreySD<-aggregate(`Density` ~ Genus + Area + Month, data = BythoPrey22, FUN = sd)
BythoPreyMEANSD2<-BythoPrey22 %>% group_by(Genus,Area,Month) %>% summarise(mean = mean(Density), std = sd(Density))


#factor

#Subset genus and area so graph separetly
Bytho<-BythoPreyMEANSD2 %>% filter(Genus %in% "Bythotrephes longimanus")
Daphnia<-BythoPreyMEANSD2 %>% filter(Genus %in% "Daphnia")
Diaphanosoma<-BythoPreyMEANSD2 %>% filter(Genus %in% "Diaphanosoma")
Leptodiaptomus<-BythoPreyMEANSD2 %>% filter(Genus %in% "Leptodiaptomus")
Leptodora<-BythoPreyMEANSD2 %>% filter(Genus %in% "Leptodora")

#factor regions
Bytho$Area<-factor(Bytho$Area, c("NC", "SB", "GB", "SMB", "NMB"))
Bytho$Month<-factor(Bytho$Month, c("June", "July", "Aug"))
Daphnia$Area<-factor(Daphnia$Area, c("NC", "SB", "GB", "SMB", "NMB"))
Daphnia$Month<-factor(Daphnia$Month, c("June", "July", "Aug"))
Diaphanosoma$Area<-factor(Diaphanosoma$Area, c("NC", "SB", "GB", "SMB", "NMB"))
Diaphanosoma$Month<-factor(Diaphanosoma$Month, c("June", "July", "Aug"))
Leptodiaptomus$Area<-factor(Leptodiaptomus$Area, c("NC", "SB", "GB", "SMB", "NMB"))
Leptodiaptomus$Month<-factor(Leptodiaptomus$Month, c("June", "July", "Aug"))
Leptodora$Area<-factor(Leptodora$Area, c("NC", "SB", "GB", "SMB", "NMB"))
Leptodora$Month<-factor(Leptodora$Month, c("June", "July", "Aug"))

#Box and Whisker graph, bytho and prey

ggplot(Bytho, aes(x=Month, y=mean, group=Month))+
  geom_boxplot()+geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
    facet_grid(.~Area)

ggplot(Daphnia, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std, color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area, scales = "free")

ggplot(Diaphanosoma, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std, color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area, scales = "free")

ggplot(Leptodiaptomus, aes(x=Month, y=Density, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area, scales = "free")

ggplot(Leptodora, aes(x=Month, y=Density, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area, scales = "free")
