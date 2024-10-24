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
library(gridExtra)
library(Hmisc)

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
                        idvar = c ("SiteID","Site","UnifiedDate","DFS","Month","Area","MonthPeriod"), 
                        timevar="Species",direction = "wide", 
                        drop=c("Count","MeterEnd","MeterMCoe","Volume"))

#go back to wide using density 153 net
Zoop153denswide<-reshape(Zoopcount153Mer, 
                        idvar = c ("SiteID","Site","UnifiedDate","DFS","Month","Area","MonthPeriod"), 
                        timevar="Species",direction = "wide", 
                        drop=c("Count","MeterEnd","MeterMCoe","Volume"))

#find difference between two datasets
setdiff(paste(Zoop153denswide$Site,Zoop153denswide$UnifiedDate),
              paste(Zoop64denswide$Site,Zoop64denswide$UnifiedDate))
#Parry Sound 18 8/9/2022 only collected for 153

#merge to make one large dataset with 64 net, 153 net, and water chemistry
HuronCSMIWidezoop<-merge(Zoop64denswide, Zoop153denswide, 
                     by=c("Site","UnifiedDate","DFS","Month","Area",
                          "MonthPeriod"),all=T)

setdiff(paste(HuronCSMIWidezoop$Site,HuronCSMIWidezoop$UnifiedDate),
        paste((subset(WaterChemistry, Depth=="Epi"))$Site,
              (subset(WaterChemistry, Depth=="Epi"))$UnifiedDate))
#Parry Sound 18 8/9/2022 only collected for 153
#Thessalon River 5 6/19/2022 only zooplankton, no water chemistry
setdiff(paste((subset(WaterChemistry, Depth=="Epi"))$Site,
          (subset(WaterChemistry, Depth=="Epi"))$UnifiedDate),
          paste(HuronCSMIWidezoop$Site,HuronCSMIWidezoop$UnifiedDate))
#Thessalon River 18 6/19/2022

setdiff(paste((subset(WaterChemistry, Depth=="Epi"))$Site,(subset(WaterChemistry, Depth=="Epi"))$Month),paste(HuronCSMIWidezoop$Site, HuronCSMIWidezoop$Month))
setdiff(paste((subset(WaterChemistry, Depth=="Epi"))$Site,(subset(WaterChemistry, Depth=="Epi"))$MonthPeriod),paste(HuronCSMIWidezoop$Site, HuronCSMIWidezoop$MonthPeriod))
setdiff(paste((subset(WaterChemistry, Depth=="Epi"))$Site,(subset(WaterChemistry, Depth=="Epi"))$DFS),paste(HuronCSMIWidezoop$Site, HuronCSMIWidezoop$DFS))
setdiff(paste((subset(WaterChemistry, Depth=="Epi"))$Site,(subset(WaterChemistry, Depth=="Epi"))$Area),paste(HuronCSMIWidezoop$Site, HuronCSMIWidezoop$Area))
#No differences.

HuronCSMIWide<-merge(HuronCSMIWidezoop, subset(WaterChemistry, Depth=="Epi"), 
                         by=c("Site","UnifiedDate","DFS","Area","MonthPeriod","Month"))


#Analysis

#Cor - include all numerical variables
names(HuronCSMIWide)
summary(HuronCSMIWide)
NHuronCSMIWide<-as.numeric(unlist(HuronCSMIWide[,c(8:51,53:91,103:116)]))
cors<-rcorr((as.numeric(HuronCSMIWide[,c(8:51,53:91,103:116)])))

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
DaphniaSB<-BythoPreyMEANSD2 %>% filter(Area %in% "SB")
DaphniaSB<-DaphniaSB %>% filter(Genus %in% "Daphnia")
DaphniaALL<-BythoPreyMEANSD2 %>% filter(Area %in% c("NC", "GB", "NMB", "SMB"))
DaphniaALL<-DaphniaALL %>% filter(Genus %in% "Daphnia")

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
DaphniaSB$Month<-factor(DaphniaSB$Month, c("June", "July", "Aug"))
DaphniaALL$Month<-factor(DaphniaALL$Month, c("June", "July", "Aug"))
DaphniaALL$Area<-factor(DaphniaALL$Area, c("NC", "GB", "SMB", "NMB"))

#Line graph, bytho and prey
ggplot(Bytho, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std, color=Area))+
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
    facet_grid(.~Area)

ggplot(DaphniaSB, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  scale_color_manual(values = "violetred3")
  
ggplot(DaphniaALL, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  labs(y=expression(paste(italic("Daphnia"), "(individuals/m^3)")))+
  theme(legend.position = "none")

ggplot(Daphnia, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  
  

DaphFig<-grid.arrange(DALL, DSB, ncol=2)

DaphFig+

ggplot(Diaphanosoma, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std, color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area, scales = "free")

ggplot(Leptodiaptomus, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std, color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area, scales = "free")

ggplot(Leptodora, aes(x=Month, y=mean, group=Area)) +
  geom_line(aes(color=Area))+
  geom_point(aes(color=Area))+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std, color=Area))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area, scales = "free")

######box and whisker plots

#subset for bw
BWBytho <- BythoPrey22 %>% filter(Genus %in% "Bythotrephes longimanus")
BWDaphnia <- BythoPrey22 %>% filter(Genus %in% "Daphnia")

#Factor
BWBytho$Month<-factor(BWBytho$Month, c("June", "July", "Aug"))
BWBytho$Area<-factor(BWBytho$Area, c("NC", "SB", "GB", "SMB", "NMB"))
BWDaphnia$Month<-factor(BWDaphnia$Month, c("June", "July", "Aug"))
BWDaphnia$Area<-factor(BWDaphnia$Area, c("NC", "SB", "GB", "SMB", "NMB"))
BWBytho$DFS<-factor(BWBytho$DFS, c("Nearshore", "Midshore", "Offshore"))
#######BWs
ggplot(BWBytho, aes(x=Month, y=Density, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  labs(y=expression(paste(italic("Bythotrephes cederstromii "), "(individuals/m^3)")))
#Subset to Remove SB from Bytho
BWBythoNOSB<-BWBytho %>% filter(Area %in% c("NC", "GB", "NMB", "SMB"))

ggplot(BWBythoNOSB, aes(x=Month, y=Density, fill = DFS))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  labs(y=expression(paste(italic("Bythotrephes cederstromii "), "(individuals/m^3)")))

ggplot(BWBythoNOSB, aes(x=DFS, y=Density, fill = DFS))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  labs(y=expression(paste(italic("Bythotrephes cederstromii "), "(individuals/m^3)")))


ggplot(BWDaphnia, aes(x=Month, y=Density, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  scale_y_log10()
 
ggplot(BWDaphnia, aes(x=Month, y=Density, fill = Month))+
   geom_boxplot()+
   scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
   theme(panel.background = element_rect(fill = "white", colour = "grey50"),
         axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
         axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
         legend.title=element_text(size=14),legend.text = element_text(size=14))+
   facet_grid(.~Area)

#subset out SB to use in own graph?
SBPrey<-BythoPrey22%>%filter(Area %in% "SB")
SBDaphPrey<-SBPrey%>%filter(Genus %in% "Daphnia")
SBDaphPreyRemoveO<-SBDaphPrey%>%filter(Species %in% c("Daphniagaleatamendotae", "Daphniaretrocurva"))
ALLR<-BythoPrey22%>%filter(Area %in% c("NC", "GB", "NMB", "SMB"))
ALLRDaph<-ALLR%>%filter(Genus %in% "Daphnia")

#factor
SBDaphPreyRemoveO$Month<-factor(SBDaphPreyRemoveO$Month, c("June", "July", "Aug"))

#BW for just SB, removed 0 valued daphnia to help graph
ggplot(SBDaphPreyRemoveO, aes(x=Month, y=Density, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme_classic()

ggplot(SBDaphPrey, aes(x=Month, y=Density, fill=Month))+
  geom_boxplot()+
  geom_jitter()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme_classic()+
  scale_y_log10()
  

#Bw all regions for daphnia
ggplot(ALLRDaph, aes(x=Month, y=Density, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)

#Dreisnnid
DreisnnidAgg<-aggregate(`Density64` ~ Area + Month, data = Dresinnid64count, FUN = mean)

ggplot(Dresinnid64count, aes(x=Month, y=Density64, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)

#remove dresinnid64 outlier -83
Dresinnidout<-Dresinnid64count[-83,]

Dresinnidout$Area<-factor(Dresinnidout$Area, c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(Dresinnidout, aes(x=Month, y=Density64, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  labs(y=expression(paste(italic("Dreissnid"), "veliger ", "(individuals/m^3)")))

DreisnnidAgg2<-aggregate(`Density64` ~ Month + Area, data = Dresinnidout, FUN = mean)
