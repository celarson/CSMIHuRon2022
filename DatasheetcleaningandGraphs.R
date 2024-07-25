#Visualizations for Lake Huron 2022 CSMI Report
#Noah Grode and Courtney Larson

#Set working directory to project directory
#all files are updated in github repository CSMIHuRon2022

################################################
#Packages and Library loading
library(ggtext)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(psych)
library(viridis)
library(readxl)
library(reshape2)
library(vegan)
library(ggpubr)
library(AICcmodavg)
library(rstatix)
library(emmeans)
########################################################
#Upload datasets

#Water chemistry
CSMIHuron2 <- read.csv("CSMIHuron2.csv", header=T)

#Water chemistry binned biweekly months
CSMIHuron2_2 <- read_csv("CSMIHuron2_2.csv")

#CTE
CCTD_H_2022 <- read.csv("2022CSMILHcombinedctddatabinned1mdepths2.csv", header=T)

#CM is our net tow meta data with STIS numbers
CM <- read.csv("CalibrationMeta.csv", header=T)

#2017 Water chem load
CSMI17WQ <- read_csv("CSMIHuron17WQ.csv")
View(CSMIHuron17WQ)

CSMI17meta <- read_csv("2017WQmeta.csv")
View(CSMI17meta)

#merge 2017 with meta
CSMI17<- merge(CSMI17WQ, CSMI17meta, by = "STIS")

#2017 and 2022 Meta and datasheet load
yearmeta <- read_csv("yearmeta.csv", col_types = cols(Year = col_character()))
C2217 <- read_csv("C2217.csv")

CSMIYALL <- merge(C2217, yearmeta, by = "STIS")
CSMIY <-merge(C2217, yearmeta, by = "STIS")

CSMIYALL$Month<-factor(CSMIYALL$Month, c("April", "May", "Early June", "June", "Late June", "Early July", "July", "Late July", "Aug"))

#Zoop153count are the zooplankton counts for the 153 net
Zoop153countwide<- read.csv("ZoopforR6.18.24.csv", header=T)
#melt from matrix to long form
Zoopcount153long<-melt(Zoop153countwide, value.name="Count", 
                   variable.name = "Species")
#Merging our metadata with longform species data
Zoopcount153Mer<-merge(Zoopcount153long,CM, by="SiteID")
#creating column density by taking count/volume columns
Zoopcount153Mer$Density<-Zoopcount153Mer$Count/Zoopcount153Mer$Volume

#Read in zooplankton meta with groups
Zoometa <- read_csv("Zoometa.csv")

Zoopcount153MerG<-merge(Zoopcount153Mer, Zoometa, by="Species")

#zooplankton counts for the 64 net
Zoop64countwide<- read.csv("Zoop64countforR.csv", header=T)
Zoopcount64long<-melt(Zoop64countwide, value.name = "Count", variable.name = "Species")
Zoopcount64Mer<-merge(Zoopcount64long, CM, by="SiteID")
Zoopcount64Mer$Density<-Zoopcount64Mer$Count/Zoopcount64Mer$Volume

#zooplankton biomass for the 153 net
Zoop153biomasswide<- read.csv("Zoop153biomassforR.csv", header=T)
Zoopbiomass153long<-melt(Zoop153biomasswide, value.name = "Count", variable.name = "Species", id.vars = "SiteID")
Zoopbiomass153Mer<-merge(Zoopbiomass153long, CM, by="SiteID")
Zoopbiomass153Mer$Biomass<-Zoopbiomass153Mer$Count/Zoopbiomass153Mer$Volume

#zooplankton biobass for the 64 net
Zoop64biomasswide<- read.csv("Zoop64biomassforR.csv", header=T)
Zoopbiomass64long<-melt(Zoop64biomasswide, value.name = "Count", variable.name = "Species")
Zoopbiomass64Mer<-merge(Zoopbiomass64long, CM, by="SiteID")
Zoopbiomass64Mer$Biomass<-Zoopbiomass64Mer$Count/Zoopbiomass64Mer$Volume

#Subset 153 for only bythrophes
Bytho153count<-subset(Zoopcount153Mer, Species == "Bythotrepheslongimanus")
Dresinnid64count<-subset(Zoopcount64Mer, Species == "DreissenidVeligers")

#order factors bytho
Bytho153count$Month<-factor(Bytho153count$Month, c("June", "July", "Aug") )
Bytho153count$Area<-factor(Bytho153count$Area, c("NC", "SB", "GB", "NMB", "SMB"))

#order factors dressinid
Dresinnid64count$Month<-factor(Dresinnid64count$Month, c("June", "July", "Aug"))
Dresinnid64count$Area<-factor(Dresinnid64count$Area, c("NC", "SB", "GB", "NMB", "SMB"))

BythoDresinnidcount<-rbind(Bytho153count, Dresinnid64count)

#Boxplot of invasives (Bythotrephes and Dreissenid)

ggplot(Bytho153count, aes(x=Month, y=Density, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  labs(y=expression(paste(italic("Bythotrephes longimanus"), "Density (count/m^3)")))

ggplot(Dresinnid64count, aes(x=Month, y=Density, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","springgreen3","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  scale_y_log10()+
  labs(y=expression(paste(italic("Dreissenid"), "veliger Density (count/m^3)")))
  

options(scipen = 999)
#############################################################
#Analysis

#Water chemistry
#NH4 ug N/L BW - June, nearshore 18, mid 46 and offshore 66 82 91

CSMIY$Area<-factor(CSMIY$Area, c("NC", "SB", "GB", "SMB", "NMB"))
CSMIY$Month<-factor(CSMIY$Month, c("June", "Aug"))

CSMI4<-CSMIHuron2
#order basin factor
CSMI4$Area<-factor(CSMI4$Area, c("NC", "SB", "GB", "SMB", "NMB"))
#order month factor
CSMI4$Month<-factor(CSMI4$Month, c("June", "August"))
#order depth factor
CSMI4$Depth<-factor(CSMI4$Depth, c("Epi", "Mid", "Bottom"))
#order distance from shore factor
CSMI4$DFS<-factor(CSMI4$DFS, c("Nearshore", "Midshore", "Offshore"))

CSMI5<-CSMIHuron2_2
#order basin factor
CSMI5$Area<-factor(CSMI5$Area, c("NC", "SB", "GB", "SMB", "NMB"))
#order month factor
CSMI5$Month<-factor(CSMI5$Month, c("E June", "L June", "E July", "L July", "Aug"))
#order depth factor
CSMI5$Depth<-factor(CSMI5$Depth, c("Epi", "Mid", "Bottom"))
#order distance from shore factor
CSMI5$DFS<-factor(CSMI5$DFS, c("Nearshore", "Midshore", "Offshore"))

#subset month is NA
CSMI4nomonthNA<-subset(CSMI4, Month!="NA")
#168 observations

#subset 2017 month 
CSMI17nomonthNA<-subset(CSMI17, Month!="NA")

#NH4 - no clear difference in 2017
ggplot(CSMI4nomonthNA, aes(x=Month, y=NH4ugNL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NH4 (μg N/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4nomonthNA, aes(x=Month, y=NH4ugNL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("NH4 (μg N/L)")

ggplot(CSMI17nomonthNA, aes(x=Month, y=NH4, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("NH4 (μg N/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=NH4, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("NH4 (μg N/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=NH4, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("NH4 (μg N/L)")

ggplot(CSMIYALL, aes(x= Year, y=NH4, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("NH4 (μg N/L)")

CSMIYALL$Area<-factor(CSMIYALL$Area, c("NC", "SB", "GB", "SMB", "NMB"))

CSMIYnomonth<-subset(CSMIY, Month!="NA")
CSMIYnomonth$Year<-factor(CSMIYnomonth$Year)

ggplot(CSMI4, aes(x=Area, y=NH4ugNL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NH4 (μg N/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=NH4ugNL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("NH4 (μg N/L)")+
  xlab("Region")

####################stopped here
####################################

#Aggregate
nh4ave<-aggregate(`NH4 ug N/L` ~ Area + Depth, data = CSMI4, FUN = mean)
nh4aveM<-aggregate(`NH4 ug N/L` ~ Area + Month, data = CSMI4, FUN = mean)
noxad<-aggregate(`NOx ug N/L` ~ Area + Depth, data = CSMI4, FUN = mean)
noxam<-aggregate(`NOx ug N/L` ~ Area + Month, data = CSMI4, FUN = mean)
SRPdepthsd<-aggregate(`SRP ug P/L`~Area + Depth, data = CSMI4, FUN = sd)
SRPdepthave<-aggregate(`SRP ug P/L`~Area + Depth, data = CSMI4, FUN = mean)
SRPmonthave<-aggregate(`SRP ug P/L`~Area + Month, data = CSMI4, FUN = mean)
KADepthave<- aggregate(`K mg/L`~ Area + Depth, data = CSMI4, FUN = mean)
KAMave<-aggregate(`K mg/L`~ Area + Month, data = CSMI4, FUN=mean)
naADave<- aggregate(`Na mg/L`~ Area + Depth, data = CSMI4, FUN = mean)
naAMave<- aggregate(`Na mg/L`~ Area + Month, data = CSMI4, FUN = mean)
CaADave<-aggregate(`Ca mg/L` ~ Area + Depth, data = CSMI4, FUN = mean)
CaAMave<-aggregate(`CamgL` ~ Area + Month, data = CSMI4, FUN = mean)
clADave<-aggregate(`ClmgL` ~ Area + Depth, data = CSMI4, FUN = mean)
clAMave<-aggregate(`ClmgL` ~ Area + Month, data = CSMI4, FUN = mean)
soAMave<-aggregate(`SO4mgL` ~ Area + Month, data = CSMI4, FUN = mean)
siR<-aggregate(`SimgSiO2L` ~ Area, data = CSMI4, FUN = mean)
tpnmonth<-aggregate(`TNugNL` ~ Area, data = CSMI5, FUN = mean)
zoo153ave<-aggregate(`Density` ~ Area, data = )

#plots for water chem
NH4AM<-ggplot(CSMI4, aes(x=Area, y=`NH4 ug N/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NH4 (μg N/L)")

NH4<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`NH4 ug N/L`))+
  geom_boxplot() +
  scale_fill_brewer()+
  theme_bw() +
  ylab("NH4 (μg N/L)")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4, aes(x=Area, y=`NH4 ug N/L`, fill=Area))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('North Channel', 'Saginaw Bay', 'Rest of Lake'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("NH4 (μg N/L)")+
  guides(fill=guide_legend(title=NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`NH4 ug N/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("NH4 (μg N/L)")+
  guides(fill=guide_legend(title=NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`NH4 ug N/L`, fill=Area))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("NH4 (μg N/L)")+
  guides(fill=guide_legend(title=NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`NH4 ug N/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("NH4 (μg N/L)")+
  guides(fill=guide_legend(title=NULL))

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`NH4 ug N/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("NH4 (μg N/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

#aggregate
Biomass153<-aggregate(Biomass ~ Area + Species, data = Zoopbiomass153Mer, FUN = mean)

###############work on this later, start with water chem first###############
Biomass153[,3][Biomass153[,3]==0]<-NA

ggplot(Biomass153, aes(x=Area, y=Species, size = `Biomass`, color = `Biomass`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()

Biomass153<-ggplot(Biomass153, aes(x=Area, y=Species, size = `Biomass`, color = `Biomass`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis(limits=c(.1,8000), breaks=seq(.1,8000, by=2500))+
  guides(color=guide_legend(), size=guide_legend())

Biomass153+scale_size_continuous(limits=c(.1,8000), breaks=seq(.1,8000, by=2500))

#Diversity
richness_zoo153C<-estimateR(Zoopcount153Mer$Density)
Shann_zoo153C<-diversity(Zoopcount153Mer, index = "shannon")





#NOx ug N/L BW - June, nearshore 18, mid 46 and offshore 66 82 91

ggplot(CSMI4%>%filter(!is.na(Month)), aes(x=Month, y=`NOx ug N/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("NOx (μg N/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

NOx+scale_x_discrete(limits=c("June", "August"))

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`NOx ug N/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("NOx (μg N/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`NOx ug N/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("NOx (μg N/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`NOx ug N/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NOx (μg N/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`NOx ug N/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("NOx (μg N/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=NOxugNL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NOx (ug N/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4, aes(x=Area, y=NOxugNL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("NOx (ug N/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=NOxugNL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("NOx (μg N/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=NOx, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("NOx (μg N/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=NOx, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("NOx (μg N/L)")

#SRP ug P/L BW - June, nearshore 18, mid 46 and offshore 66 82 91

SRP<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`SRP ug P/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("SRP (μg P/L)")+
  guides(fill=guide_legend(title = NULL))

SRP+scale_x_discrete(limits=c("June", "August"))

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`SRP ug P/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("SRP (μg P/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))


ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`SRP ug P/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("SRP (μg P/L)")+
  guides(fill=guide_legend(title=NULL))

ggplot(CSMI4nomonthNA, aes(x=Month, y=SRPugPL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("SRP (μg P/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4, aes(x=Area, y=SRPugPL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("SRP (μg P/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=SRPugPL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("SRP (μg P/L)")

#used in final version

ggplot(CSMI4, aes(x=Area, y=`SRP ug P/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("SRP (μg P/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`SRP ug P/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("SRP (μg P/L)")+
  xlab("Region")

ggplot(CSMIYnomonth, aes(x= Year, y=SRP, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("SRP (μg P/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=SRP, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("SRP (μg P/L)")

#K+ mg/l bw - June,  nearshore 18, mid 46 and offshore 66 82 91

K<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`K mg/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("K+ (mg/L)")+
  guides(fill=guide_legend(title = NULL))

K+scale_x_discrete(limits=c("June", "August"))

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`K mg/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("K+ (mg/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`K mg/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("K+ (mg/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`K mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("K+ (mg/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`K mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("K+ (mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=KmgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("K+ (mg/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4, aes(x=Area, y=KmgL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("K+ mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=KmgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("K+ (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=K, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("K+ (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=K, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("K+ (mg/L)")

#Na mg/L bw- June,  nearshore 18, mid 46 and offshore 66 82 91

na<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`Na mg/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("Na+ (mg/L)")+
  guides(fill=guide_legend(title = NULL))

na+scale_x_discrete(limits=c("June", "August"))                     

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`Na mg/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Na+ (mg/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`Na mg/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Na+ (mg/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`Na mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Na+ (mg/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`Na mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("Na+ (mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=NamgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Na+ (mg/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4, aes(x=Area, y=NamgL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("Na+ mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=NamgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Na+ (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=Na, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Na+ (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=Na, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Na+ (mg/L)")

#Ca mg/L bw june,nearshore 18, mid 46 and offshore 66 82 91

Ca<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`Ca mg/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("Ca++ (mg/L)")+
  guides(fill=guide_legend(title = NULL))

Ca+scale_x_discrete(limits=c("June", "August"))                     

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`Ca mg/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Ca++ (mg/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`Ca mg/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Ca++ (mg/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`Ca mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Ca++ (mg/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`Ca mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("Ca++ (mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=CamgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Ca++ (mg/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4, aes(x=Area, y=CamgL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("Ca++ mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=CamgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Ca++ (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=Ca, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Ca++ (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=Ca, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Ca++ (mg/L)")

#Mg mg/L bw  june,nearshore 18, mid 46 and offshore 66 82 91

Mg<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`Mg mg/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("Mg++ (mg/L)")+
  guides(fill=guide_legend(title = NULL))

Mg+scale_x_discrete(limits=c("June", "August")) 

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`Mg mg/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Mg++ (mg/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`Mg mg/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Mg++ (mg/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`Mg mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Mg++ (mg/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`Mg mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("Mg++ (mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=MgmgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Mg++ (mg/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4, aes(x=Area, y=MgmgL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("Mg++ mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=MgmgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Mg++ (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=Mg, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Mg++ (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=Mg, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Mg++ (mg/L)")

#Cl mg/L bw june, ,nearshore 18, mid 46 and offshore 66 82 91

Cl<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`Cl mg/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("Cl- (mg/L)")+
  guides(fill=guide_legend(title = NULL))

Cl+scale_x_discrete(limits=c("June", "August")) 

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`Cl mg/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Cl- (mg/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`Cl mg/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Cl- (mg/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`Cl mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Cl- (mg/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`Cl mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("Cl- (mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=ClmgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Cl- (mg/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4, aes(x=Area, y=ClmgL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("Cl- mg/L)")+
  xlab("Region")


ggplot(CSMI4nomonthNA, aes(x=Month, y=ClmgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Cl- (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=Cl, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Cl- (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=Cl, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Cl- (mg/L)")

#SO4 mg/L bw June, ,nearshore 18, mid 46 and offshore 66 82 91

SO4<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`SO4 mg/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("SO4 (mg/L)")+
  guides(fill=guide_legend(title = NULL))

SO4+scale_x_discrete(limits=c("June", "August")) 

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`SO4 mg/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("SO4 (mg/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`SO4 mg/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("SO4 (mg/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`SO4 mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("SO4 (mg/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`SO4 mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("SO4 (mg/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=SO4mgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("SO4 (mg/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4nomonthNA, aes(x=Month, y=SO4mgL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("SO4 (mg/L)")


ggplot(CSMI4, aes(x=Area, y=SO4mgL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("SO4 mg/L)")+
  xlab("Region")

ggplot(CSMIYnomonth, aes(x= Year, y=SO4, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("SO4 (mg/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=SO4, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("SO4 (mg/L)")

#DOC mg/Cl bw june ,nearshore 18, mid 46 and offshore 66 82 91

DOC<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`DOC mg C/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("DOC (mg C/L)")+
  guides(fill=guide_legend(title = NULL))

DOC+scale_x_discrete(limits=c("June", "August")) 

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`DOC mg C/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("DOC (mg C/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`DOC mg C/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("DOC (mg C/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`DOC mg C/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("DOC (mg C/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`DOC mg C/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("DOC (mg C/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=DOCmgCL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("DOC (mg C/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4nomonthNA, aes(x=Month, y=DOCmgCL, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("DOC (mg C/L)")

ggplot(CSMI4, aes(x=Area, y=DOCmgCL, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("DOC mgCL)")+
  xlab("Region")

ggplot(CSMIYnomonth, aes(x= Year, y=DOC, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("DOC (mg C/L)")

ggplot(CSMIYnomonth, aes(x= Year, y=DOC, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("DOC (mg C/L)")

#Si bw june,nearshore 18, mid 46 and offshore 66 82 91

Si<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`Si mg SiO2/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("Si (mg SiO2/L)")+
  guides(fill=guide_legend(title = NULL))

Si+scale_x_discrete(limits=c("June", "August")) 

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`Si mg SiO2/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Si (mg SiO2/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`Si mg SiO2/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Si (mg SiO2/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`Si mg SiO2/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Si (mg SiO2/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`Si mg SiO2/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("Si (mg SiO2/L)")+
  xlab("Region")

ggplot(CSMI4nomonthNA, aes(x=Month, y=SimgSiO2L, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Si (mg SiO2/L)")+
  facet_wrap(.~Area)

ggplot(CSMI4nomonthNA, aes(x=Month, y=SimgSiO2L, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  ylab("Si (mg SiO2/L)")

ggplot(CSMI4, aes(x=Area, y=SimgSiO2L, fill = Month))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme_classic()+
  ylab("Si (mg SiO2/L)")+
  xlab("Region")

#chla june august nearshore 18, mid 46 and offshore 66 82 91

chla<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`chla    ug/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("chl-a (μg/L)")+
  guides(fill=guide_legend(title = NULL))

chla+scale_x_discrete(limits=c("Early June","June","Late June","Early July","July", "Late July", "August"))

chla<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`chla    ug/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("chl-a (μg/L)")+
  guides(fill=guide_legend(title = NULL))

chla+scale_x_discrete(limits=c("Early June","Late June","Early July", "Late July")) 

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`chla    ug/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("chl-a (μg/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`chla    ug/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("chl-a (μg/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`chla    ug/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("chl-a (μg/L)")+
  xlab("Region")

ggplot(CSMI4%>%filter(!is.na(Month)), aes(x=Area, y=`chla    ug/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("chl-a (μg/L)")+
  xlab("Region")

ggplot(CSMI4, aes(x=Area, y=`chla    ug/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Reds 3")+
  theme_classic()+
  ylab("chl-a (μg/L)")+
  xlab("Region")

#subset month is NA
CSMI5nomonthNA<-subset(CSMI5, Month!="NA")
CSMI5nomid<-subset(CSMIlon)

CSMI5long<-melt(CSMI5)

CSMI5nomid<-subset(CSMI5long, Depth!="Mid")

CSMI5nomid2<-subset(CSMI5, Depth!="Mid")

ggplot(subset(CSMI5nomid,variable==c("chlaugL","TNugNL","TPugPL")),
       aes(x=Month, y=value, fill = Depth))+
  geom_boxplot()+
  scale_color_brewer()+
  xlab("Region")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))+
  facet_grid(variable~Area, scales = "free")

ggplot(subset(CSMI5nomid,variable=="chlaugL")+
       aes(x=Month, y=value, fill = Depth))+
  geom_boxplot()+
  scale_color_brewer()+
  xlab("Region")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=16),legend.text = element_text(size=14))+
  facet_grid(.~Area, scales = "free")

ggplot(CSMI5nomid2, aes(x=Month, y=chlaugL, fill = Depth))+
  geom_boxplot()+
  scale_color_brewer()+
  theme_classic()+
  ylab("Chl-a (μg/L)")+
  xlab("Region")+
  facet_wrap(.~Area, scales = "free_y")

#tn june august  nearshore 18, mid 46 and offshore 66 82 91

tn<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`TN ug N/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("TN (μg N/L)")+
  guides(fill=guide_legend(title = NULL))

tn+scale_x_discrete(limits=c("Early June","June","Late June","Early July","July", "Late July", "August"))

tn+scale_x_discrete(limits=c("Early June","Late June","Early July", "Late July")) 

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`TN ug N/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TN (μg N/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`TN ug N/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TN (μg N/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`TN ug N/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("TN (μg N/L)")+
  xlab("Region")

ggplot(CSMI4%>%filter(!is.na(Month)), aes(x=Area, y=`TN ug N/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TN (μg N/L)")

ggplot(CSMI5nomid2, aes(x=Month, y=TNugNL, fill = Depth))+
  geom_boxplot()+
  scale_color_brewer()+
  theme_classic()+
  ylab("TN (μg N/L)")+
  xlab("Region")+
  facet_wrap(.~Area, scales = "free_y")

#TP june august nearshore 18, mid 46 and offshore 66 82 91

tp<-ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=Month, y=`TP ug P/L`, fill = DFS))+
  geom_boxplot() +
  scale_fill_discrete(breaks=c('Nearshore', 'Midshore', 'Offshore'))+
  scale_fill_brewer()+
  theme_bw() +
  ylab("TP (μg P/L)")+
  guides(fill=guide_legend(title = NULL))

tp+scale_x_discrete(limits=c("Early June","June","Late June","Early July","July", "Late July", "August"))

tp+scale_x_discrete(limits=c("Early June","Late June","Early July", "Late July")) 

ggplot(CSMI4%>%filter(!is.na(DFS)), aes(x=DFS, y=`TP ug P/L`, fill=DFS))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TP (μg P/L)") +
  xlab("Distance from Shore")+
  guides(fill=guide_legend(title = NULL))

ggplot(CSMI4%>%filter(!is.na(Depth)), aes(x=Depth, y=`TP ug P/L`, fill=Depth))+
  geom_boxplot()+
  scale_fill_discrete(breaks=c('Epi', 'Mid', 'Bottom'))+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TP (μg P/L)")+
  guides(fill=guide_legend(title=NULL))

#used in final version

ggplot(CSMI4, aes(x=Area, y=`TP ug P/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer("Reds 3")+
  theme_classic()+
  ylab("TP (μg P/L)")+
  xlab("Region")

ggplot(CSMI4%>%filter(!is.na(Month)), aes(x=Area, y=`TP ug P/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TP (μg P/L)")

ggplot(CSMI5nomid2, aes(x=Month, y=TPugPL, fill = Depth))+
  geom_boxplot()+
  scale_color_brewer()+
  theme_classic()+
  ylab("TP (μg P/L)")+
  xlab("Region")+
  facet_wrap(.~Area, scales = "free_y")

#summary stats
#NH4 summary stats and CI
describeBy(CSMI3$`NH4 ug N/L`, CSMI3$DFS)
tNHnear<-qt(0.95, df=13-1)
19.86+c(-1,1)*tNHnear*2.49
19.86-15.42
tNHmid<-qt(0.95, df=102-1)
16.22+c(-1,1)*tNHmid*0.61
16.22-15.20
tNHoff<-qt(0.95, df=45-1)
20.72+c(-1,1)*tNHoff*2.31
20.72-16.83

describeBy(CSMI3$`NH4 ug N/L`, CSMI3$Month)

#NOx summary stats and CI
describeBy(CSMI3$`NOx ug N/L`, CSMI3$DFS)
259.03+c(-1,1)*tNHnear*18.1
259.02-226.77
283.95+c(-1,1)*tNHmid*3.77
283.95-277.69
301.11+c(-1,1)*tNHoff*4.94
301.11-292.81

#SRP summary stats and CI
describeBy(CSMI3$`SRP ug P/L`, CSMI3$DFS)
3.1+c(-1,1)*tNHnear*0.21
3.1-2.72
2.95+c(-1,1)*tNHmid*0.07
2.95-2.83
2.91+c(-1,1)*tNHoff*0.1
2.91-2.74

#K summary stats and CI
describeBy(CSMI3$`K mg/L`, CSMI3$DFS)
0.84+c(-1,1)*tNHnear*0.04
0.84-0.77
0.85+c(-1,1)*tNHmid*0.01
0.85-0.83
0.88+c(-1,1)*tNHoff*0.01
0.88-0.86

#na summary stats and CI
describeBy(CSMI3$`Na mg/L`, CSMI3$DFS)
4.33+c(-1,1)*tNHnear*0.17
4.33-4.03
4.21+c(-1,1)*tNHmid*0.04
4.21-4.14
4.24+c(-1,1)*tNHoff*0.05
4.24-4.15

#ca summary stats and CI
describeBy(CSMI3$`Ca mg/L`, CSMI3$DFS)
20.67+c(-1,1)*tNHnear*1.1
20.67-18.71
21.75+c(-1,1)*tNHmid*0.2
21.75-21.42
22.33+c(-1,1)*tNHoff*0.23
22.33-21.94

#mg summary stats and CI
describeBy(CSMI3$`Mg mg/L`, CSMI3$DFS)
6.3+c(-1,1)*tNHnear*0.4
6.3-5.6
6.6+c(-1,1)*tNHmid*0.1
6.6-6.4
6.9+c(-1,1)*tNHoff*0.09
6.9-6.7

#cl summary stats and CI
describeBy(CSMI3$`Cl mg/L`, CSMI3$DFS)
7.17+c(-1,1)*tNHnear*0.4
7.17-6.46
6.96+c(-1,1)*tNHmid*0.12
6.96-6.76
7.24+c(-1,1)*tNHoff*0.12
7.24-7.04

#DOC summary stats and CI
describeBy(CSMI3$`DOC mg C/L`, CSMI3$DFS)
2.58+c(-1,1)*tNHnear*0.16
2.58-2.29
2.14+c(-1,1)*tNHmid*0.03
2.14-2.09
2.04+c(-1,1)*tNHoff*0.03
2.04-1.99

#si summary stats and CI
describeBy(CSMI3$`Si mg SiO2/L`, CSMI3$DFS)
4.5+c(-1,1)*tNHnear*0.41
4.5-3.8
4.5+c(-1,1)*tNHmid*0.06
4.72+c(-1,1)*tNHoff*0.07
4.72-4.60

#chla summary stats and CI
describeBy(CSMI3$`chla    ug/L`, CSMI3$DFS)
tchlamid<-qt(0.95, df=155-1)
tchlaoff<-qt(0.95, df=72-1)
1.16+c(-1,1)*tNHnear*0.18
1.16-0.84
0.79+c(-1,1)*tchlamid*0.03
0.79-0.74
0.71+c(-1,1)*tchlaoff*0.04
0.71-0.64

#tn summary stats and CI
describeBy(CSMI3$`TN ug N/L`, CSMI3$DFS)
tTNmid<-qt(0.95, df=156-1)
407.46+c(-1,1)*tNHnear*19.02
407.46-373.56
419.94+c(-1,1)*tTNmid*4.17
419.94-413.03
416.85+c(-1,1)*tchlaoff*3.48
416.85-411.05

#tp summary stats and CI
describeBy(CSMI3$`TP ug P/L`, CSMI3$DFS)
6.4+c(-1,1)*tNHnear*0.4
6.4-5.7
5.4+c(-1,1)*tTNmid*0.14
5.4-5.2
5.6+c(-1,1)*tchlaoff*0.3
5.6-5.1

#so4 summary stats and CI
describeBy(CSMI3$`SO4 mg/L`, CSMI3$DFS)
13.91+c(-1,1)*tNHnear*0.72
13.91-12.63
13.61+c(-1,1)*tNHmid*0.2
13.61-13.28
14.19+c(-1,1)*tNHoff*0.2
14.19-13.

describeBy(CSMI3$`NH4 ug N/L`, CSMI3$Month)
taugNH4<-qt(0.95, df=79-1)
24.03+c(-1,1)*taugNH4*1.34
24.03-21.79
tjuneNH4<-qt(0.95, df=89-1)
12.62+c(-1,1)*tjuneNH4*0.39
12.62-11.97

describeBy(CSMI3$`NOx ug N/L`, CSMI3$Month)
282.44+c(-1,1)*tjuneNH4*4.16
282.44-275.52
287.92+c(-1,1)*taugNH4*6.17
287.92-277.65

describeBy(CSMI3$`SRP ug P/L`, CSMI3$Month)
2.93+c(-1,1)*taugNH4*0.07
2.93-2.81
2.93+c(-1,1)*tjuneNH4*0.08

describeBy(CSMI3$`K mg/L`, CSMI3$Month)
0.9+c(-1,1)*taugNH4*0.01
0.84+c(-1,1)*tjuneNH4*0.01

describeBy(CSMI3$`Na mg/L`, CSMI3$Month)
4.45+c(-1,1)*taugNH4*0.12
4.45-4.25
4.24+c(-1,1)*tjuneNH4*0.04
4.24-4.17

describeBy(CSMI3$`Ca mg/L`, CSMI3$Month)
22.51+c(-1,1)*taugNH4*0.24
22.51-22.11
21.59+c(-1,1)*tjuneNH4*0.25
21.59-21.17

describeBy(CSMI3$`Mg mg/L`, CSMI3$Month)
6.93+c(-1,1)*taugNH4*0.11
6.93-6.74
6.57+c(-1,1)*tjuneNH4*0.12
6.57-6.37

describeBy(CSMI3$`Cl mg/L`, CSMI3$Month)
7.02+c(-1,1)*tjuneNH4*0.13
7.02-6.80
7.59+c(-1,1)*taugNH4*0.23
7.59-7.21

describeBy(CSMI3$`DOC mg C/L`, CSMI3$Month)
taugDOC<-qt(0.95, df=78-1)
2.28+c(-1,1)*taugDOC*0.05
2.28-2.20
2.21+c(-1,1)*tjuneNH4*0.04
2.21-2.14

describeBy(CSMI3$`Si mg SiO2/L`, CSMI3$Month)
4.52+c(-1,1)*taugNH4*0.08
4.52-4.38
4.63+c(-1,1)*tjuneNH4*0.08
4.63-4.49

describeBy(CSMI3$`SO4 mg/L`, CSMI3$Month)
13.82+c(-1,1)*taugNH4*0.19
13.82-13.50
14.03+c(-1,1)*tjuneNH4*0.23
14.03-13.65

#bubble plots

Zoo643<-Zoo64
Zoo643$Month<-factor(Zoo643$Month, c("June", "July", "August"))


Zoo64M <- read_csv("CSMI/Zoo64MonthCode.csv")
View(Zoo64M)

Zoo3<-Zoo2
Zoo3$Month<-factor(Zoo3$Month, c("June","July", "August"))

Zoo64M2<-Zoo64M
Zoo64M2$Month<-factor(Zoo64M2$Month, c("June","July", "August"))

Zoo64DFS <- read_csv("CSMI/Zoo64DFS.csv")
View(Zoo64DFS)

Zoo64DFS2M<-Zoo64DFS
Zoo64DFS2M$DFS<-factor(Zoo64DFS2M$DFS, c("Nearshore", "Midshore", "Offshore"))

Zoo153M <- read_csv("CSMI/Zoo153Month.csv")
View(Zoo153M)

Zoo153M2<-Zoo153M
Zoo153M2$Month<-factor(Zoo153M2$Month, c("June","July", "August"))

Zoo153DFS <- read_excel("CSMI/Zoo153DFS.xlsx")
View(Zoo153DFS)

Zoo153DFS2<-Zoo153DFS
Zoo153DFS2$DFS<-factor(Zoo153DFS2$DFS, c("Nearshore", "Midshore", "Offshore"))


Zoo153biomassArea <- read_csv("~/CSMI/Zoo153biomassArea.csv")
Zoo153biomassArea$Area<-factor(Zoo153biomassArea$Area, c("NC", "SB", "GB", "SMB", "NMB"))

#64 biomass month
ggplot(Zoo64M2, aes(x=Month, y=Groups, size = `Average Biomass (ug/m^3)`, color = `Average Biomass (ug/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#64 biomass DFS
ggplot(Zoo64DFS2M, aes(x=DFS, y=Groups, size = `Average Biomass (ug/m^3)`, color = `Average Biomass (ug/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#153 biomass month
ggplot(Zoo153M2, aes(x=Month, y=Group, size = `Average Biomass (ug/m^3)`, color = `Average Biomass (ug/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#153 biomass dfs
ggplot(Zoo153DFS2, aes(x=DFS, y=Group, size = `Average Biomass (ug/m^3)`, color = `Average Biomass (ug/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#153 count DFS
Zoo153DensDFS <- read_csv("CSMI/Zoo153DensDFS.csv")

Zoo153DDFS2<-Zoo153DensDFS
Zoo153DDFS2$DFS<-factor(Zoo153DDFS2$DFS, c("Nearshore", "Midshore", "Offshore"))

ggplot(Zoo153DDFS2, aes(x=DFS, y=Group, size = `Average Density (Count/m^3)`, color = `Average Density (Count/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#153 density month

Zoo153DenMonth <- read_csv("CSMI/Zoo153DenMonth.csv")
View(Zoo153DenMonth)

Zoo153DenM2<-Zoo153DenMonth
Zoo153DenM2$Month<-factor(Zoo153DenM2$Month, c("June", "July", "August"))

ggplot(Zoo153DenM2, aes(x=Month, y=Group, size = `Average Density (Count/m^3)`, color = `Average Density (Count/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#153 count area
Zoo153countarea <- read_csv("~/CSMI/Zoo153countarea.csv")
View(Zoo153countarea)

zoo153countarea2<-Zoo153countarea
zoo153countarea2$Area<-factor(zoo153countarea2$Area, c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(zoo153countarea2, aes(x=Area, y=Group, size = `Average Density (Count/m^3)`, color = `Average Density (Count/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type = "viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")

#153 biomass area

Zoo153biomassArea <- read_csv("~/CSMI/Zoo153biomassArea.csv")
View(Zoo153biomassArea)

z153bioarea2<-Zoo153biomassArea
z153bioarea2$Area<-factor(z153bioarea2$Area,c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(Zoo153biomassArea, aes(x=Area, y=Group, size = `Biomass (ug/m^3)`, color = `Biomass (ug/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type = "viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")

ggplot(z153bioarea2, aes(x=Area, y=Group, size = `Biomass (ug/m^3)`)) +
  geom_point(shape = `Biomass (ug/m^3)`, colour = `Biomass (ug/m^3)`)+
  theme_classic()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#64 biomass area

Zoo64bioarea <- read_csv("~/CSMI/Zoo64bioarea.csv")
View(Zoo64bioarea)

zoo64bioarea2<-Zoo64bioarea
zoo64bioarea2$Area<-factor(zoo64bioarea2$Area, c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(zoo64bioarea2, aes(x=Area, y=Group, size = `Average Biomass (ug/m^3)`, color = `Average Biomass (ug/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type = "viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")

#64 count area

Zoo64countarea <- read_csv("~/CSMI/Zoo64countarea.csv")
View(Zoo64countarea)

zoo64counarea2<-Zoo64countarea
zoo64counarea2$Month<-factor(zoo64counarea2$Month, c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(zoo64counarea2, aes(x=Month, y=Group, size = `Average Density (Count/m^3)`, color = `Average Density (Count/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type = "viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")

#64countmonth

Zoo64countmonth <- read_csv("CSMI/Zoo64countmonth.csv")
View(Zoo64countmonth)

zoo64countmonth2<-Zoo64countmonth
zoo64countmonth2$Month<-factor(zoo64countmonth2$Month, c("June", "July", "August"))

ggplot(zoo64countmonth2, aes(x=Month, y=Group, size = `Average Density (Count/m^3)`, color = `Average Density (Count/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)


#64 count DFS
zoo64countDFS <- read_csv("CSMI/zoo64countDFS.csv")
View(zoo64countDFS)

zoo64countdfs<-zoo64countDFS
zoo64countdfs$DFS<-factor(zoo64countdfs$DFS, c("Nearshore", "Midshore", "Offshore"))

ggplot(zoo64countdfs, aes(x=DFS, y=Group, size = `Average Density (Count/m^3)`, color = `Average Density (Count/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#Zooplankton 153 species by group
#subset
Zoop153AdultCalnoid<-subset(Zoopcount153MerG, Group == "AdultCalanoid")
Zoop153AdImmCalanoid<-subset(Zoopcount153MerG, Group == c("AdultCalanoid", "ImmatureCalanoid"))


ggplot(Ag153DenAdultCalanoidNA, aes(x=Area, y=Species, size = `Density`, color = `Density`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type = "viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")+
  ylab("Adult Calanoid")

Ag153DenAdultCalanoid<-aggregate(`Density` ~ Area + Species, data = Zoop153AdultCalnoid, FUN = mean)
Ag153DenAdultCalanoidNA<-subset(Ag153DenAdultCalanoid, Density!= "0")

#immature and adult calanoid
Ag153AdImmCalanoid<-aggregate(`Density` ~ Area + Species, data = Zoop153AdImmCalanoid, FUN = mean)
Ag153AdImmCalanoidNA<-subset(Ag153AdImmCalanoid, Density!="0")

ggplot(Ag153AdImmCalanoidNA, aes(x=Area, y=Species, size = `Density`, color = `Density`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type = "viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")+
  ylab("Adult and Immature Calanoids")

#Zoo 153 species by group by month
#June
Ag153DenAdultCalanoidmonth<-aggregate(`Density` ~ Area + Species + Month, data = Zoop153AdultCalnoid, FUN = mean)
Ag153DenAdultCalanoidJuneNA<-subset(Ag153DenAdultCalanoidmonth, Density!="0")
Ag153DenAdultCalanoidJuneNA<-subset(Ag153DenAdultCalanoidJuneNA, Month=="June")

ggplot(Ag153DenAdultCalanoidJuneNA, aes(x=Area, y=Species, size = `Density`, color = `Density`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type="viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")+
  ylab("Adult Calanoid")

#July
Ag153DenAdultCalanoidJulyNA<-subset(Ag153DenAdultCalanoidmonth, Density!="0")
Ag153DenAdultCalanoidJulyNA<-subset(Ag153DenAdultCalanoidJulyNA, Month=="July")

ggplot(Ag153DenAdultCalanoidJulyNA, aes(x=Area, y=Species, size=`Density`,  color = `Density`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type="viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")+
  ylab("Adult Calanoid")

#August
Ag153DenAdultCalanoidAugNA<-subset(Ag153DenAdultCalanoidmonth, Density!="0")
Ag153DenAdultCalanoidAugNA<-subset(Ag153DenAdultCalanoidAugNA, Month=="Aug")

ggplot(Ag153DenAdultCalanoidAugNA, aes(x=Area, y=Species, size=`Density`,  color = `Density`))+
  geom_point()+
  theme_bw()+
  scale_color_continuous(guide="legend", type="viridis")+
  scale_y_discrete(limits=rev)+
  xlab("Region")+
  ylab("Adult Calanoid")

#############################Linear regression

#NH4
nhr<-lm(NH4 ~ Area + Year + Month, data = CSMIYALL)
residualsnhr<-nhr$residuals
hist(residualsnhr)
qqnorm(residualsnhr)
qqline(residualsnhr)
summary(nhr)

nhAY<-lm(NH4 ~ Area + Year, data = CSMIYALL)
nhA<-lm(NH4 ~ Area, data = CSMIYALL)
nhYM<-lm(NH4 ~ Year + Month, data = CSMIYALL)
nhAM<-lm(NH4 ~ Area + Month, data = CSMIYALL)
nhY<-lm(NH4 ~ Year, data = CSMIYALL)
nhM<-lm(NH4 ~ Month, data = CSMIYALL)

nhmod<-list(nhr, nhAY, nhA, nhYM, nhAM, nhY, nhM)
nhnames<-c('nhr', 'nhAY', 'nhA', 'nhYM', 'nhAM', 'nhY','nhM')
aictab(cand.set = nhmod, modnames = nhnames)

summary(nhYM)
nhIYM<-lm(NH4 ~ Year*Month, data = CSMIYALL)
nhIYMmod<-list(nhIYM)
nhIYMnam<-c('nhIYMmod')
aictab(cand.set = nhIYMmod, modnames = nhIYMnam)
summary(nhIYM)

#NOx
nor<-lm(NOx ~ Area + Year + Month, data = CSMIYALL)
rNOx<-nor$residuals
hist(rNOx)
qqnorm(rNOx)
qqline(rNOx)
summary(nor)

noAY<-lm(NOx ~ Area + Year, data = CSMIYALL)
noA<-lm(NOx ~ Area, data = CSMIYALL)
noYM<-lm(NOx ~ Year + Month, data = CSMIYALL)
noAM<-lm(NOx ~ Area + Month, data = CSMIYALL)
noY<-lm(NOx ~ Year, data = CSMIYALL)
noM<-lm(NOx ~ Month, data = CSMIYALL)

nomod<-list(nor,noAY,noA, noYM,noAM,noY,noM)
nonames<-c('nor','noAY', 'noA', 'noYM','noAM','noY','noM')
aictab(cand.set = nomod, modnames = nonames)

noIALL<-lm(NOx ~ (Year+Month+Area)^2, data = CSMIYALL)
noIMA<-lm(NOx ~ Year+Month*Area, data = CSMIYALL)
noIAY<-lm(NOx ~ Month + Area*Year, data = CSMIYALL)
noIYM<-lm(NOx ~ Month*Year + Area, data = CSMIYALL)

noImod<-list(noIALL, noIMA, noIAY, noIYM)
noInames<-c('noIALL', 'noIMA', 'noIAY', 'noIYM')
aictab(cand.set = noImod, modnames = noInames)

#SRP
srp<-lm(SRP ~ Area + Year + Month, data = CSMIYALL)
rsrp<-srp$residuals
hist(rsrp)
qqnorm(rsrp)
qqline(rsrp)
summary(srp)

srpAY<-lm(SRP ~ Area + Year, data = CSMIYALL)
srpA<-lm(SRP ~ Area, data = CSMIYALL)
srpYM<-lm(SRP ~ Year + Month, data = CSMIYALL)
srpAM<-lm(SRP ~ Area + Month, data = CSMIYALL)
srpY<-lm(SRP ~ Year, data = CSMIYALL)
srpM<-lm(SRP ~ Month, data = CSMIYALL)

srpmod<-list(srp, srpAY, srpA,srpYM,srpAM,srpY,srpM)
srpnam<-c('srp','srpAY', 'srpA', 'srpYM', 'srpAM', 'srpY','srpM')
aictab(cand.set = srpmod, modnames = srpnam)

srpIALL<-lm(SRP ~ (Year+Month+Area)^2, data = CSMIYALL)
srpIMA<-lm(SRP ~ Year+Month*Area, data = CSMIYALL)
srpIAY<-lm(SRP ~ Month + Area*Year, data = CSMIYALL)
srpIYM<-lm(SRP ~ Month*Year + Area, data = CSMIYALL)

srpImod<-list(srpIALL, srpIMA, srpIAY, srpIYM)
srpInames<-c('srpIALL', 'srpIMA', 'srpIAY', 'srpIYM')
aictab(cand.set = srpImod, modnames = srpInames)

#TN
tn<-lm(TN ~ Area + Year + Month, data = CSMIYALL)
rtn<-tn$residuals
hist(rtn)
qqnorm(rtn)
qqline(rtn)
summary(tn)

tnAY<-lm(TN ~ Area + Year, data = CSMIYALL)
tnA<-lm(TN ~ Area, data = CSMIYALL)
tnYM<-lm(TN ~ Year + Month, data = CSMIYALL)
tnAM<-lm(TN ~ Area + Month, data = CSMIYALL)
tnY<-lm(TN ~ Year, data = CSMIYALL)
tnM<-lm(TN ~ Month, data = CSMIYALL)

tnmod<-list(tn, tnAY, tnA, tnYM,tnAM,tnY,tnM)
tnnam<-c('tn', 'tnAY', 'tnA','tnYM','tnAM','tnY','tnM')
aictab(cand.set = tnmod, modnames = tnnam)

summary(tnAM)

tnIAM<-lm(TN ~ Area*Month, data = CSMIYALL)
tnImod<-list(tnIAM)
tnInam<-c('tnIAM')
aictab(cand.set = tnImod, modnames = tnInam)

#TP
tp<-lm(TP ~ Area + Year + Month, data = CSMIYALL)
rtp<-tp$residuals
hist(rtp)
qqnorm(rtp)
qqline(rtp)
summary(tp)

tpAY<-lm(TP ~ Area + Year, data = CSMIYALL)
tpA<-lm(TP ~ Area, data = CSMIYALL)
tpYM<-lm(TP ~ Year + Month, data = CSMIYALL)
tpAM<-lm(TP ~ Area + Month, data = CSMIYALL)
tpY<-lm(TP ~ Year, data = CSMIYALL)
tpM<-lm(TP ~ Month, data = CSMIYALL)

tpmod<-list(tp, tpAY, tpA, tpYM,tpAM,tpY,tpM)
tpnam<-c('tp', 'tpAY','tpA','tpYM','tpAM','tpY','tpM')
aictab(cand.set = tpmod, modnames = tpnam)

tpIALL<-lm(TP ~ (Year+Month+Area)^2, data = CSMIYALL)
tpIMA<-lm(TP ~ Year+Month*Area, data = CSMIYALL)
tpIAY<-lm(TP ~ Month + Area*Year, data = CSMIYALL)
tpIYM<-lm(TP ~ Month*Year + Area, data = CSMIYALL)

tpImod<-list(tpIALL, tpIMA, tpIAY, tpIYM)
tpInames<-c('tpIALL', 'tpIMA', 'tpIAY', 'tpIYM')
aictab(cand.set = tpImod, modnames = tpInames)

#K
k<-lm(K ~ Area + Year + Month, data = CSMIYALL)
rk<-k$residuals
hist(rk)
qqnorm(rk)
qqline(rk)
summary(k)

kAY<-lm(K ~ Area + Year, data = CSMIYALL)
kA<-lm(K ~ Area, data = CSMIYALL)
kYM<-lm(K ~ Year + Month, data = CSMIYALL)
kAM<-lm(K ~ Area + Month, data = CSMIYALL)
kY<-lm(K ~ Year, data = CSMIYALL)
kM<-lm(K ~ Month, data = CSMIYALL)

kmod<-list(k,kAY,kA,kYM,kAM,kY,kM)
knam<-c('k', 'kAY', 'kA','kYM','kAM','kY','kM')
aictab(cand.set = kmod, modnames = knam)

kIALL<-lm(K ~ (Year+Month+Area)^2, data = CSMIYALL)
kIMA<-lm(K ~ Year+Month*Area, data = CSMIYALL)
kIAY<-lm(K ~ Month + Area*Year, data = CSMIYALL)
kIYM<-lm(K ~ Month*Year + Area, data = CSMIYALL)

kImod<-list(kIALL, kIMA, kIAY, kIYM)
kInames<-c('kIALL', 'kIMA', 'kIAY', 'kIYM')
aictab(cand.set = kImod, modnames = kInames)

#na
na<-lm(Na ~ Area + Year + Month, data = CSMIYALL)
rna<-na$residuals
hist(rna)
qqnorm(rna)
qqline(rna)
summary(na)

naAY<-lm(Na ~ Area + Year, data = CSMIYALL)
naA<-lm(Na ~ Area, data = CSMIYALL)
naYM<-lm(Na ~ Year + Month, data = CSMIYALL)
naAM<-lm(Na ~ Area + Month, data = CSMIYALL)
naY<-lm(Na ~ Year, data = CSMIYALL)
naM<-lm(Na ~ Month, data = CSMIYALL)

namod<-list(na, naAY, naA,naYM,naAM,naY,naM)
naname<-c('na', 'naAY', 'naA','naYM','naAM','naY','naM')
aictab(cand.set = namod, modnames = naname)

naIALL<-lm(Na ~ (Year+Month+Area)^2, data = CSMIYALL)
naIMA<-lm(Na ~ Year+Month*Area, data = CSMIYALL)
naIAY<-lm(Na ~ Month + Area*Year, data = CSMIYALL)
naIYM<-lm(Na ~ Month*Year + Area, data = CSMIYALL)

naImod<-list(naIALL, naIMA, naIAY, naIYM)
naInames<-c('naIALL', 'naIMA', 'naIAY', 'naIYM')
aictab(cand.set = naImod, modnames = naInames)

#ca
ca<-lm(Ca ~ Area + Year + Month, data = CSMIYALL)
rca<-ca$residuals
hist(rca)
qqnorm(rca)
qqline(rca)
summary(ca)

caAY<-lm(Ca ~ Area + Year, data = CSMIYALL)
caA<-lm(Ca ~ Area, data = CSMIYALL)
caYM<-lm(Ca ~ Year + Month, data = CSMIYALL)
caAM<-lm(Ca ~ Area + Month, data = CSMIYALL)
caY<-lm(Ca ~ Year, data = CSMIYALL)
caM<-lm(Ca ~ Month, data = CSMIYALL)

camod<-list(ca, caAY, caA,caYM,caAM,caY,caM)
caname<-c('ca', 'caAY', 'caA','caYM','caAM','caY','caM')
aictab(cand.set = camod, modnames = caname)

caIALL<-lm(Ca ~ (Year+Month+Area)^2, data = CSMIYALL)
caIMA<-lm(Ca ~ Year+Month*Area, data = CSMIYALL)
caIAY<-lm(Ca ~ Month + Area*Year, data = CSMIYALL)
caIYM<-lm(Ca ~ Month*Year + Area, data = CSMIYALL)

caImod<-list(caIALL, caIMA, caIAY, caIYM)
caInames<-c('caIALL', 'caIMA', 'caIAY', 'caIYM')
aictab(cand.set = caImod, modnames = caInames)

#Mg
mg<-lm(Mg ~ Area + Year + Month, data = CSMIYALL)
mg2<-lm(Mg ~ Area + Year + Month + DFS, data = CSMIYALL)
rmg<-mg$residuals
hist(rmg)
qqnorm(rmg)
qqline(rmg)
summary(mg)
rmg2<-mg2$residuals
hist(rmg2)
qqnorm(rmg2)

mgAY<-lm(Mg ~ Area + Year, data = CSMIYALL)
mgA<-lm(Mg ~ Area, data = CSMIYALL)
mgYM<-lm(Mg ~ Year + Month, data = CSMIYALL)
mgAM<-lm(Mg ~ Area + Month, data = CSMIYALL)
mgY<-lm(Mg ~ Year, data = CSMIYALL)
mgM<-lm(Mg ~ Month, data = CSMIYALL)

mgmod<-list(mg,mgAY,mgA,mgYM,mgAM,mgY,mgM)
mgnam<-c('mg', 'mgAY', 'mgA','mgYM','mgAM','mgY','mgM')
aictab(cand.set = mgmod, modnames = mgnam)

mgIALL<-lm(Mg ~ (Year+Month+Area)^2, data = CSMIYALL)
mgIMA<-lm(Mg ~ Year+Month*Area, data = CSMIYALL)
mgIAY<-lm(Mg ~ Month + Area*Year, data = CSMIYALL)
mgIYM<-lm(Mg ~ Month*Year + Area, data = CSMIYALL)

mgImod<-list(mgIALL, mgIMA, mgIAY, mgIYM)
mgInames<-c('mgIALL', 'mgIMA', 'mgIAY', 'mgIYM')
aictab(cand.set = mgImod, modnames = mgInames)

#cl
cl<-lm(Cl ~ Area + Year + Month, data = CSMIYALL)
rcl<-cl$residuals
hist(rcl)
qqnorm(rcl)
qqline(rcl)
summary(cl)

clAY<-lm(Cl ~ Area + Year, data = CSMIYALL)
clA<-lm(Cl ~ Area, data = CSMIYALL)
clYM<-lm(Cl ~ Year + Month, data = CSMIYALL)
clAM<-lm(Cl ~ Area + Month, data = CSMIYALL)
clY<-lm(Cl ~ Year, data = CSMIYALL)
clM<-lm(Cl ~ Month, data = CSMIYALL)

clmod<-list(cl,clAY,clA,clYM,clAM,clY,clM)
clnam<-c('cl', 'clAY', 'clA','clYM','clAM','clY','clM')
aictab(cand.set = clmod, modnames = clnam)

summary(clAY)

clIAY<-lm(Cl ~ Area*Year, data = CSMIYALL)
clImod<-list(clIAY)
clInam<-c('clIAY')
aictab(cand.set = clImod, modnames = clInam)

#so
so<-lm(SO4 ~ Area + Year + Month, data = CSMIYALL)
rso<-so$residuals
hist(rso)
qqnorm(rso)
qqline(rso)
summary(so)

soAY<-lm(SO4 ~ Area + Year, data = CSMIYALL)
soA<-lm(SO4 ~ Area, data = CSMIYALL)
soYM<-lm(SO4 ~ Year + Month, data = CSMIYALL)
soAM<-lm(SO4 ~ Area + Month, data = CSMIYALL)
soY<-lm(SO4 ~ Year, data = CSMIYALL)
soM<-lm(SO4 ~ Month, data = CSMIYALL)

somod<-list(so,soAY,soA,soYM,soAM,soY,soM)
sonam<-c('so', 'soAY', 'soA','soYM','soAM','soY','soM')
aictab(cand.set = somod, modnames = sonam)

soIAY<-lm(SO4 ~ Area*Year, data = CSMIYALL)
soImod<-list(soIAY)
soInam<-c('soIAY')
aictab(cand.set = soImod, modnames = soInam)

#chla
chla<-lm(chla ~ Area + Year + Month, data = CSMIYALL)
rchla<-chla$residuals
hist(rchla)
qqnorm(rchla)
qqline(rchla)
summary(chla)

chAY<-lm(chla ~ Area + Year, data = CSMIYALL)
chA<-lm(chla ~ Area, data = CSMIYALL)
chYM<-lm(chla ~ Year + Month, data = CSMIYALL)
chAM<-lm(chla ~ Area + Month, data = CSMIYALL)
chY<-lm(chla ~ Year, data = CSMIYALL)
chM<-lm(chla ~ Month, data = CSMIYALL)

chmod<-list(chla,chAY,chA,chYM,chAM,chY,chM)
chnam<-c('chla', 'chAY', 'chA','chYM','chAM','chY','chM')
aictab(cand.set = chmod, modnames = chnam)

chIALL<-lm(chla ~ (Year+Month+Area)^2, data = CSMIYALL)
chIMA<-lm(chla ~ Year+Month*Area, data = CSMIYALL)
chIAY<-lm(chla ~ Month + Area*Year, data = CSMIYALL)
chIYM<-lm(chla ~ Month*Year + Area, data = CSMIYALL)

chImod<-list(chIALL, chIMA, chIAY, chIYM)
chInames<-c('chIALL', 'chIMA', 'chIAY', 'chIYM')
aictab(cand.set = chImod, modnames = chInames)

#DOC
doc<-lm(DOC ~ Area + Year + Month, data = CSMIYALL)
rdoc<-doc$residuals
hist(rdoc)
qqnorm(rdoc)
qqline(rdoc)
summary(doc)

docAY<-lm(DOC ~ Area + Year, data = CSMIYALL)
docA<-lm(DOC ~ Area, data = CSMIYALL)
docYM<-lm(DOC ~ Year + Month, data = CSMIYALL)
docAM<-lm(DOC ~ Area + Month, data = CSMIYALL)
docY<-lm(DOC ~ Year, data = CSMIYALL)
docM<-lm(DOC ~ Month, data = CSMIYALL)

docmod<-list(doc,docAY,docA,docYM,docAM,docY,docM)
docnam<-c('doc', 'docAY', 'docA','docYM','docAM','docY','docM')
aictab(cand.set = docmod, modnames = docnam)

docIALL<-lm(DOC ~ (Year+Month+Area)^2, data = CSMIYALL)
docIMA<-lm(DOC ~ Year+Month*Area, data = CSMIYALL)
docIAY<-lm(DOC ~ Month + Area*Year, data = CSMIYALL)
docIYM<-lm(DOC ~ Month*Year + Area, data = CSMIYALL)

docImod<-list(docIALL, docIMA, docIAY, docIYM)
docInames<-c('docIALL', 'docIMA', 'docIAY', 'docIYM')
aictab(cand.set = docImod, modnames = docInames)

summary(docIALL)

summary(doc)
summary(docAY)
summary(docA)
summary(docYM)
summary(AM)

################################################################

#ANOVA
aov.NH4<-CSMIYALL%>%anova_test(NH4 ~ Area*Month)

###NOx

#run anova
CSMIYALL$Year<-as.factor(CSMIYALL$Year)
CSMIYALL$Month<-as.factor(CSMIYALL$Month)

aov.NOx<-CSMIYALL%>%anova_test(NOx ~ (Year+Month+Area)^2)
#See interactions each has
CSMIYALL %>% group_by(Year) %>%
  anova_test(NOx ~ Month, error = noIALL)
CSMIYALL %>% group_by(Area) %>%
  anova_test(NOx ~ Month, error = noIALL)
#pairwise analysis
CSMIYALL %>% group_by(Year) %>% emmeans_test(NOx ~ Month, p.adjust.method = "bonferroni", model = noIALL)
CSMIYALL %>% group_by(Month) %>% emmeans_test(NOx ~ Year, p.adjust.method = "bonferroni", model = noIALL)
droplevels(CSMIYALL$Month)
CSMIYALL$Year
is.factor(CSMIYALL$Year)
levels(CSMIYALL$Year)
is.factor(CSMIYALL$Month)
levels(CSMIYALL$Month)
CSMIYALL$Month<-droplevels(CSMIYALL$Month)
CSMIYALLNEW<-subset(CSMIYALL, Month!=c("April","May", "Early June", "Late June", "Early July", "Late July"))

CSMIYnomonth<-subset(CSMIY, Month!="NA")