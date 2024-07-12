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



