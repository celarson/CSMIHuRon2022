#Packages and Library loading
install.packages("tidyverse")
install.packages("psych")
install.packages("viridis")
install.packages("reshape2")
install.packages("xlsx")
install.packages("writexl")
install.packages("vegan")
install.packages("ggpubr")
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

#Change to long form zooplankton

#CM is our meta data with STIS numbers
CM <- read_csv("~/CSMI/CalibrationMeta.csv")
View(CM)

Zoop153count <- read_csv("~/CSMI/ZoopforR6.18.24.csv")
View(ZoopforR6_18_24)

Zoopcount153<-melt(Zoop153count, value.name="Count", 
                   variable.name = "Species")

#Merging our metadata with longform species data
Zoopcount153Mer<-merge(Zoopcount153,CM, by="SiteID")
#creating column density by taking count/volume columns
Zoopcount153Mer$Density<-Zoopcount153Mer$Count/Zoopcount153Mer$Volume


Zoop64count <- read_csv("~/CSMI/Zoop64countforR.csv")
View(Zoop64count)

Zoopcount64<-melt(Zoop64count, value.name = "Count", variable.name = "Species")
Zoopcount64Mer<-merge(Zoopcount64, CM, by="SiteID")
Zoopcount64Mer$Density<-Zoopcount64Mer$Count/Zoopcount64Mer$Volume

Zoop153biomass <- read_csv("~/CSMI/Zoop153biomassforR.csv", col_types = cols(Epischuralacustris = col_double()))
View(Zoop153biomass)

Zoopbiomass153<-melt(Zoop153biomass, value.name = "Count", variable.name = "Species")
Zoopbiomass153Mer<-merge(Zoopbiomass153, CM, by="SiteID")
Zoopbiomass153Mer$Biomass<-Zoopbiomass153Mer$Count/Zoopbiomass153Mer$Volume

Zoop64biomass <- read_csv("~/CSMI/Zoop64biomassforR.csv")
View(Zoop64biomass)

Zoopbiomass64<-melt(Zoop64biomass, value.name = "Count", variable.name = "Species")
Zoopbiomass64Mer<-merge(Zoopbiomass64, CM, by="SiteID")
Zoopbiomass64Mer$Biomass<-Zoopbiomass64Mer$Count/Zoopbiomass64Mer$Volume

#aggregate
Biomass153<-aggregate(Biomass ~ Area + Species, data = Zoopbiomass153Mer, FUN = mean)

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

#read data
CSMIHuron2 <- read_csv("~/CSMI/CSMIHuron2.csv")
View(CSMIHuron2)

CCTD_H_2022 <- read_csv("CSMI/2022 CSMI LH combined ctd data binned 1m depths 2.csv")
View(CCTD_H_2022)


#NH4 ug N/L BW - June, nearshore 18, mid 46 and offshore 66 82 91

CSMI4<-CSMIHuron2
CSMI4$Area<-factor(CSMI4$Area, c("NC", "SB", "GB", "SMB", "NMB"))
CSMI4$Month<-factor(CSMI4$Month, c("June", "August"))
CSMI4$Depth<-factor(CSMI4$Depth, c("Epi", "Mid", "Bottom"))
CSMI4$DFS<-factor(CSMI4$DFS, c("Nearshore", "Midshore", "Offshore"))

NH4AD<-ggplot(CSMI4, aes(x=Area, y=`NH4 ug N/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NH4 (μg N/L)")

ggplot(CSMI4, aes(x=Area, y=`NH4 ug N/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NH4 (μg N/L)")

ggplot(CSMI4, aes(x=Area, y=`NH4 ug N/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NH4 (μg N/L)")

nh4ave<-aggregate(`NH4 ug N/L` ~ Area + Depth, data = CSMI4, FUN = mean)
nh4aveM<-aggregate(`NH4 ug N/L` ~ Area + Month, data = CSMI4, FUN = mean)

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

#NOx ug N/L BW - June, nearshore 18, mid 46 and offshore 66 82 91

noxad<-aggregate(`NOx ug N/L` ~ Area + Depth, data = CSMI4, FUN = mean)

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

ggplot(CSMI4, aes(x=Area, y=`NOx ug N/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NOx (μg N/L)")

ggplot(CSMI4, aes(x=Area, y=`NOx ug N/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("NOx (μg N/L)")

ggarrange(NH4AD, NH4AM, NOxAD, NOxAM)
ggarrange(NH4AD, NH4AM)

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

ggplot(CSMI4, aes(x=Area, y=`SRP ug P/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("SRP (μg P/L)")

ggplot(CSMI4, aes(x=Area, y=`SRP ug P/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
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

ggplot(CSMI4, aes(x=Area, y=`K mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("K+ (mg/L)")

ggplot(CSMI4, aes(x=Area, y=`K mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
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

ggplot(CSMI4, aes(x=Area, y=`Na mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Na+ (mg/L)")

ggplot(CSMI4, aes(x=Area, y=`Na mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
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

ggplot(CSMI4, aes(x=Area, y=`Ca mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
  ylab("Ca++ (mg/L)")

ggplot(CSMI4, aes(x=Area, y=`Ca mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_classic()+
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

ggplot(CSMI4, aes(x=Area, y=`Mg mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Mg++ (mg/L)")

ggplot(CSMI4, aes(x=Area, y=`Mg mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
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

ggplot(CSMI4, aes(x=Area, y=`Cl mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Cl- (mg/L)")

ggplot(CSMI4, aes(x=Area, y=`Cl mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
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

ggplot(CSMI4, aes(x=Area, y=`SO4 mg/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("SO4 (mg/L)")

ggplot(CSMI4, aes(x=Area, y=`SO4 mg/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
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

ggplot(CSMI4, aes(x=Area, y=`DOC mg C/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("DOC (mg C/L)")

ggplot(CSMI4, aes(x=Area, y=`DOC mg C/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
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

ggplot(CSMI4, aes(x=Area, y=`Si mg SiO2/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Si (mg SiO2/L)")

ggplot(CSMI4, aes(x=Area, y=`Si mg SiO2/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("Si (mg SiO2/L)")

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

ggplot(CSMI4, aes(x=Area, y=`chla    ug/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("chl-a (μg/L)")

ggplot(CSMI4%>%filter(!is.na(Month)), aes(x=Area, y=`chla    ug/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("chl-a (μg/L)")


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

ggplot(CSMI4, aes(x=Area, y=`TN ug N/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TN (μg N/L)")

ggplot(CSMI4%>%filter(!is.na(Month)), aes(x=Area, y=`TN ug N/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TN (μg N/L)")


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

ggplot(CSMI4, aes(x=Area, y=`TP ug P/L`, fill = Depth))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TP (μg P/L)")

ggplot(CSMI4%>%filter(!is.na(Month)), aes(x=Area, y=`TP ug P/L`, fill = Month))+
  geom_boxplot()+
  scale_fill_brewer()+
  theme_bw()+
  ylab("TP (μg P/L)")

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
Zoo153countarea <- read_csv("CSMI/Zoo153countarea.csv")
View(Zoo153countarea)

zoo153countarea2<-Zoo153countarea
zoo153countarea2$Area<-factor(zoo153countarea2$Area, c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(zoo153countarea2, aes(x=Area, y=Group, size = `Average Density (Count/m^3)`, color = `Average Density (Count/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#153 biomass area

Zoo153biomassArea <- read_csv("CSMI/Zoo153biomassArea.csv")
View(Zoo153biomassArea)

z153bioarea2<-Zoo153biomassArea
z153bioarea2$Area<-factor(z153bioarea2$Area,c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(z153bioarea2, aes(x=Area, y=Group, size = `Biomass (ug/m^3)`, color = `Biomass (ug/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#64 biomass area

Zoo64bioarea <- read_csv("CSMI/Zoo64bioarea.csv")
View(Zoo64bioarea)

zoo64bioarea2<-Zoo64bioarea
zoo64bioarea2$Area<-factor(zoo64bioarea2$Area, c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(zoo64bioarea2, aes(x=Area, y=Group, size = `Average Biomass (ug/m^3)`, color = `Average Biomass (ug/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

#64 count area

Zoo64countarea <- read_csv("CSMI/Zoo64countarea.csv")
View(Zoo64countarea)

zoo64counarea2<-Zoo64countarea
zoo64counarea2$Month<-factor(zoo64counarea2$Month, c("NC", "SB", "GB", "SMB", "NMB"))

ggplot(zoo64counarea2, aes(x=Month, y=Group, size = `Average Density (Count/m^3)`, color = `Average Density (Count/m^3)`))+
  geom_point()+
  theme_bw()+
  scale_color_viridis()+
  scale_y_discrete(limits=rev)

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
