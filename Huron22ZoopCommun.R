##################################################libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(readr)
library(Hmisc)
library(EnvStats)
library(vegan)
library(reshape2)
library(lubridate)
library(readxl)
library(data.table)
library(leaps)
library(olsrr)
library(MASS)
library(fossil)


WQepianalysis2 <- read_excel("WQepianalysis2.xlsx")
ZoopcombforWQ <- read_excel("Zoopcomb.xlsx")

############################################################Checking magnitude of species counts between 2017 and 2022

#subset unifieddate and julian out
zoopcombmagnitude<-subset(zoopcomb, select = -c(UnifiedDate, Julian))
zoopcomblong<-melt(zoopcombmagnitude, value.name = "Density",
                     variable.name = "Species")

ggplot(zoopcombmagnitude, aes(x=Year, y=Bythotrepheslongimanus, fill = Year))+
  geom_boxplot()+
  theme_classic()

############################################################community analysis

######set up datasets to run vegan

#create key 1. comb 2. 2022 3. 2017
ZoopKey<-subset(zoopcomb, select = c("SiteID","Site","UnifiedDate","CruiseSeason","Year", "MonthPeriod", "DFS", "Month", "Area", "Julian"))
ZoopKey22<-subset(zoop22comb, select = c("SiteID","Site","UnifiedDate","CruiseSeason","Year", "MonthPeriod", "DFS", "Month", "Area", "Julian"))
ZoopKey17<-subset(zoop17comb, select = c("SiteID","Site","UnifiedDate","CruiseSeason","Year", "MonthPeriod", "DFS", "Month", "Area", "Julian"))
ZoopKeycomb2<-subset(zoopcomb2, select = c("SiteID","Site","UnifiedDate","CruiseSeason","Year", "MonthPeriod", "DFS", "Month", "Area", "Julian"))

#dataset with just density stuff 1. comb 2. 2022 3. 2017
ZoopDensitycomb<-subset(zoopcomb, select = -c(SiteID, Site,UnifiedDate, CruiseSeason,Year, MonthPeriod, DFS, Month, Area, Julian))
ZoopDensity22<-subset(zoop22comb, select = -c(SiteID, Site,UnifiedDate, CruiseSeason,Year, MonthPeriod, DFS, Month, Area, Julian))
ZoopDensity17<-subset(zoop17comb, select = -c(SiteID, Site,UnifiedDate, CruiseSeason,Year, MonthPeriod, DFS, Month, Area, Julian))
ZoopDensitycomb2<-subset(zoopcomb2, select = -c(SiteID, Site,UnifiedDate, CruiseSeason,Year, MonthPeriod, DFS, Month, Area, Julian))

#replaceNA's with 0
ZoopDensitycombnona<-ZoopDensitycomb
ZoopDensitycombnona22<-ZoopDensity22
ZoopDensitycombnona17<-ZoopDensity17
ZoopDensitycomb2nona<-ZoopDensitycomb2
ZoopDensitycombnona[is.na(ZoopDensitycombnona)]<-0
ZoopDensitycombnona22[is.na(ZoopDensitycombnona22)]<-0.000000
ZoopDensitycombnona17[is.na(ZoopDensitycombnona17)]<-0
ZoopDensitycomb2nona[is.na(ZoopDensitycomb2nona)]<-0

#datamatrix
Zoopdensitymat<-as.matrix(ZoopDensitycomb)
Zoopdensitymatnona<-as.matrix(ZoopDensitycombnona)
Zoopdensitymatnona22<-as.matrix(ZoopDensitycombnona22)
Zoopdensitymatnona17<-as.matrix(ZoopDensitycombnona17)
Zoopdensitymatnonacomb2<-as.matrix(ZoopDensitycomb2nona)

#adonis permanova
set.seed(19980222)
adonis2(Zoopdensitymatnonacomb2~(DFS+Area+Year+CruiseSeason)^2, data=ZoopKeycomb2, permutations = 999, method = "bray", na.rm = TRUE)

#NMDS
adoniscomb2<-metaMDS(Zoopdensitymatnonacomb2, distance="bray", k=2, trymax = 100, na.rm=TRUE)
adoniscomb2Jacc<-metaMDS(Zoopdensitymatnonacomb2, distance="jaccard", k=2, trymax = 100, na.rm=TRUE)

#Set up enfit for NMDS
MDA.enfitcomb2<-envfit(adoniscomb2, ZoopKeycomb2, permutations = 999)
Jacc.enfitcomb2<-envfit(adoniscomb2Jacc, ZoopKeycomb2, permutations = 999)

#set up zoopscores data
zoopscorescomb2<-as.data.frame(scores(adoniscomb2, display = "sites"))
zoopscorescomb2Jacc<-as.data.frame((scores(adoniscomb2Jacc, display = "sites")))
#bind with key for zoop scores 1. comb 2. 22 3. 17
zoopscorescomb2<-cbind(zoopscorescomb2, Area = ZoopKeycomb2$Area)
zoopscorescomb2<-cbind(zoopscorescomb2, DFS = ZoopKeycomb2$DFS)
zoopscorescomb2<-cbind(zoopscorescomb2, CruiseSeason = ZoopKeycomb2$CruiseSeason)
zoopscorescomb2<-cbind(zoopscorescomb2, Year = ZoopKeycomb2$Year)

zoopscorescomb2Jacc<-cbind(zoopscorescomb2Jacc, Area = ZoopKeycomb2$Area)
zoopscorescomb2Jacc<-cbind(zoopscorescomb2Jacc, DFS = ZoopKeycomb2$DFS)
zoopscorescomb2Jacc<-cbind(zoopscorescomb2Jacc, CruiseSeason = ZoopKeycomb2$CruiseSeason)
zoopscorescomb2Jacc<-cbind(zoopscorescomb2Jacc, Year = ZoopKeycomb2$Year)

#enfitdata
env.scorecomb<-as.data.frame(scores(MDA.enfitcomb2, display = "vectors"))
env.scorecomb<-cbind(env.scorecomb, env.variables = rownames(env.scorecomb))

jaccenv.scorecomb<-as.data.frame(scores(Jacc.enfitcomb2, display = "vectors"))
jaccenv.scorecomb<-cbind(jaccenv.scorecomb, env.variables = rownames(jaccenv.scorecomb))
#ordination plot comb
#factor levels
zoopscorescomb2$Area<-factor(zoopscorescomb2$Area, levels = c("NC", "SB", "GB", "NMB", "SMB"))
zoopscorescomb2$DFS<-factor(zoopscorescomb2$DFS, levels = c("Nearshore", "Midshore", "Offshore"))
zoopscorescomb2$CruiseSeason<-factor(zoopscorescomb2$CruiseSeason, levels = c("Late Spring", "Early Summer", "Summer", "Late Summer"))

zoopscorescomb2Jacc$Area<-factor(zoopscorescomb2Jacc$Area, levels =c("NC", "SB", "GB", "NMB", "SMB"))
zoopscorescomb2Jacc$DFS<-factor(zoopscorescomb2Jacc$DFS, levels =c("Nearshore", "Midshore", "Offshore"))
zoopscorescomb2Jacc$CruiseSeason<-factor(zoopscorescomb2Jacc$CruiseSeason, levels = c("Late Spring", "Early Summer", "Summer", "Late Summer"))

ggplot(zoopscorescomb2, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, colour = factor(zoopscorescomb2$Area)))+
  coord_fixed()+
  stat_ellipse(aes(color=Area))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Area")

ggplot(zoopscorescomb2, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, colour = factor(zoopscorescomb2$DFS)))+
  coord_fixed()+
  stat_ellipse(aes(color=DFS))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "DFS")

ggplot(zoopscorescomb2, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, colour = factor(CruiseSeason)))+
  coord_fixed()+
  stat_ellipse(aes(color=CruiseSeason))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "CruiseSeason")

ggplot(zoopscorescomb2, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, colour = factor(Year)))+
  coord_fixed()+
  stat_ellipse(aes(color=Year))+
  scale_color_manual(values = c("lightgreen", "darkgreen"))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Year")

ggplot(zoopscorescomb2Jacc, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, colour = factor(zoopscorescomb2$Area)))+
  coord_fixed()+
  stat_ellipse(aes(color=Area))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Area")

ggplot(zoopscorescomb2Jacc, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, colour = factor(zoopscorescomb2$DFS)))+
  coord_fixed()+
  stat_ellipse(aes(color=DFS))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "DFS")

ggplot(zoopscorescomb2Jacc, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, colour = factor(CruiseSeason)))+
  coord_fixed()+
  stat_ellipse(aes(color=CruiseSeason))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "CruiseSeason")

ggplot(zoopscorescomb2Jacc, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, colour = factor(Year)))+
  coord_fixed()+
  stat_ellipse(aes(color=Year))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Year")

#########species richness

#prep datasets for Richness
Zoopcombforrich <- read_excel("Zoopcombforrich.xlsx")
Zooprich<-subset(Zoopcombforrich, select = c("SiteID","Site","UnifiedDate","CruiseSeason","Year", "MonthPeriod", "DFS", "Month", "Area", "Julian", "Richness"))

#Check to see if richness is evenly distributed 
hist(Zooprich$Richness)
qqnorm(Zooprich$Richness)
qqline(Zooprich$Richness)
#Anova
sppr_aov<- aov(Zooprich$Richness ~ (DFS+Area+CruiseSeason+Year)^2, data = Zooprich)
summary(sppr_aov)

Zooprich$Area<-factor(Zooprich$Area, levels = c("NC", "SB", "GB", "NMB", "SMB"))
Zooprich$DFS<-factor(Zooprich$DFS, levels = c("Nearshore", "Midshore", "Offshore"))
Zooprich$CruiseSeason<-factor(Zooprich$CruiseSeason, levels = c("Late Spring", "Early Summer", "Summer", "Late Summer"))

ggplot(Zooprich, aes(x=Year, y=Richness, fill = Year))+
  geom_boxplot()+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  facet_grid(.~Area)+
  labs(y="Species Richness", x=NULL)

ggplot(Zooprich, aes(x=CruiseSeason, y=Richness, fill = Year))+
  geom_boxplot()+scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  labs(y="Species Richness", x=NULL)

ggplot(Zooprich, aes(x=CruiseSeason, y=Richness, fill = CruiseSeason))+
  geom_boxplot()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  labs(y="Species Richness", x=NULL)

ggplot(Zooprich, aes(x=DFS, y=Richness, fill = Year))+
  geom_boxplot()+scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  labs(y="Species Richness", x=NULL)

ggplot(Zooprich, aes(x=DFS, y=Richness, fill = DFS))+
  geom_boxplot()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  labs(y="Species Richness", x=NULL)

ggplot(Zooprich, aes(x=Area, y=Richness, fill = Year))+
  geom_boxplot()+scale_fill_manual(values=c("lightgreen","darkgreen"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  labs(y="Species Richness", x=NULL)

ggplot(Zooprich, aes(x=Area, y=Richness, fill = Area))+
  geom_boxplot()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  labs(y="Species Richness", x=NULL)

ggplot(Zooprich, aes(x=CruiseSeason, y=Richness, fill = DFS))+
  geom_boxplot()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=14),
        legend.title=element_text(size=14),legend.text = element_text(size=14))+
  labs(y="Species Richness", x=NULL)

######WQ

#import dataset of WQ data
C2217 <- read_csv("C2217.csv")
View(C2217)

#import metadata
yearmeta2 <- read_csv("yearmeta2.csv")
View(yearmeta2)

Date <- read_csv("Date20222017.csv")

#merge datasets
WQALL <- merge(C2217, yearmeta2, by = "STIS")
WQALL<-merge(WQALL, Date, by = "STIS")

#creat Julian Date
WQALL$Date<-as.Date(WQALL$Date, "%m/%d/%y")
WQALL$Julian<-yday(WQALL$Date)
head(WQALL$Julian)

##################################################Correlation test
#check to see what each column is vectored as
names(WQALL)
summary(WQALL)
str(WQALL)
#change year column to character vector
WQALL$Year<-as.character(WQALL$Year)
summary(WQALL)

#select only numeric columns?
WQALLnum<-select_if(WQALL, is.numeric)

#remove STIS, Staioncode, Date and Log variables
WQALLcor<-select(WQALL, select=-c(STIS, Date, StationCode, K_log, Ca_log, Mg_log, Na_log, Cl_log, SO4_log))

#run correlation
corsWQ<-rcorr(as.matrix(WQALLnum),type="spearman")
corsdf=data.frame(corsWQ$r)
corsdfp=data.frame(corsWQ$P)

#run correlation with categorical included
corsWQcat<-rcorr(as.matrix(WQALLcor), type="spearman")

#write CSV
write.csv(corsdf, "corsdata2WQ.csv")
write.csv(corsdfp, "corsdata2WQp.csv")

#k and Ca 0.764539
#mg and k 0.81902
#so4 and k 0.780153
#Na and Cl 0.73889
#Mg and Ca 0.911285
#SO4 and Ca 0.839194
#So4 and Mg 0.908743

#important variables
#NH4
#NOx
#TP
#SRP
#K
#Ca
#SO4

##################################################Scatterplot correlations
#K and Ca
ggplot(WQALL, aes(x=K, y=Ca))+
  geom_point()+
  theme_classic()

ggplot(WQALL, aes(x=K, y=Ca, color = Year))+
  geom_point()+
  theme_classic()

#K and Mg
ggplot(WQALL, aes(x=K, y=Mg))+
  geom_point()+
  theme_classic()

ggplot(WQALL, aes(x=K, y=Mg, color = Year))+
  geom_point()+
  theme_classic()

#k and SO4
ggplot(WQALL, aes(x=K, y=SO4))+
  geom_point()+
  theme_classic()

ggplot(WQALL, aes(x=K, y=SO4, color = Year))+
  geom_point()+
  theme_classic()

#Na and Cl
ggplot(WQALL, aes(x=Na, y=Cl))+
  geom_point()+
  theme_classic()

ggplot(WQALL, aes(x=Na, y=Cl, color = Year))+
  geom_point()+
  theme_classic()

#Mg and Ca
ggplot(WQALL, aes(x=Mg, y=Ca))+
  geom_point()+
  theme_classic()

ggplot(WQALL, aes(x=Mg, y=Ca, color = Year))+
  geom_point()+
  theme_classic()

#SO4 and Ca
ggplot(WQALL, aes(x=SO4, y=Ca))+
  geom_point()+
  theme_classic()

ggplot(WQALL, aes(x=SO4, y=Ca, color = Year))+
  geom_point()+
  theme_classic()

#SO4 and Mg
ggplot(WQALL, aes(x=SO4, y=Mg))+
  geom_point()+
  theme_classic()

ggplot(WQALL, aes(x=SO4, y=Mg, color = Year))+
  geom_point()+
  theme_classic()


##################################################Normal Distributions
###Will use a histogram, qqplot, boxplot and shapiro-wilk test to determine normality

##WQ
#####Made a mistake - need to run full dataset instead of 22 and 17 sep

##New WQ

#K
hist(WQALL$K)
qqnorm(WQALL$K)
qqline(WQALL$K)
ggplot(WQALL, aes(y=K))+
  geom_boxplot()+
  theme_classic()
shapiro.test(WQALL$K)
#0.95529 p=2.237e-10
WQALL$K_log<-log(WQALL$K)
hist(WQALL$K_log)
qqnorm(WQALL$K_log)
qqline(WQALL$K_log)
ggplot(WQALL, aes(y=K_log))+
  geom_boxplot()+
  theme_classic()
shapiro.test(WQALL$K_log)
#W=0.96 p=1.183e-09
#rosner test to detect outliers based on boxplot
#5 outliers
rosnerTest(WQALL$K, k=5)

#k on just Epi
hist(WQepianalysis2$K)
qqnorm(WQepianalysis2$K)
qqline(WQepianalysis2$K)
shapiro.test(WQepianalysis2$K)
#p=5.26e-07


#Ca
hist(WQALL$Ca)
qqnorm(WQALL$Ca)
qqline(WQALL$Ca)
ggplot(WQALL, aes(y=Ca))+
  geom_boxplot()+
  theme_classic()
shapiro.test(WQALL$Ca)
#W=0.98381 p=7.012e-05
WQALL$Ca_log<-log(WQALL$Ca)
hist(WQALL$Ca_log)
qqnorm(WQALL$Ca_log)
qqline(WQALL$Ca_log)
ggplot(WQALL, aes(y=Ca_log))+
  geom_boxplot()+
  theme_classic()
shapiro.test(WQALL$Ca_log)
#w=0.94967 p=3.522e-11

#Ca just epi
hist(WQepianalysis2$Ca)
qqnorm(WQepianalysis2$Ca)
qqline(WQepianalysis2$Ca)
shapiro.test(WQepianalysis2$Ca)
#p 0.01302

#Mg
hist(WQALL$Mg)
qqnorm(WQALL$Mg)
qqline(WQALL$Mg)
ggplot(WQALL, aes(y=Mg))+
  geom_boxplot()+
  theme_classic()
shapiro.test(WQALL$Mg)
#W=0.93557 p=5.793e-13
WQALL$Mg_log<-log(WQALL$Mg)
hist(WQALL$Mg_log)
qqnorm(WQALL$Mg_log)
qqline(WQALL$Mg_log)
ggplot(WQALL, aes(y=Mg_log))+
  geom_boxplot()+
  theme_classic()
shapiro.test(WQALL$Mg_log)
#w=0.31464 p=2.2e-16

#mg just epi
hist(WQepianalysis2$Mg)
qqnorm(WQepianalysis2$Mg)
qqline(WQepianalysis2$Mg)
shapiro.test(WQepianalysis2$Mg)
#p=0.01499

#Na
hist(WQALL$Na)
qqnorm(WQALL$Na)
qqline(WQALL$Na)
ggplot(WQALL, aes(y=Na))+
  geom_boxplot()+
  theme_classic()
shapiro.test(WQALL$Na)
#W = 0.67804, p-value < 2.2e-16
WQALL$Na_log<-log(WQALL$Na)
hist(WQALL$Na_log)
qqnorm(WQALL$Na_log)
qqline(WQALL$Na_log)
ggplot(WQALL, aes(y=Na_log))+
  geom_boxplot()+
  theme_classic()
shapiro.test(WQALL$Na_log)
#W = 0.79611, p-value < 2.2e-16

#Na just epi
hist(WQepianalysis2$Na)
qqnorm(WQepianalysis2$Na)
qqline(WQepianalysis2$Na)
shapiro.test(WQepianalysis2$Na)
#p=2.2e-16

#Cl
hist(WQALL$Cl)
qqnorm(WQALL$Cl)
qqline(WQALL$Cl)
shapiro.test(WQALL$Cl)
#W = 0.74054, p-value < 2.2e-16
WQALL$Cl_log<-log(WQALL$Cl)
hist(WQALL$Cl_log)
qqnorm(WQALL$Cl_log)
qqline(WQALL$Cl_log)
shapiro.test(WQALL$Cl_log)
#W = 0.59541, p-value < 2.2e-16

#Cl just epi
WQepianalysis2$Cl<-as.numeric(WQepianalysis2$Cl)
hist(WQepianalysis2$Cl)
qqnorm(WQepianalysis2$Cl)
qqline(WQepianalysis2$Cl)
shapiro.test(WQepianalysis2$Cl)
#5.36e-15

#SO4
hist(WQALL$SO4)
qqnorm(WQALL$SO4)
qqline(WQALL$SO4)
shapiro.test(WQALL$SO4)
#W = 0.89236, p-value < 2.2e-16
WQALL$SO4_log<-log(WQALL$SO4)
hist(WQALL$SO4)
qqnorm(WQALL$SO4_log)
qqline(WQALL$SO4_log)
shapiro.test(WQALL$SO4_log)
#W = 0.46498, p-value < 2.2e-16

#SO4 epi
hist(WQepianalysis2$SO4)
qqnorm(WQepianalysis2$SO4)
qqline(WQepianalysis2$SO4)
shapiro.test(WQepianalysis2$SO4)
#1.17e-10

#NH4 epi
hist(WQepianalysis2$NH4)
qqnorm(WQepianalysis2$NH4)
qqline(WQepianalysis2$NH4)
shapiro.test(WQepianalysis2$NH4)
#p=3.232e-12

#NOx epi
hist(WQepianalysis2$NOx)
qqnorm(WQepianalysis2$NOx)
qqline(WQepianalysis2$NOx)
shapiro.test(WQepianalysis2$NOx)
#p=0.1422

#SRP epi
hist(WQepianalysis2$SRP)
qqnorm(WQepianalysis2$SRP)
qqline(WQepianalysis2$SRP)
shapiro.test(WQepianalysis2$SRP)
#p=7.059e-05

#TN epi
hist(WQepianalysis2$TN)
qqnorm(WQepianalysis2$TN)
qqline(WQepianalysis2$TN)
shapiro.test(WQepianalysis2$TN)
#2.484e-07

#TP
hist(WQepianalysis2$TP)
qqnorm(WQepianalysis2$TP)
qqline(WQepianalysis2$TP)
shapiro.test(WQepianalysis2$TP)
#p=1.21e-10

#Log transfrom
WQepianalysis$Cl<-as.numeric(WQepianalysis$Cl)

WQepianalysis2<-WQepianalysis

WQepianalysis2$Cl_log<-log(WQepianalysis2$Cl)
WQepianalysis2$NH4_log<-log(WQepianalysis2$NH4)
WQepianalysis2$SO4_log<-log(WQepianalysis2$SO4)
WQepianalysis2$Na_log<-log(WQepianalysis2$Na)
WQepianalysis2$K_log<-log(WQepianalysis2$K)
WQepianalysis2$TP_log<-log(WQepianalysis2$TP)
WQepianalysis2$TN_log<-log(WQepianalysis2$TN)
WQepianalysis2$SRP_log<-log(WQepianalysis2$SRP)
WQepianalysis2$chla_log<-log(WQepianalysis2$chla)
WQepianalysis2$DOC_log<-log(WQepianalysis2$DOC)


#SO4 epi transformed
hist(WQepianalysis2$SO4_log)
qqnorm(WQepianalysis2$SO4_log)
qqline(WQepianalysis2$SO4_log)
shapiro.test(WQepianalysis2$SO4_log)

#NH4 epi transformed
hist(WQepianalysis2$NH4_log)
qqnorm(WQepianalysis2$NH4_log)
qqline(WQepianalysis2$NH4_log)
shapiro.test(WQepianalysis2$NH4_log)

#SRP epi transformed
hist(WQepianalysis2$SRP_log)
qqnorm(WQepianalysis2$SRP_log)
qqline(WQepianalysis2$SRP_log)
shapiro.test(WQepianalysis2$SRP_log)

#TN epi transformed
hist(WQepianalysis2$TN_log)
qqnorm(WQepianalysis2$TN_log)
qqline(WQepianalysis2$TN_log)
shapiro.test(WQepianalysis2$TN_log)

#TP epi transformed
hist(WQepianalysis2$TP_log)
qqnorm(WQepianalysis2$TP_log)
qqline(WQepianalysis2$TP_log)
shapiro.test(WQepianalysis2$TP_log)

#K epi transformed
hist(WQepianalysis2$K_log)
qqnorm(WQepianalysis2$K_log)
qqline(WQepianalysis2$K_log)
shapiro.test(WQepianalysis2$K_log)

#Na epi transformed
hist(WQepianalysis2$Na_log)
qqnorm(WQepianalysis2$Na_log)
qqline(WQepianalysis2$Na_log)
shapiro.test(WQepianalysis2$Na_log)

#Cl epi transformed
hist(WQepianalysis2$Cl_log)
qqnorm(WQepianalysis2$Cl_log)
qqline(WQepianalysis2$Cl_log)
shapiro.test(WQepianalysis2$Cl_log)

#chla
hist(WQepianalysis2$chla)
qqnorm(WQepianalysis2$chla)
qqline(WQepianalysis2$chla)
shapiro.test(WQepianalysis2$chla)

hist(WQepianalysis2$chla_log)
qqnorm(WQepianalysis2$chla_log)
qqline(WQepianalysis2$chla_log)
shapiro.test(WQepianalysis2$chla_log)

#DOC
hist(WQepianalysis2$DOC_log)
qqnorm(WQepianalysis2$DOC_log)
qqline(WQepianalysis2$DOC_log)
shapiro.test(WQepianalysis2$DOC_log)

#corr.test with log transformed

#remove STIS, Staioncode, Date and Log variables
WQALLcor2<-select(WQepianalysis2, select=-c(STIS, Date, StationCode, Depth, DFS, Month, Area, Year, Julian, NH4, SO4, Cl, Na, K, TP, TN, SRP, chla, DOC))

#run correlation
corsWQ<-rcorr(as.matrix(WQALLcor2),type="spearman")
corsdf=data.frame(corsWQ$r)
corsdfp=data.frame(corsWQ$P)

#run correlation with categorical included
corsWQcat<-rcorr(as.matrix(WQALLcor), type="spearman")

#write CSV
write.csv(corsdf, "corsdata2WQ.csv")
write.csv(corsdfp, "corsdata2WQp.csv")


#variable reduction based on biological relevance and correlation
#remove SO4_log, Ca_log, Mg_log, K_log and Na_log
WQepianalysis2<-subset(WQepianalysis2, select=-c(Depth, SO4, SO4_log, Ca, Mg,
                                                 K, K_log, Na, Na_log))

WQepianalysis2<-subset(WQepianalysis2, select=-c(NH4, TP, TN, SRP, DOC, chla, Cl))

writexl::write_xlsx(zoopcomb2, "Zoopcomb.xlsx")
ZoopcombforWQ <- read_excel("Zoopcomb.xlsx")

WQepianalysis3 <- WQepianalysis2

WQZOOP<-right_join(WQepianalysis3, ZoopcombforWQ, by = "STIS")

WQZOOP<-subset(WQZOOP, select=-c(SiSiO2, StationCode, UnifiedDate, Date, Year.y, DFS.y, Month.y, Area.y, Julian.y))
names(WQZOOP)[names(WQZOOP)=='Month.x']<-'Month'
names(WQZOOP)[names(WQZOOP)=='DFS.x']<-'DFS'
names(WQZOOP)[names(WQZOOP)=='Area.x']<-'Area'
names(WQZOOP)[names(WQZOOP)=='Year.x']<-'Year'
names(WQZOOP)[names(WQZOOP)=='Julian.x']<-'Julian'

WQZOOP<-WQZOOP[-c(200),]

#####linear regressions test with species richness and variables
#create dataset with richness and WQ

Zooprich2<-Zooprich
writexl::write_xlsx(Zooprich2, "Zooprich2.xlsx")
Zooprich2 <- read_excel("Zooprich2.xlsx")

RichWQ<-merge(Zooprich2, WQepianalysis3, by = "STIS")

#linear regression

bestfit<-regsubsets(Richness ~ NOx + Cl_log + NH4_log + TP_log + TN_log + SRP_log + chla_log + DOC_log, data = RichWQ)
summary(bestfit)

NOX.mod<-lm(Richness ~ NOx, data = RichWQ)
summary(NOX.mod)
#p<0.05

CL_log.mod<-lm(Richness ~ Cl_log, data = RichWQ)
summary(CL_log.mod)
#p>0.05

Fullmodel.mod<-lm(Richness ~ NOx + NH4_log + TP_log + TN_log + SRP_log + chla_log, data = RichWQ)
summary(Fullinteractions.mod)

all.model = ols_step_all_possible(Fullmodel.mod)
all.model

stepAIC(Fullmodel.mod, direction = "both")
stepAIC(Fullmodel.mod, direction = "backward")


ggplot(RichWQ, aes(Richness, chla_log))+
  geom_point()+
  geom_smooth(method = 'lm')

ggplot(RichWQ, aes(Richness, NH4_log))+
  geom_point()+
  geom_smooth(method = 'lm')

#Permanova with WQ data
WQKey<-subset(WQZOOP, select = c("STIS","NOx","DFS","Month","Area", "Year","Julian",
                                   "Cl_log", "NH4_log", "TP_log", "TN_log", "SRP_log", "chla_log",
                                   "DOC_log", "Site", "MonthPeriod", "CruiseSeason"))

WQMat<-subset(WQZOOP, select = -c(STIS,NOx,DFS,Month,Area,Year,Julian,
                                  Cl_log, NH4_log, TP_log, TN_log, SRP_log, chla_log,
                                  DOC_log, Site, MonthPeriod, CruiseSeason))
WQKey$DFS<-as.factor(WQKey$DFS)
WQKey$CruiseSeason<-as.factor(WQKey$CruiseSeason)
WQKey$Area<-as.factor(WQKey$Area)
WQKey$Year<-as.factor(WQKey$Year)

WQMat<-as.matrix(WQMat)
WQMatnona<-WQMat
WQMatnona[is.na(WQMatnona)]<-0
rowSums(WQMatnona)

adonis2(WQMatnona ~ (NOx + NH4_log +TP_log + SRP_log + chla_log + DFS + CruiseSeason + Area + Year)^2, data = WQKey, permutations = 999, method = "bray", na.rm = TRUE)

######################################
#################### Data manipulation 
#import dataset of WQ data
C2217 <- read_csv("C2217.csv")
View(C2217)


#import for just epi dataset
CSMIHuron2217epi <- read_excel("CSMIHuron2_22forzoop.xlsx")

#import metadata
yearmeta2 <- read_csv("yearmeta2.csv")
View(yearmeta2)

Date <- read_csv("Date20222017.csv")

#merge datasets
WQALL <- merge(C2217, yearmeta2, by = "STIS")
WQALL<-merge(WQALL, Date, by = "STIS")

WQepianalysis<-merge(CSMIHuron2217epi, yearmeta2, by = "STIS")
WQepianalysis<-merge(WQepianalysis, Date, by = "STIS")

#creat Julian Date
WQALL$Date<-as.Date(WQALL$Date, "%m/%d/%y")
WQALL$Julian<-yday(WQALL$Date)
head(WQALL$Julian)


WQepianalysis$Date<-as.Date(WQepianalysis$Date, "%m/%d/%y")
WQepianalysis$Julian<-yday(WQepianalysis$Date)
writexl::write_xlsx()


#subset epi values
WQepianalysis2<- subset(WQepianalysis, Depth == "Epi")

######zooplankton

###########2017 Zoop Do NOT run this code!/Use it

#import 2017 dataset
Zoop2017_2 <- read_csv("Zoop2017_2.csv")
#delete duplicate rows
Zoop2017dataset<-distinct(Zoop2017_2)
#calculate meter end
Zoop2017dataset$MeterEnd<-Zoop2017dataset$FLOW_END-Zoop2017dataset$FLOW_START
#Calculate volume
Zoop2017dataset$Volume<-Zoop2017dataset$MeterEnd*Zoop2017dataset$DISTANCE_PARM*0.19634954

#write out so you can delete non-zooplankton rows and edit datasheet - see metadata sheet for edit notes and add STIS numbers
writexl::write_xlsx(Zoop2017dataset, "Zoop2017dataset.xlsx")

#import 2017 zoop to reshape and melt with calibration meta 2017
CalibrationMeta2017 <- read_excel("CalibrationMeta2017.xlsx")
Zoop2017cleaned <- read_excel("Zoop2017dataset.xlsx")

#calculate Julian date
CalibrationMeta2017$UnifiedDate<-as.Date(CalibrationMeta2017$UnifiedDate, "%m/%d/%y")
CalibrationMeta2017$Julian<-yday(CalibrationMeta2017$UnifiedDate)

#merge 2017 zoop with calibration meta sheet
zoopdensity2017_2 <- merge(Zoop2017cleaned, CalibrationMeta2017, by="SiteID")

#reshape to wide format
zoopdensity2017wide<-reshape(zoopdensity2017_2, 
                             idvar = c ("SiteID","Site","UnifiedDate","CruiseSeason","Year", "DFS","Month","Area","MonthPeriod","Julian"), 
                             timevar="REF_NAME",direction = "wide", 
                             drop=c("Count","MeterEnd","LIFE_STAGE_NAME","LIFE_STAGE","ITIS","Density_calc","TOTAL_COUNT","SAMPLE_VOLUME",
                                    "ALIQUOT_VOLUME","SUB_SAMPLE","STRAT_MAX","STRAT_MIN", "TOD",
                                    "Volume", "DISTANCE_PARM", "FLOW_START", "FLOW_END",
                                    "OP_DATE", "VESSEL", "SERIAL", "PORT", "STATION_DEPTH","CRUISE"))

#write out just to look at it quick
writexl::write_xlsx(zoopdensity2017_3, "Zoopdensitywide2017.xlsx")

##############################################Run this code for zoop2017
Zoopdensitywide2017 <- read_excel("Zoopdensitywide2017.xlsx")

###2022 zoop

#import zoop153
Zoop153_22 <- read_csv("ZoopforR6.21.24.csv")
View(Zoop153_22)

#import zoop64
Zoop64_22 <- read_csv("Zoop64countforR.csv")
View(Zoop64_22)

CalibrationMeta <- read_csv("CalibrationMeta.csv")

#calculate Julian date for calibration meta
CalibrationMeta$UnifiedDate<-as.Date(CalibrationMeta$UnifiedDate, "%m/%d/%y")
CalibrationMeta$Julian<-yday(CalibrationMeta$UnifiedDate)
head(CalibrationMeta$Julian)

#melt datasets in order to calculate density
Zoopcount153long<-melt(Zoop153_22, value.name="Count", 
                       variable.name = "Species")

Zoopcount64long<-melt(Zoop64_22, value.name = "Count",
                      variable.name = "Species")

#Merging our metadata with longform species data
Zoopcount153Mer<-merge(Zoopcount153long,CalibrationMeta, by="SiteID")

Zoopcount64Mer<-merge(Zoopcount64long, CalibrationMeta, by="SiteID")

#Create one long total zooplankton dataset 2022
Zoopcount<-rbind(Zoopcount153Mer,Zoopcount64Mer)

#create density variable 
Zoopcount$Density<-Zoopcount$Count/Zoopcount$Volume

#reshape density from long to wide
Zoopdensity2022wide<-reshape(Zoopcount, 
                             idvar = c ("SiteID","Site","UnifiedDate","Year","CruiseSeason","DFS","Month","Area","MonthPeriod","Julian"), 
                             timevar="Species",direction = "wide", 
                             drop=c("Count","MeterEnd","MeterMCoe","Volume"))

writexl::write_xlsx(Zoopdensitywide, "Zoop2022densitywide.xlsx")

#####Make taxonomic rankings match between 2017 and 2022

#merge columns in 2022
#Ascomorpha ovalis and Ascomorpha eucaudis -> Ascomorpha
zoop2022<-Zoopdensity2022wide
zoop2022$Ascomorpha<-zoop2022$Density.Ascomorphaeucaudis+zoop2022$Density.Ascomorphaovalis
#Bosmina longirostris and Bosmina spp -> Bosmina
zoop2022$Bosmina<-zoop2022$Density.Bosminalongirostris+zoop2022$Density.Bosminaspp
#Keratella multiple species to -> Keratella
zoop2022$Keratella<-zoop2022$Density.Keratellacochlearis+zoop2022$Density.Keratellacrassa+zoop2022$Density.Keratellaearlinae+zoop2022$Density.Keratellahiemalis+zoop2022$Density.Keratellaquadrata
#Lecane multiple species to -> Lecane
zoop2022$Lecane<-zoop2022$Density.Lecaneflexilis+zoop2022$Density.Lecaneludwigi+zoop2022$Density.Lecanespp
#Monostlya multiple species to -> Monostyla
zoop2022$Monostyla<-zoop2022$Density.Monostylacopeis+zoop2022$Density.Monostylaspp+zoop2022$Density.Monostylalunaris
#Notholca multiple species to -> Notholca
zoop2022$Notholca<-zoop2022$Density.Notholcaacuminata+zoop2022$Density.Notholcafoliacea+zoop2022$Density.Notholcalaurentiae+zoop2022$Density.Notholcasquamula
#Ploesoma multiple species to -> Ploesoma
zoop2022$Ploesoma<-zoop2022$Density.Ploesomahudsoni+zoop2022$Density.Ploesomalenticulare+zoop2022$Density.Ploesomatruncatum
#Trichocerca multiple species to -> Trichocerca
zoop2022$Trichocerca<-zoop2022$Density.Trichocercacylindrica+zoop2022$Density.Trichocercamulticrinis+zoop2022$Density.Trichocercaporcellus+zoop2022$Density.Trichocercarousseleti+zoop2022$Density.Trichocercaspp

#rename 2022 column names
zoop2022_2<-zoop2022
names(zoop2022_2)[names(zoop2022_2)=='Density.Alonaspp']<-'Alona'
names(zoop2022_2)[names(zoop2022_2)=='Density.Asplanchnapriodonta']<-'Asplanchna'
names(zoop2022_2)[names(zoop2022_2)=='Density.Bythotrepheslongimanus']<-'Bythotrepheslongimanus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Calanoidcopepodites']<-'Calanoida'
names(zoop2022_2)[names(zoop2022_2)=='Density.Collothecaspp']<-'Collotheca'
names(zoop2022_2)[names(zoop2022_2)=='Density.Conochilusunicornis']<-'Conochilus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Copepodnauplii']<-'Copepoda'
names(zoop2022_2)[names(zoop2022_2)=='Density.Cyclopscopepodites']<-'Cyclopoids'
names(zoop2022_2)[names(zoop2022_2)=='Density.Daphniagaleatamendotae']<-'Daphniagaleatamendotae'
names(zoop2022_2)[names(zoop2022_2)=='Density.Daphniaparvula']<-'Daphniaparvula'
names(zoop2022_2)[names(zoop2022_2)=='Density.Daphniaretrocurva']<-'Daphniaretrocurva'
names(zoop2022_2)[names(zoop2022_2)=='Density.Daphniaspp']<-'Daphniaspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Diacyclopsthomasi']<-'Diacyclopsthomasi'
names(zoop2022_2)[names(zoop2022_2)=='Density.Daphnialongiremis']<-'Daphnialongiremis'
names(zoop2022_2)[names(zoop2022_2)=='Density.Chydorussp']<-'Chydorussp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Diacyclopsnanus']<-'Diacyclopsnanus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Diaphanosomaspp']<-'Diaphanosomaspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.DreissenidVeligers']<-'Dreissenidveligers'
names(zoop2022_2)[names(zoop2022_2)=='Density.Eubosminacoregoni']<-'Eubosminacoregoni'
names(zoop2022_2)[names(zoop2022_2)=='Density.Epischuralacustris']<-'Epischuralacustris'
names(zoop2022_2)[names(zoop2022_2)=='Density.Euchlanisspp']<-'Euchlanisspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Eurycercuslamellatus']<-'Eurycercuslamellatus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Eurytemoracarolleea']<-'Eurytemoracarolleea'
names(zoop2022_2)[names(zoop2022_2)=='Density.Gastropusstylifer']<-'Gastropus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Harpacticoidcopepodites']<-'Harpacticoida'
names(zoop2022_2)[names(zoop2022_2)=='Density.Holopediumgibberum']<-'Holopediumgibberum'
names(zoop2022_2)[names(zoop2022_2)=='Density.Kellicottialongispina']<-'Kellicottia'
names(zoop2022_2)[names(zoop2022_2)=='Density.Lepadellaspp']<-'Lepadellaspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Leptodiaptomusashlandi']<-'Leptodiaptomusashlandi'
names(zoop2022_2)[names(zoop2022_2)=='Density.Leptodiaptomusminutus']<-'Leptodiaptomusminutus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Leptodiaptomussicilis']<-'Leptodiaptomussicilis'
names(zoop2022_2)[names(zoop2022_2)=='Density.Leptodorakindti']<-'Leptodorakindti'
names(zoop2022_2)[names(zoop2022_2)=='Density.Limnocalanuscopepodites']<-'Limnocalanuscopepodites'
names(zoop2022_2)[names(zoop2022_2)=='Density.Limnocalanusmacrurus']<-'Limnocalanusmacrurus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Mesocyclopsedax']<-'Mesocyclopsedax'
names(zoop2022_2)[names(zoop2022_2)=='Density.Microcyclopsvaricans']<-'Microcyclopsvaricans'
names(zoop2022_2)[names(zoop2022_2)=='Density.Mysisdiluviana']<-'Mysisspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Polyarthraspp']<-'Polyarthra'
names(zoop2022_2)[names(zoop2022_2)=='Density.Polyphemuspediculus']<-'Polyphemuspediculus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Pompholyxsulcata']<-'Pompholyxsulcata'
names(zoop2022_2)[names(zoop2022_2)=='Density.Ptygura']<-'Ptygura'
names(zoop2022_2)[names(zoop2022_2)=='Density.Rotifera']<-'Rotifera'
names(zoop2022_2)[names(zoop2022_2)=='Density.Senecellacalanoides']<-'Senecellacalanoides'
names(zoop2022_2)[names(zoop2022_2)=='Density.Senecellacopepodites']<-'Senecellacopepodites'
names(zoop2022_2)[names(zoop2022_2)=='Density.Skistodiaptomusoregonensis']<-'Skistodiaptomusoregonensis'
names(zoop2022_2)[names(zoop2022_2)=='Density.Stephanoceros']<-'Stephanoceros'
names(zoop2022_2)[names(zoop2022_2)=='Density.Synchaetaspp']<-'Synchaetaspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Trichotriaspp']<-'Trichotriaspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Tropocyclopsprasinusmexicanus']<-'Tropocyclopsprasinusmexicanus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Tylotrochamonopus']<-'Tylotrochamonopus'
names(zoop2022_2)[names(zoop2022_2)=='Density.Cercopagispengoi']<-'Cercopagispengoi'
names(zoop2022_2)[names(zoop2022_2)=='Density.Ceriodaphniaspp']<-'Ceriodaphniaspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Alonella']<-'Alonella'
names(zoop2022_2)[names(zoop2022_2)=='Density.Bdelloidea']<-'Bdelloidea'
names(zoop2022_2)[names(zoop2022_2)=='Density.Cephalodellaspp']<-'Cephalodellaspp'
names(zoop2022_2)[names(zoop2022_2)=='Density.Epischuracopepodites']<-'Epischuracopepodites'


#remove columns that were merged to another column in 2022
zoop2022_3<-zoop2022_2
zoop2022_3<-zoop2022_3 %>% subset(-c(Density.Ascomorphaeucaudis,Density.Ascomorphaovalis))
zoop2022_3<-zoop2022_3 %>% subset(-c(Density.Bosminalongirostris,Density.Bosminaspp))
zoop2022_3<-zoop2022_3 %>% subset(-c(Density.Keratellacochlearis,Density.Keratellacrassa,Density.Keratellaearlinae,Density.Keratellahiemalis,Density.Keratellaquadrata))
zoop2022_3<-zoop2022_3 %>% subset(-c(Density.Lecaneflexilis,Density.Lecaneludwigi,Density.Lecanespp,
                                     Density.Monostylacopeis,Density.Monostylaspp,Density.Monostylalunaris
                                     ,Density.Notholcaacuminata,Density.Notholcafoliacea,Density.Notholcalaurentiae,Density.Notholcasquamula
                                     ,Density.Ploesomahudsoni,Density.Ploesomalenticulare,Density.Ploesomatruncatum,
                                     Density.Trichocercacylindrica,Density.Trichocercamulticrinis,Density.Trichocercaporcellus,Density.Trichocercarousseleti,Density.Trichocercaspp))

writexl::write_xlsx(zoop2022_3, "Zoop2022_315364.xlsx")
#reimport, moved 64 net to 153 and deleted parry sound 18
Zoop2022_4 <- read_excel("Zoop2022_315364.xlsx")


#merge columns in 2017
#Daphnia multiple species to -> Daphnia
zoop2017_2$Daphnia<-zoop2017_2$`Density.Daphnia galeata mendotae`+zoop2017_2$`Density.Daphnia retrocurva`+zoop2017_2$`Density.Daphnia retrocurva`
zoop2017_2<-zoop2017_2 %>% select(-c(`Density.Daphnia galeata mendotae`, `Density.Daphnia pulex`, `Density.Daphnia retrocurva`))

#zoop 2017 rename columns
zoop2017_2<-zoopdensity2017wide
names(zoop2017_2)[names(zoop2017_2)=='Density.Acanthocyclops']<-'Acanthocyclops'
names(zoop2017_2)[names(zoop2017_2)=='Density.Asplanchna']<-'Asplanchna'
names(zoop2017_2)[names(zoop2017_2)=='Density.Alona']<-'Alona'
names(zoop2017_2)[names(zoop2017_2)=='Density.Bosmina']<-'Bosmina'
names(zoop2017_2)[names(zoop2017_2)=='Density.Bythotrephes longimanus or spiny water flea']<-'Bythotrepheslongimanus'
names(zoop2017_2)[names(zoop2017_2)=='Density.Calanoida or calanoids']<-'Calanoida'
names(zoop2017_2)[names(zoop2017_2)=='Density.Conochilus']<-'Conochilus'
names(zoop2017_2)[names(zoop2017_2)=='Density.Cephalodellaspp']<-'Cephalodellaspp'
names(zoop2017_2)[names(zoop2017_2)=='Density.Copepoda or copepods']<-'Copepoda'
names(zoop2017_2)[names(zoop2017_2)=='Density.Cyclopoida or cyclopoids']<-'Cyclopoids'
names(zoop2017_2)[names(zoop2017_2)=='Density.Diacyclops thomasi']<-'Diacyclopsthomasi'
names(zoop2017_2)[names(zoop2017_2)=='Density.Gastropus']<-'Gastropus'
names(zoop2017_2)[names(zoop2017_2)=='Density.Kellicottia']<-'Kellicottia'
names(zoop2017_2)[names(zoop2017_2)=='Density.Keratella']<-'Keratella'
names(zoop2017_2)[names(zoop2017_2)=='Density.Leptodiaptomus ashlandi']<-'Leptodiaptomusashlandi'
names(zoop2017_2)[names(zoop2017_2)=='Density.Leptodiaptomus minutus']<-'Leptodiaptomusminutus'
names(zoop2017_2)[names(zoop2017_2)=='Density.Leptodiaptomus sicilis']<-'Leptodiaptomussicilis'
names(zoop2017_2)[names(zoop2017_2)=='Density.Limnocalanus macrurus']<-'Limnocalanusmacrurus'
names(zoop2017_2)[names(zoop2017_2)=='Density.Notholca']<-'Notholca'
names(zoop2017_2)[names(zoop2017_2)=='Density.Polyarthra']<-'Polyarthra'
names(zoop2017_2)[names(zoop2017_2)=='Density.Skistodiaptomus oregonensis']<-'Skistodiaptomusoregonensis'
names(zoop2017_2)[names(zoop2017_2)=='Density.Synchaeta']<-'Synchaeta'
names(zoop2017_2)[names(zoop2017_2)=='Density.Tropocyclops prasinus mexicanus']<-'Tropocyclopsprasinusmexicanus'
names(zoop2017_2)[names(zoop2017_2)=='Density.Dreissena']<-'Dreissenidveligers'
names(zoop2017_2)[names(zoop2017_2)=='Density.Leptodiaptomus']<-'Leptodiaptomusspp'
names(zoop2017_2)[names(zoop2017_2)=='Density.Mesocyclops edax']<-'Mesocyclopsedax'
names(zoop2017_2)[names(zoop2017_2)=='Density.Senecella  calanoides']<-'Senecellacalanoides'
names(zoop2017_2)[names(zoop2017_2)=='Density.Harpacticoida or harpacticoids']<-'Harpacticoida'
names(zoop2017_2)[names(zoop2017_2)=='Density.Mysis spp.']<-'Mysisspp'
names(zoop2017_2)[names(zoop2017_2)=='Density.Collotheca']<-'Collotheca'
names(zoop2017_2)[names(zoop2017_2)=='Density.Chydorus sphaericus']<-'Chydorussp'
names(zoop2017_2)[names(zoop2017_2)=='Density.Daphnia galeata mendotae']<-'Daphniagaleatamendotae'
names(zoop2017_2)[names(zoop2017_2)=='Density.Lecane']<-'Lecane'
names(zoop2017_2)[names(zoop2017_2)=='Density.Cephalodella']<-'Cephalodella'
names(zoop2017_2)[names(zoop2017_2)=='Density.Lepadella']<-'Lepadella'
names(zoop2017_2)[names(zoop2017_2)=='Density.Macrotrachela']<-'Macrotrachela'
names(zoop2017_2)[names(zoop2017_2)=='Density.Monostyla']<-'Monostyla'
names(zoop2017_2)[names(zoop2017_2)=='Density.Holopedium gibberum']<-'Holopediumgibberum'
names(zoop2017_2)[names(zoop2017_2)=='Density.Ascomorpha']<-'Ascomorpha'
names(zoop2017_2)[names(zoop2017_2)=='Density.Eubosmina coregoni (formerly Bosmina coregoni)']<-'Eubosminacoregoni'
names(zoop2017_2)[names(zoop2017_2)=='Density.Wierzejskiella']<-'Wierzejskiella'
names(zoop2017_2)[names(zoop2017_2)=='Density.Daphnia retrocurva']<-'Daphniaretrocurva'
names(zoop2017_2)[names(zoop2017_2)=='Density.Epischura lacustris']<-'Epischuralacustris'
names(zoop2017_2)[names(zoop2017_2)=='Density.Leptodora kindtii']<-'Leptodorakindti'
names(zoop2017_2)[names(zoop2017_2)=='Density.Ploesoma']<-'Ploesoma'
names(zoop2017_2)[names(zoop2017_2)=='Density.Polyarthra']<-'Polyarthra'
names(zoop2017_2)[names(zoop2017_2)=='Density.Pompholyx']<-'Pompholyx'
names(zoop2017_2)[names(zoop2017_2)=='Density.Trichocerca']<-'Trichocerca'
names(zoop2017_2)[names(zoop2017_2)=='Density.Daphnia pulex']<-'Daphniapulex'
names(zoop2017_2)[names(zoop2017_2)=='Density.Proales']<-'Proales'
names(zoop2017_2)[names(zoop2017_2)=='Density.Rotifera or rotifers']<-'Rotifera'
names(zoop2017_2)[names(zoop2017_2)=='Density.Euchlanis']<-'Euchlanis'
names(zoop2017_2)[names(zoop2017_2)=='Density.Stephanoceros']<-'Stephanoceros'
names(zoop2017_2)[names(zoop2017_2)=='Density.Encentrum']<-'Encentrum'
names(zoop2017_2)[names(zoop2017_2)=='Density.Diaphanosoma']<-'Diaphanosomaspp'

writexl::write_xlsx(zoop2017_2, "Zoop2017names.xlsx")
writexl::write_xlsx(zoop2022_3, "Zoop2022names.xlsx")


#Daphnia multiple species to -> Daphnia
Zoop2022_4$Daphnia<-Zoop2022_4$Daphniaspp+Zoop2022_4$Daphniaparvula+Zoop2022_4$Daphniaretrocurva+Zoop2022_4$Daphniagaleatamendotae+Zoop2022_4$Daphnialongiremis
Zoop2022_4$Euchlanis<-
  Zoop2022_4<-Zoop2022_4 %>% select(-c(Daphniaspp,Daphniaparvula,Daphniaretrocurva,Daphniagaleatamendotae,Daphnialongiremis))
names(Zoop2022_4)[names(Zoop2022_4)=='Euchlanisspp']<-'Euchlanis'
names(Zoop2022_4)[names(Zoop2022_4)=='Synchaetaspp']<-'Synchaeta'
names(Zoop2022_4)[names(Zoop2022_4)=='Cephalodellaspp']<-'Cephalodella'
names(Zoop2022_4)[names(Zoop2022_4)=='Lepadellaspp']<-'Lepadella'

#combine 2022 and 2017 datasets
zoop22comb<-zoop2022_3
zoop17comb<-zoop2017_2
zoop22comb2<-Zoop2022_4
zoop22comb2$UnifiedDate<-as.Date(zoop22comb2$UnifiedDate, tryFormats = "%Y-%m-%d")
zoop17comb$UnifiedDate<-as.Date(zoop17comb$UnifiedDate, tryFormats = "%Y-%m-%d")
zoopcomb2<-rbindlist(list(zoop22comb2, zoop17comb), fill =TRUE)

writexl::write_xlsx(zoopcomb2, "Zoopcombforrich.xlsx")
