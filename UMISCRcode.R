<<<<<<< HEAD
#Install Packages
install.packages("ggplot2")
install.packages("ggpubr")
library(ggplot2)
library(ggpubr)
library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(lubridate)

#Upload datasets
ZoopforR6_18_24 <- read_csv("~/CSMI/ZoopforR6.18.24.csv")
View(ZoopforR6_18_24)

Zoop64countforR <- read_csv("~/CSMI/Zoop64countforR.csv")
View(Zoop64countforR)

CalibrationMeta <- read_excel("~/CSMI/CalibrationMeta.xlsx")
CM <- (CalibrationMeta)
View(CalibrationMeta)

Zoopcount153long<-melt(ZoopforR6_18_24, value.name="Count", 
                       variable.name = "Species")
#Merging our metadata with longform species data
Zoopcount153Mer<-merge(Zoopcount153long,CM, by="SiteID")
#creating column density by taking count/volume columns
Zoopcount153Mer$Density<-Zoopcount153Mer$Count/Zoopcount153Mer$Volume

Zoop64count<-melt(Zoop64countforR, value.name = "Count", variable.name = "Species")
Zoop64countMer<-merge(Zoop64count, CM, by="SiteID")
Zoop64countMer$Density<-Zoop64countMer$Count/Zoop64countMer$Volume

#Cor
cor.test(Zoopcount153Mer$Species, Zoopcount153Mer$Density)



=======
#Install Packages
install.packages("ggplot2")
install.packages("ggpubr")
library(ggplot2)
library(ggpubr)
library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(lubridate)

#Upload datasets
ZoopforR6_18_24 <- read_csv("~/CSMI/ZoopforR6.18.24.csv")
View(ZoopforR6_18_24)

Zoop64countforR <- read_csv("~/CSMI/Zoop64countforR.csv")
View(Zoop64countforR)

CalibrationMeta <- read_excel("~/CSMI/CalibrationMeta.xlsx")
CM <- (CalibrationMeta)
View(CalibrationMeta)

Zoopcount153long<-melt(ZoopforR6_18_24, value.name="Count", 
                       variable.name = "Species")
#Merging our metadata with longform species data
Zoopcount153Mer<-merge(Zoopcount153long,CM, by="SiteID")
#creating column density by taking count/volume columns
Zoopcount153Mer$Density<-Zoopcount153Mer$Count/Zoopcount153Mer$Volume

Zoop64count<-melt(Zoop64countforR, value.name = "Count", variable.name = "Species")
Zoop64countMer<-merge(Zoop64count, CM, by="SiteID")
Zoop64countMer$Density<-Zoop64countMer$Count/Zoop64countMer$Volume

#Cor
cor.test(Zoopcount153Mer$Species, Zoopcount153Mer$Density)



>>>>>>> c4ef7cfcec9d426c72279303a92826ff5cb98a6a
