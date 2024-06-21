#R code for Lake Huron 2022 CSMI survey report EPA ORD
#Courtney Larson and Noah Grode

install.packages("reshape2")

#libraries
library(reshape2)
#programs

#color vectors

#Upload files
Zoopreshap <- read.csv("ZoopforR6.21.24.csv", header=T)

#melt to long format
Zoopcount153<-melt(Zoopreshap, value.name="Count", 
               variable.name = "Species")
