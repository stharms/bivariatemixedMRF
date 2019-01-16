rm(list = ls(all = T))
library("tidyverse")
library("raster")
library("sp")
library("foreign")
library("viridisLite") #for better color palettes
library(ggplot2)
library(FastImputation)
source("ModelwCovs.R")
source("ReducedBivariateLLPs.R")
source("BinaryUpperBounds.R")

#read in the data
wheat <- read.csv("WheatTrial.csv", header=T)
wheatsp <- SpatialPointsDataFrame(coords = wheat[,c(6,7)], data = wheat[,c(8:11)])
wheatdf <- wheat[,c(6,7,8:11)] ; names(wheatdf)<- c("row","col", "DayHead","GrossYld","KerSpk","TotKerWt")

#imputing the 5 missing values
ttrain <- TrainFastImputation(wheatdf, idvars=c("row","col"))
impwheat <- FastImputation(wheatdf, ttrain, verbose = TRUE)
wheatdf <- impwheat

#re-coding the data
wheatdf$DayHead[which(wheat$DH_SR_FI>=79)] <- 1; wheatdf$DayHead[which(wheat$DH_SR_FI<79)] <-0;
wheatdf$DayHead[which(is.na(wheat$DH_SR_FI))] <- 0
wheatr <- rasterFromXYZ(wheatdf)

#using a regular rectangular grid
wheatsub <- crop(x=wheatr,
                 y = extent(wheatr, 1, 20, 1, 40))

########################################################################
#plots of the data
plot(wheatsub$DayHead, col=gray.colors(2), axes=F, main = "Days to Head")
plot(wheatr$GrossYld/1000, col=viridis(30), axes=F, main = "Gross Yield")
plot(wheatr$KerSpk)
plot(wheatr$TotKerWt, col=magma(30), axes=F, main = "Thousand-Kernel Weight")

ggplot(data = wheatdf, mapping = aes(col,row) ) + geom_raster(aes(fill = as.factor(DayHead)), show.legend=T)+
  scale_fill_grey(name="Days >= 79")+ ggtitle("Days to Head")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), 
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=, face = "bold"), legend.key.size = unit(1, "cm"))
pdf("wheatDayHead.pdf")

ggplot(data = wheatdf, mapping = aes(col,row) ) + geom_raster(aes(fill = GrossYld))+
  scale_fill_gradientn(name= "Bushels",colours = magma(8)) + ggtitle("Gross Yield")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=, face = "bold"), legend.key.size = unit(1, "cm"))
pdf("wheatGrossYield.pdf")

ggplot(data = wheatdf, mapping = aes(col,row) ) + geom_raster(aes(fill = TotKerWt))+
  scale_fill_gradientn(name= "Weight (g)",colours = viridis(8)) + ggtitle("Thousand-Kernel Weight")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=, face = "bold"), legend.key.size = unit(1, "cm"))
pdf("wheatTKW.pdf")

wheatsubpts <- data.frame(rasterToPoints(wheatsub,spatial=F))
################################################################################
##############

