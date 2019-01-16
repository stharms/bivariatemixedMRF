library(sp);library(tidyverse);library(ggplot2);library(raster);
library(ggmap)
library(maps)
library(mapdata)
library(rgeos)

#read in data
chickens <- read.csv(file="chickens.csv")
geoproj <- " +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
chickens2 <- chickens #%>% filter(year >= 1990 , month <=5, month >=3)

#keep the variables we need and convert to spatial points data frame
chickens2 <- data.frame( long=chickens2$decimalLongitude,lat=chickens2$decimalLatitude, month=chickens2$month, year=chickens2$year, presence = 1)
spchicks2 <- SpatialPointsDataFrame(coords=cbind(chickens2$long,chickens2$lat),data=chickens2[,3:5],proj4string = CRS(geoproj))

states <- map_data("state")
greatplains <- c("south dakota", "north dakota", "kansas","nebraska","oklahoma","iowa",
                 "colorado","wyoming","minnesota")
greatplains <- c("nebraska")
midwest <- subset(states, region %in% greatplains)#subset(states, region %in% ks)
counties <- map_data("county")
mw_county <- subset(counties, region %in% greatplains)#subset(counties, region %in% ks)

ggplot(data = midwest, aes(x = long, y = lat)) + 
  geom_polygon(color = "black", aes(group=group))+ coord_fixed(1.3) +
  geom_polygon(data = mw_county, aes(group=group), fill = NA, color = "white")+
  geom_point(data = chickens2, mapping = aes(x = long, y = lat), color = "red") + xlim(-105,-95) + ylim(38,45)

