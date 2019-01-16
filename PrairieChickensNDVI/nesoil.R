################################################
library(tidyverse); library(raster);library(sp);
library(scoringRules); library(RColorBrewer)
library("foreign")
library("viridisLite")



#overlay code
#takes a while
get_CDLoverMU_tables <- function(areasym, mulayer, cdlfile){
  ## NEB Mapunit shape
  mu.ne <- rgdal::readOGR(dsn=".", layer="soilmu_a_ne111")
  mukeyname <- levels(mu.ne@data$MUKEY);
  ## table CDL pixels by MUKEY
  result <- matrix(nrow = length(mukeyname), ncol = length(count.col))
  
  foldername <- paste0("Figures/CDLOverMU/", areasym)
  dir.create(foldername)
  
  #convert that set to points
  cdlpoint <- raster::rasterToPoints(cdlfile, spatial = T,
                                     fun = function(x) x %in% count.col)

  #project those points onto the new coordinate projection
  cdlpoint.proj <- sp::spTransform(cdlpoint, CRS(projargs = geoproj))
  mupoint <- cdlpoint.proj
  iter=0
  repeat{
    iter <- iter + 1
    #get polygon for mapunit
    musub <- subset(mu.ne, MUKEY == mukeyname[iter])
    #overlay the points onto the mapunit polygon
    overlay <- sp::over(cdlpoint.proj, geometry(musub))
    mupoint[!is.na(overlay),1 ] <- as.numeric(mukeyname[iter])
    
    if(iter == length(mukeyname)) {break}
  }
  #stop here for now
  return(mupoint)
  ## save result ####
  colnames(result) <- paste("Category", count.col, sep = "_")
  result <- result %>% as_tibble() %>% 
    mutate(MUKEY = mukeyname) %>%
    dplyr::select(MUKEY, everything())
  filename <- paste0("_MU_CDL_Table.csv")
  write.csv(result, file = filename, row.names = F)
  #return(result)
}

###########################################
geoproj <- " +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

#this is the CDL data and associated dbf (legend)
cdl.corn <- raster("nebCDL.tif",crs = geoproj)
plot(cdl.corn)
cdl.neb.dbf <- foreign::read.dbf("nebCDL.dbf")

#NDVI
ndvi.neb <-raster("ndvi4172018.tif", crs=geoproj) %>% projectRaster(to=cdl.corn, crs=crs(cdl.corn))
vci.neb <- raster("vci4032018.tif", crs=geoproj) %>% projectRaster(to=cdl.corn, crs=crs(cdl.corn))

dim(ndvi.neb)
dim(cdl.corn)
dim(vci.neb)

#ndvi.neb <- disaggregate(ndvi.neb, c(1452/174,1084/177))
#crs(ndvi.neb)<- crs(cdl.corn)

#tabulating the CDL data
count.col <- cdl.neb.dbf %>% filter( !is.na(CLASS_NAME)) %>%
  dplyr::select(VALUE) %>% unlist %>% unname
count.col <- unique(vci.neb$vci4032018)
format(object.size(cdl.corn), units = "Kb")

#here's the overlaying. The mapunit info (downloaded from Web Soil Survey) is in "soilmu_a_ia095)
neb<- get_CDLoverMU_tables("NE111", mulayer="soilmu_a_ne111", cdl.corn)
nebndvi <- get_CDLoverMU_tables("NE111", mulayer="soilmu_a_ne111", ndvi.neb)
nebvci <- get_CDLoverMU_tables("NE111", mulayer="soilmu_a_ne111", vci.neb)
#convert it to a raster layer
cornpt <- raster::rasterToPoints(cdl.corn, spatial = T)
ndvipts <- rasterToPoints(ndvi.neb, spatial = T)
vcipts <- rasterToPoints(vci.neb, spatial = T)


#name and reogranize the layers
cbnd <- data.frame(cdl=cornpt,mapu=neb$nebCDL)[,-4]; names(cbnd)<- c("cdl","x","y","mapunit")
cbndvi <- data.frame(ndvi=ndvipts, mapu=nebndvi$ndvi4172018)[,-4];names(cbndvi)<- c("NDVI","x","y","mapunit")
cbndvc <- data.frame(vci=vcipts, mapu=nebvci$vci4032018)[,-4];names(cbndvc)<- c("VCI","x","y","mapunit")
#there are a few NAs, which can be removed(?) for now
cbnd <- cbnd[which(cbnd$mapunit>1000),c(2,3,1,4)]
#Change all of the pixels to 1(=corn) or 0 (=not corn) and we have a binary variable
cbnd$corn[which(cbnd$corn!=1)]<-0
cbndvi <- cbndvi[which(cbndvi$mapunit>1000),c(2,3,1,4)]
cbndvc <- cbndvc[which(cbndvc$mapunit>1000),c(2,3,1,4)]
combined <- merge(cbndvi, cbnd, by=c("x","y"))[,c(1:2,4,5,3)]; names(combined)[3]<-"mapunit"
combined <- merge(combined, cbndvc, by=c("x","y"))[,c(1:6)];names(combined)[3]<-"mapunit";

rm(cbndvi); rm(cbndvc)

#This is the soil data
csr <- read.csv("NE111soil.csv", header = T)[,-2]; head(csr)
names(csr)<- c("musym","mapunit", "WindGroup","pH","pctsand","pctsilt","watersupply")
#Merge the soil data with our CDL data based on map unit
mt <- merge(combined,csr[,-1], by="mapunit")
#reorganize the columns
mt <- mt[,c(2:3,1,4:11)]
head(mt)

#convert it back to a raster
rprast <-rasterFromXYZ(mt, crs=crs(cdl.corn))
#plot the Corn Suitability Rating
plot(rprast$NDVI, col=inferno(100))
plot(rprast$VCI, col = inferno(100), axes=F)
#aggregate the data to a lower resolution
aggd <- aggregate(rprast, fun=mean, na.rm=T, fact=80)
dim(aggd)

#some plots of our data
#plot(aggc$corn)
plot(aggd$NDVI, col = inferno(100))
plot(aggd$VCI, col = inferno(100), axes=F)


#some plots to look at
plot(aggd$WindGroup, col=inferno(80));
plot(aggd$pH, col=inferno(80))
plot(aggd$pctsand, col=inferno(80))

#overlaying the prairie chickens
aggdpts <- rasterToPoints(aggd, spatial=T)
aggd2 <- rasterToPoints(aggd, spatial=F)
geoproj = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
chickensp <- SpatialPoints(coords=chickens2[,1:2], proj4string = CRS("+init=epsg:4326"))
chickenstrans <- spTransform(chickensp, crs(cdl.corn))
chickenslink <- crop(chickenstrans, extent(aggdpts))
dists <- apply(gDistance(aggdpts, chickenslink, byid=TRUE), 1, which.min)

aggdpts$chicken <- 0
i=0
repeat{
  i <- i+1
  if(i %in% dists){aggdpts$chicken[i]=1}
  if(i>=length(aggdpts$mapunit)){break}
}

aggd2<- cbind(aggd2, chicken=aggdpts$chicken)
summary(aggd2)
newrast <- rasterFromXYZ(aggd2)
plot(newrast$chicken, col=rev(gray.colors(2)))
#cropping the data just to look at a subset for now
nesub <- crop(x = newrast,
              y = extent(newrast,3,33,6, 36));dim(nesub);

#as a data frame for model estimation
nesubpts <- data.frame(rasterToPoints(nesub, spatial=F))
head(nesubpts)




##############################################################
#code for the plots
plot(nesub$chicken.2, axes=F, col = c("grey","dark red"))
plot(nesub$NDVI, col=rev(brewer.pal(9,name="Greens")), axes=F)
plot(nesub$WindGroup, col=rev(gray.colors(30)), axes = F)
plot(nesub$pctsand, col=rev(gray.colors(30)), axes = F)
plot(nesub$watersupply, col = rev(gray.colors(30)), axes = F, useRaster=F, bty="n")

g <- ggplot(data=nesubpts)
g+ geom_histogram(aes(x=VCI), color="black", fill="dark green")+labs(title="NDVI",x="NDVI", y = "Count")+
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        axis.text.x = element_text(size = 20), axis.title.y= element_text(size=14))
g+ geom_histogram(aes(x=WindGroup), color="black", fill="black")+labs(title="Wind Erodibility",x="", y = "Count")+
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        axis.text.x = element_text(size = 15), axis.title.y= element_text(size=14))
g+ geom_histogram(aes(x=pctsand), color="black", fill="black")+labs(title="Percent Sand Composition",x="", y = "Count")+
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        axis.text.x = element_text(size = 15), axis.title.y= element_text(size=14))
g+ geom_histogram(aes(x=watersupply), color="black", fill="black")+labs(title="Water Supply Rating",x="", y = "Count")+
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        axis.text.x = element_text(size = 15), axis.title.y= element_text(size=14))

nerast <- rasterFromXYZ(nesubpts)

