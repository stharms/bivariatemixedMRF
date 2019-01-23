library(tidyverse); library(raster);library(sp);
library(scoringRules)
library("viridisLite"); library(RColorBrewer); library(ggplot2)
source("ReducedBivariateLLPs.R")
source("ScoringRules.R")
source("Sstatistic.R");
source("BinaryUpperBounds.R")
source("ModelwCovs.R")
source("bootfcts.R")

##################################
#here we can test out the model on our new data
k1=31;k2=31; nbss <-rectanggrid.4nbrs(k1,k2) #neighborhood and lattice size
newys <- nesubpts$chicken; newzs <- (nesubpts$NDVI)#these are the variables
ycov <- as.matrix.data.frame(nesubpts[,c(7,9,11)]);
zcov <- NULL
#get some initial values for likelihood estimation
ips <- initialparams.cov(newys , newzs , nbss , bcovs = ycov , gcovs=zcov)
#S-value for binary field
Dbinary(newys,nbss,k1=31,k2=31)[1]
#estimate the model parameters for the full model using pseudolikelihood
pps <-  optim(par=ips,f=logplikbivar.covs, k1=k1, k2=k2,ys=newys, zs=newzs,
              bcovs=ycov, gcovs=zcov,
              nbrs=nbss,hessian=T,control=list(maxit=1000))
round(pps$par,5)

#CRPS scores
brierscore(pps$par[1:3],newys, newzs, nbss,
           ycov, bcovpars=pps$par[4:7], zcov, gcovpars=pps$par[8], pps$par[9])
GaussScore(pps$par[1:3],newys, newzs, nbss,
           ycov, bcovpars=pps$par[4:7], zcov, gcovpars=pps$par[8], pps$par[9])

chickenboots.t2 <- bootstrapbivar(pps$par[1:3], pps$par[4:7], pps$par[8], pps$par[9], newys, newzs, ycov,zcov, k1=31,k2=31, nbrs=nbss,
                               B=500, M=500, S=10)


#bootstrap estimates
round(apply(chickenboots.t2$ests,2,mean),3)
round(apply(chickenboots.t2$ests,2,median),3)
round(apply(chickenboots.t2$ests,2,sd),3)
#CIs
CIs <- cbind(2*pps$par - apply(chickenboots.t2$ests,2,quantile,probs=.95),
             2*pps$par - apply(chickenboots.t2$ests,2,quantile,probs=.05))
round(CIs, 3)

#Plots of the original data
plot(nesub$chicken, axes=F, col = c("grey","dark red"), main = "Greater Prairie Chickens", useRaster=F)
plot(nesub$NDVI, col=rev(brewer.pal(9,name="Greens")), axes=F, main = "NDVI", useRaster= F, bty="n")
plot(nesub$WindGroup, col=rev(gray.colors(30)), axes = F)
plot(nesub$pctsand, col=rev(gray.colors(30)), axes = F)
plot(nesub$watersupply, col = rev(gray.colors(30)), axes = F, useRaster=F, bty="n")

#histograms
g <- ggplot(data=nesubpts)
g+ geom_histogram(aes(x=NDVI), color="black", fill="dark green")+labs(title="NDVI",x="NDVI", y = "Count")+
  theme(plot.title = element_text(size = 22, face = "bold",hjust = 0.5),
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