library(tidyverse);
library(doParallel);
library(foreach)
source("ModelwCovs.R")

#Do this in parallel (computation is intense!)
cl <- makeCluster(18)
registerDoParallel(cl)
getDoParWorkers()
getDoParRegistered()

#make a neighborhood
k1=30;k2=30;kk=k1*k2
NBHD4=rectanggrid.4nbrs(k1,k2)
nbss=rectanggrid.4nbrs(k1,k2)
locs=1:kk

#initial values
sig2 = 100
bmeans <- c(0.5); gmeans <- c(1); spatpars <- c(0.5,0.9,3)
truths=c(spatpars,bmeans,gmeans,sig2)
iout=1
NSIM_out =1000
NSIM_in = 500
truths=c(spatpars,bmeans,gmeans,sig2)
ys <- rbinom(prob = exp(bmeans[1])/(1+exp(bmeans[1])), size = 1, n = k1*k2) 
zs <- rnorm(k1*k2, mean = gmeans, sd = sqrt(sig2))

#Use this function in parallel
#The actual bootstrap function is at the bottom of this script
parflex <- function(){
  ys <- rbinom(prob = exp(bmeans[1])/(1+exp(bmeans[1])), size = 1, n = k1*k2) 
  zs <- rnorm(k1*k2, mean = gmeans, sd = sqrt(sig2))
  boots <- bootstrapbivar(spatpars, bmeans, gmeans, sig2=sig2, ys, zs, k1,k2, nbss, B=500, M=NSIM_in, S=20)
  c(boots$MCest, c( boots$CI[,1] <= truths & truths <= boots$CI[,2]))
}

#Parameter Set 1
sig2=1
spatpars <- c(.2,.5,1)
truths=c(spatpars,bmeans,gmeans,sig2)
iout=1
#Run it in parallel (24 hours+)
CIweak <- foreach(icount(NSIM_out), .combine=rbind) %dopar% {parflex()}
apply(CIweak, 2, mean) #MC means and CI coverage
apply(CIweak, 2, sd) #MC SEs

#Parameter Set 2
sig2=100
spatpars <- c(5,.9,3)
iout=1
truths=c(spatpars,bmeans,gmeans,sig2)
CIstrong <- foreach(icount(NSIM_out), .combine=rbind) %dopar% {parflex()}
apply(CIstrong, 2, mean)
apply(CIstrong, 2, sd)

####################################################
####Making Plots and Summarizing Data###############
library(ggplot2); library(RColorBrewer)
cis <- read.csv(file="try2str.1-1000.csv",header=T)[,-1]
head(cis)
round(apply(cis,2,FUN=sd),3)[1:6]
round(apply(cis,2,mean),3)
g <- ggplot(data=cis)
g+ geom_histogram(aes(x=mu), fill="black", color="white") +
  labs(x=expression(mu)) + ggtitle(expression("Histogram of "*mu*" Estimates"))+
  theme(plot.title = element_text(size = 24, face = "bold", hjust=0.5),
        axis.title.y = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text = element_text(size=16)) + scale_x_continuous(breaks = round(seq(-40, 15, by = 5)))

g+geom_point(aes(x=delta,y=mu, color=as.factor(CIdelta+(CImu)*(CIdelta+5)))) +
  labs(x=expression(delta), y=expression(mu)) +
  ggtitle(expression("Estimated Global Mean Parameters"))+
  theme(plot.title = element_text(size = 24, face = "bold", hjust=0.5),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12), legend.key.size = unit(1, "cm"),
        axis.title.y = element_text(size=24, angle=0, vjust=0.5),
        axis.title.x = element_text(size=24),
        axis.text = element_text(size=16)) +
  scale_x_continuous(breaks = round(seq(-4, 2, by = 1)))+
  scale_color_manual(name=expression("95% CIs contain true "*mu*" and "*delta*"?"),
                     labels = c("Neither", expression(delta*" only"),
                    expression(mu*" only"),"Both"),values = brewer.pal(name="RdGy", n=4))


summary(lm(data=cis, formula=mu~delta))
table(cis$CIdelta,cis$CImu)

cisw <- read.csv(file="try2.1-1000.csv", header=T)[,-1]
names(cisw)[4:5] <- c("mu","delta")
g <- ggplot(data=cisw)
g+geom_histogram(aes(x=mu))
head(cisw)
summary(lm(data=cisw, formula=mu~delta))

####################################################################
####################################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
##Bootstrap Simulation Function to estimate CIs#############

bootstrapbivar <- function(spatials,bmeanpars, gmeanpars, sig2, startys, startzs, k1, k2, nbrs, B, M, S){
  ##Compute bootstrap standard errors 
  ##spatials/bmeanpars/gmeanpars - vector of estimated parameters
  ##startys - arbitrary vector of 0's and 1's
  ##startzs - arbitrary vector of gaussian rv
  ##B - number of iterations until burn-in
  ##M - number of bootstrap samples
  ##S - thinning value
  locs=1:(k1*k2);n=length(locs)
  ests <- matrix(0, nrow=M,ncol=6)
  dataset <- spatial.genfieldbivar.covs(startys,startzs,nbss,M=B,
                                        spatials,bmeanpars, gmeanpars, sigmasq=sig2, bcovs=NULL, gcovs=NULL)
  
  datasetys <- dataset[,1]
  datasetzs <- dataset[,2]
  ps <- initialparams.cov(datasetys,datasetzs,nbss,bcovs=NULL,gcovs=NULL)
  initialests <- optim(par=ps,f=logplikbivar.covs, k1=k1, k2=k2,ys=datasetys,
                      zs=datasetzs, bcovs=NULL, gcovs=NULL, nbrs=nbss,hessian=T,control=list(maxit=500))$par
  cntzeros <- 0
  cnt <- 1
  repeat{
    newobs <- spatial.genfieldbivar.covs(datasetys,datasetzs,nbss,M=S,
                                         initialests[1:3],initialests[4], initialests[5], sigmasq=initialests[6],
                                  bcovs=NULL, gcovs=NULL)
    newys <- newobs[,1]
    newzs <- newobs[,2]
    estsfull <- try(optim(par=initialests,f=logplikbivar.covs, k1=k1, k2=k2,ys=newys,
                      zs=newzs, bcovs=NULL, gcovs=NULL, nbrs=nbss,hessian=T,control=list(maxit=500)))
    if(length(estsfull)>1){
    ests[cnt,] <- estsfull$par;
    datasetys <- newys
    datasetzs <- newzs
    if(cnt%%5==0){print(paste(cnt,".", iout))}
    cnt <- cnt+1
    }
    else{cntzeros <- cntzeros+1}
    
    if(cnt==(M+1)){break}
    if(cntzeros == 20){break}
  }
  CIupper <- 2*initialests - apply(ests, MARGIN = 2, FUN = quantile, prob= .025)
  CIlower <- 2*initialests - apply(ests, MARGIN = 2, FUN = quantile, prob= .975)
  return(list(CI= cbind(CIlower,CIupper), MCest = initialests))
}

