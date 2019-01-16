library(tidyverse);library(viridisLite); library(scoringRules)
source("ReducedBinaryLLPs.R"); source("ReducedGaussLLPs.R");
source("ReducedBivariateLLPs.R")
source("ScoringRules.R")
source("Sstatistic.R");
source("BinaryUpperBounds.R")
source("ModelwCovs.R")

####################################################################
####################################################################

##################################################
##################################################
bootstrapbivar <- function(spatials,bmeanpars, gmeanpars, sig2, startys, startzs, bcovt, gcovt, k1, k2, nbrs, B, M){
  ##Compute bootstrap standard errors 
  ##spatials/bmeanpars/gmeanpars - vector of estimated parameters
  ##startys - arbitrary vector of 0's and 1's
  ##startzs - arbitrary vector of gaussian rv
  ##B - number of iterations until burn-in
  ##M - number of bootstrap samples
  locs=1:(k1*k2);n=length(locs)
  estsfull <-matrix(nrow=M, ncol=length(c(spatials,bmeanpars,gmeanpars))+1)
  #colnames(estsfull)<- c("rho","eta.z","eta.y","b0.y","b1.y","b0.z","b1.z","sig.sq")
  ests.zcov<-ests.ycov<-ests.nocov<-estsfull
  bscore <- array(dim=c(M,1,3)); gscore<-bscore
  S=10
  
  cntzeros <- NULL
  cnt <- 1
  newys<-startys; newzs<-startzs
  ps <- initialparams.cov(newys,newzs,nbss,bcovt,gcovt)
  ps[3] <- Dbinary(newys,nbrs,k1=k1,k2=k2)[1]
  repeat{
    dataset <- spatial.genfieldbivar.covs(newys,newzs,nbss,M=S,
                                          spatials,bmeanpars, gmeanpars, sigmasq=sig2, bcovs=bcovt, gcovs=gcovt)
    newys <- dataset[,1]
    newzs <- dataset[,2]
    ests <- try(optim(par=ps,f=logplikbivar.covs, k1=k1, k2=k2,ys=newys,
                  zs=newzs, bcovs=bcovt,gcovs=gcovt, nbrs=nbss,hessian=T,control=list(maxit=1000)),silent=TRUE)
    if(length(ests)>1){
    estsfull[cnt,] <- c(ests$par)
    bscore[cnt,1,]<- brierscore(estsfull[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=estsfull[cnt,4:7], gcovt, gcovpars=estsfull[cnt,8], estsfull[cnt,9])
    gscore[cnt,1,]<- GaussScore(estsfull[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=estsfull[cnt,4:7], gcovt, gcovpars=estsfull[cnt,8], estsfull[cnt,9])
    cnt<-cnt+1
    }
    print(paste(cnt-1, "full model"))
    if(cnt==(M+1)){break}
  }
  
  ests<-cbind(estsfull)
  #hess[,cnt] <- vechalf(ests$hessian)
  #conv[cnt] <- ests$convergence
  colnames(gscore)<- colnames(bscore)<-c("Full")
  return(list(ests = ests, CRPSg = gscore, CRPSb = bscore))
}
