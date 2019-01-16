library(tidyverse);library(viridisLite); library(scoringRules)
source("ReducedBivariateLLPs.R")
source("ScoringRules.R")
source("Sstatistic.R");
source("BinaryUpperBounds.R")
source("ModelwCovs.R")

k1=20;k2=40; nbss <-rectanggrid.4nbrs(k1,k2) #neighborhood and lattice size
newys <- wheatsubpts$DayHead; newzs <- wheatsubpts$GrossYld/1000#these are the variables
#covariates
#ycov <- as.matrix.data.frame(wheatsubpts[,c(1,2,5,6)]); zcov <- as.matrix.data.frame(wheatsubpts[,c(1,2,5,6)])
ycov <- zcov<- wheatsubpts[,6]

#get some initial values for likelihood estimation
ips <- initialparams.cov(newys , newzs , nbss , bcovs = ycov , gcovs=zcov)
#estimate the full model
pps <- optim(par=ips,f=logplikbivar.covs, k1=k1, k2=k2,ys=newys, zs=newzs,
             bcovs=ycov, gcovs=zcov,
             nbrs=nbss,hessian=T,control=list(maxit=1000))
#parameter estimates
round(pps$par,5)

#S-values (moment estimators) of spatial dependence
Dgaussian(wheatsubpts$GrossYld/1000,nbss, nbins=25, sigma2=pps$par[length(pps$par)])
Dbinary(wheatsubpts$DayHead, nbss, k1, k2)
#########################################################################
#bootstrap confidence interval estimates (huge function is at the bottom of this script)
wheatboots <- bootstrapbivar(pps$par[1:3], pps$par[4:5], pps$par[6:7], pps$par[8], newys, newzs, ycov,zcov, k1=20,k2=40, nbrs=nbss,
                                              B=200, M=500)
#########################################################################
####################
#Estimate model parameters for all models considered
ipsf <- initialparams.cov(newys , newzs , nbss , bcovs = ycov , gcovs=zcov)
########################
#full
full <- optim(par=ipsf,f=logplikbivar.covs, k1=k1, k2=k2,ys=newys, zs=newzs,
             bcovs=ycov, gcovs=zcov,
             nbrs=nbss,hessian=T,control=list(maxit=1000))$par
######################
#no covariates
estsnc <- optim(par=initialparams.cov(newys,newzs,nbss),f=logplikbivar.covs, k1=k1, k2=k2,ys=newys,
                 zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
estsnc <-c(estsnc[1:4],0,estsnc[5],0,estsnc[6])
#####################
#bivariate non-spatial
bns <- optim(par=c(ipsf[1], ipsf[4:8]),f=logplikbivar.ns, k1=k1, k2=k2, gcovs=zcov, bcovs=ycov,
      zs=newzs, ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par
bns <-c(bns[1],0,0,bns[2:6])
#######################
#define another covariate for the univariate models
bcovg = cbind(ycov, newzs-mean(newzs))
gcovb = cbind(zcov, newys-mean(newys))
ps <- initialparams.cov(newys,newzs,nbss,bcovg,gcovb)
###############################################################
#univariate spatial
bestsunis <- optim(par=c(ps[3:6]),f=logplikuni.b, k1=k1, k2=k2, bcovs=ycov,
                     ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par
gestsunis <- optim(par=c(ps[2],ps[7:10]),f=logplikuni.g, k1=k1, k2=k2, gcovs=zcov,
                     zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
unisp <- c(gestsunis[4],gestsunis[1],bestsunis[1:3],gestsunis[c(2,3,5)])
######################################################
#univariate non-spatial
bestsunins <- optim(par=c(ipsf[4:5]),f=logplik.b.ns, k1=k1, k2=k2, bcovs=ycov,
                      ys=newys, nbrs=nbss,control=list(maxit=500), hessian=T)$par
gestsunins <- optim(par=ipsf[6:8],f=logplikuni.g.ns, k1=k1, k2=k2, gcovs=zcov,
                      zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
unins <- c(0,0,0,bestsunins[1:2],gestsunins[c(1,2,3)])
############################################################################
#############################################################################
#######Summarize the analysis for the report
#all of the estimates
wheatests <- data.frame(matrix(c(full,estsnc, unisp, unins, bns),nrow=5,ncol=8, byrow=T));
colnames(wheatests)<- c("rho","eZ","eY","b0y","b1y","b0z","b1z","sigsq")
rownames(wheatests)<- c("full", "no covariates", "univ. spatial", "univ. non-sp", "biv.ns")
wheatests
##############################################################
#all of the confidence intervals
CIub <- 2*wheatests-matrix(apply(wheatboots$ests,2,quantile, probs=c(.025)),nrow=5,ncol=8, byrow=T)
CIlb <- 2*wheatests-matrix(apply(wheatboots$ests,2,quantile, probs=c(.975)),nrow=5,ncol=8, byrow=T)
CIs <- data.frame(cbind(CIlb[,1:8],CIub[,1:8]),cbind(CIlb[,9:16],CIub[,9:16])); #rownames(CIs)<- rep(c("rho","eZ","eY","b0y","b1y","b0z","b1z","sigsq"),times=5)
###########################################################
# Continuous Ranked Probability Scores
fullbs<- brierscore(full[1:3],newys, newzs, nbss,
                            ycov, bcovpars=full[4:5], zcov, gcovpars=full[6:7], full[8])
fullgs<- GaussScore(full[1:3],newys, newzs, nbss,
                    ycov, bcovpars=full[4:5], zcov, gcovpars=full[6:7], full[8])
ncbs<- brierscore(estsnc[1:3],newys, newzs, nbss,
                    NULL, bcovpars=estsnc[4:5], NULL, gcovpars=estsnc[6:7], estsnc[8])
ncgs<- GaussScore(estsnc[1:3],newys, newzs, nbss,
                    NULL, bcovpars=estsnc[4:5], NULL, gcovpars=estsnc[6:7], estsnc[8])
bnbs<- brierscore(bns[1:3],newys, newzs, nbss,
                    ycov, bcovpars=bns[4:5], zcov, gcovpars=bns[6:7], bns[8])
bngs<- GaussScore(bns[1:3],newys, newzs, nbss,
                    ycov, bcovpars=bns[4:5], zcov, gcovpars=bns[6:7], bns[8])
unispbs<- brierscore(unisp[1:3],newys, newzs, nbss,
                    ycov, bcovpars=unisp[4:5], zcov, gcovpars=unisp[6:7], unisp[8])
unispgs<- GaussScore(unisp[1:3],newys, newzs, nbss,
                    ycov, bcovpars=unisp[4:5], zcov, gcovpars=unisp[6:7], unisp[8])
uninsbs<- brierscore(unins[1:3],newys, newzs, nbss,
                    ycov, bcovpars=unins[4:5], zcov, gcovpars=unins[6:7], unins[8])
uninsgs<- GaussScore(unins[1:3],newys, newzs, nbss,
                    ycov, bcovpars=unins[4:5], zcov, gcovpars=unins[6:7], unins[8])

wheatscore.b <- rbind(fullbs,ncbs,unispbs,uninsbs,bnbs)
wheatscore.g <- rbind(fullgs,ncgs,unispgs,uninsgs,bngs)
##############
modelests <- round(data.frame(wheatests, BinaryScore = wheatscore.b[,2],
                              BinarySkill=wheatscore.b[,1], GaussianScore = wheatscore.g[,2], GaussianSkill=wheatscore.g[,1]),4)
xtable(modelests)
















###################################################################
####################################################################
####################################################################
iout=1
bootstrapbivar <- function(spatials,bmeanpars, gmeanpars, sig2, startys, startzs, bcovt, gcovt, k1, k2, nbrs, B, M){
  ##Compute bootstrap standard errors 
  ##spatials/bmeanpars/gmeanpars - vector of estimated parameters
  ##startys - arbitrary vector of 0's and 1's
  ##startzs - arbitrary vector of gaussian rv
  ##B - number of iterations until burn-in
  ##M - number of bootstrap samples
  locs=1:(k1*k2);n=length(locs)
  estsfull <-matrix(nrow=M, ncol=8)
  colnames(estsfull)<- c("rho","eta.z","eta.y","b0.y","b1.y","b0.z","b1.z","sig.sq")
  uni.ns<-uni.sp <-ests.nocov<-ests.bns<-estsfull
  bscore <- array(dim=c(M,5,3)); gscore<-bscore
  
  cntzeros <- NULL
  cnt <- 1
  newys<-startys; newzs<-startzs
  repeat{
    dataset <- spatial.genfieldbivar.covs(newys, newzs,nbss,M=B,
                                          spatials,bmeanpars, gmeanpars, sigmasq=sig2, bcovs=bcovt, gcovs=gcovt)
    newys <- dataset[,1]
    newzs <- dataset[,2]
    ps <- initialparams.cov(newys,newzs,nbss,bcovt,gcovt)
    ests <- optim(par=ps,f=logplikbivar.covs, k1=k1, k2=k2,ys=newys,
                      zs=newzs, bcovs=bcovt,gcovs=gcovt, nbrs=nbss,hessian=T,control=list(maxit=1000))
    estsfull[cnt,] <- c(ests$par)
    bscore[cnt,1,]<- brierscore(estsfull[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=estsfull[cnt,4:5], gcovt, gcovpars=estsfull[cnt,6:7], estsfull[cnt,8])
    gscore[cnt,1,]<- GaussScore(estsfull[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=estsfull[cnt,4:5], gcovt, gcovpars=estsfull[cnt,6:7], estsfull[cnt,8])
    cnt<-cnt+1
    print(paste(cnt-1, "full model"))
    if(cnt==(M+1)){break}
  }
  cnt<-1
  newys<-startys; newzs<-startzs
  ests.nc <- optim(par=initialparams.cov(newys,newzs,nbss),f=logplikbivar.covs, k1=k1, k2=k2,ys=newys,
                zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
  ests.nc <-c(ests.nc[1:4],0,ests.nc[5],0,ests.nc[6])
  repeat{ 
    dataset <- spatial.genfieldbivar.covs(newys, newzs,nbss,M=B,
                                          ests.nc[1:3],ests.nc[4], ests.nc[6], sigmasq=ests.nc[8], bcovs=NULL, gcovs=NULL)
    newys <- dataset[,1]
    newzs <- dataset[,2]
    ests <- optim(par=initialparams.cov(newys,newzs,nbss),f=logplikbivar.covs, k1=k1, k2=k2,ys=newys,
                      zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    ests.nocov[cnt,] <-c(ests[1:4],0,ests[5],0,ests[6])
    bscore[cnt,2,]<- brierscore(ests.nocov[cnt,1:3],newys, newzs, nbss,
                                NULL, bcovpars=ests.nocov[cnt,4], NULL, gcovpars=ests.nocov[cnt,6], ests.nocov[cnt,8])
    gscore[cnt,2,]<- GaussScore(ests.nocov[cnt,1:3],newys, newzs, nbss,
                                NULL, bcovpars=ests.nocov[cnt,4], NULL, gcovpars=ests.nocov[cnt,6], ests.nocov[cnt,8])
    print(paste(cnt, ".nocovariates"))
    cnt<-cnt+1
    if(cnt==(M+1)){break}
  }
  
  cnt<-1
  newys<-startys; newzs<-startzs
  ps <- initialparams.cov(newys,newzs,nbss,bcovt,gcovt)
  ests.ns <- optim(par=c(ps[1], ps[4:8]),f=logplikbivar.ns, k1=k1, k2=k2, gcovs=gcovt, bcovs=bcovt,
                   zs=newzs, ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par
  ests.ns <-c(ests.ns[1],0,0,ests.ns[2:6])
  repeat{   
    dataset <- spatial.genfieldbivar.covs(newys, newzs,nbss,M=B,
                                          ests.ns[1:3],ests.ns[4:5], ests.ns[6:7], sigmasq=ests.ns[8], bcovs=bcovt, gcovs=gcovt)
    newys <- dataset[,1]
    newzs <- dataset[,2]
    ps <- initialparams.cov(newys,newzs,nbss,bcovt,gcovt)
    ests <- optim(par=c(ps[1], ps[4:8]),f=logplikbivar.ns, k1=k1, k2=k2, gcovs=gcovt, bcovs=bcovt,
                  zs=newzs, ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    ests.bns[cnt,] <-c(ests[1],0,0,ests[2:6])
    bscore[cnt,5,]<- brierscore(ests.bns[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=ests.bns[cnt,4:5], gcovt, gcovpars=ests.bns[cnt,6:7], ests.bns[cnt,8])
    gscore[cnt,5,]<- GaussScore(ests.bns[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=ests.bns[cnt,4:5], gcovt, gcovpars=ests.bns[cnt,6:7], ests.bns[cnt,8])
    cnt<-cnt+1
    print(paste(cnt-1, "bivariate non-spatial"))
    if(cnt==(M+1)){break}
  }
  
  cnt<-1
  newys<-startys; newzs<-startzs
  bcovg = cbind(bcovt, newzs-mean(newzs))
  gcovb = cbind(gcovt, newys-mean(newys))
  ps <- initialparams.cov(newys,newzs,nbss,bcovg,gcovb)
  bests.uni.s <- optim(par=c(ps[3:6]),f=logplikuni.b, k1=k1, k2=k2, bcovs=bcovg,
                       ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par
  gests.uni.s <- optim(par=c(ps[2],ps[7:10]),f=logplikuni.g, k1=k1, k2=k2, gcovs=gcovb,
                       zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
  uni.spp <- c(gests.uni.s[4],gests.uni.s[1],bests.uni.s[1:3],gests.uni.s[c(2,3,5)])
  repeat{   
    dataset <- spatial.genfieldbivar.covs(newys, newzs,nbss,M=B,
                                          uni.spp[1:3],uni.spp[4:5], uni.spp[6:7], sigmasq=uni.spp[8], bcovs=bcovt, gcovs=gcovt)
    newys <- dataset[,1]
    newzs <- dataset[,2]
    
    bcovg = cbind(bcovt, newzs-mean(newzs))
    gcovb = cbind(gcovt, newys-mean(newys))
    ps <- initialparams.cov(newys,newzs,nbss,bcovg,gcovb)
    bests.uni.s <- optim(par=c(ps[3:6]),f=logplikuni.b, k1=k1, k2=k2, bcovs=bcovg,
                         ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    gests.uni.s <- optim(par=c(ps[2],ps[7:10]),f=logplikuni.g, k1=k1, k2=k2, gcovs=gcovb,
                         zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    uni.sp[cnt,] <- c(gests.uni.s[4],gests.uni.s[1],bests.uni.s[1:3],gests.uni.s[c(2,3,5)])
    bscore[cnt,3,]<- brierscore(uni.sp[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=uni.sp[cnt,4:5], gcovt, gcovpars=uni.sp[cnt,6:7], uni.sp[cnt,8])
    gscore[cnt,3,]<- GaussScore(uni.sp[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=uni.sp[cnt,4:5], gcovt, gcovpars=uni.sp[cnt,6:7], uni.sp[cnt,8])
    cnt<-cnt+1
    print(paste(cnt-1, "univariate spatial"))
    if(cnt==(M+1)){break}
  }
  
  cnt<-1
  newys<-startys; newzs<-startzs
  bcovg = cbind(bcovt, newzs-mean(newzs))
  gcovb = cbind(gcovt, newys-mean(newys))
  ps <- initialparams.cov(newys,newzs,nbss,bcovt,gcovt)
  bests.uni.ns <- optim(par=c(ps[4:5]),f=logplik.b.ns, k1=k1, k2=k2, bcovs=bcovt,
                        ys=newys, nbrs=nbss,control=list(maxit=500), hessian=T)$par
  gests.uni.ns <- optim(par=ps[6:8],f=logplikuni.g.ns, k1=k1, k2=k2, gcovs=gcovt,
                        zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
  uni.nsn <- c(0,0,0,bests.uni.ns[1:2],gests.uni.ns[c(1:3)])
  repeat{   
    dataset <- spatial.genfieldbivar.covs(newys, newzs,nbss,M=B,
                                          uni.nsn[1:3],uni.nsn[4:5], uni.nsn[6:7], sigmasq=uni.nsn[8], bcovs=bcovt, gcovs=gcovt)
    newys <- dataset[,1]
    newzs <- dataset[,2]
    bcovg = cbind(bcovt, newzs-mean(newzs))
    gcovb = cbind(gcovt, newys-mean(newys))
    ps <- initialparams.cov(newys,newzs,nbss,bcovt,gcovt)
    bests.uni.ns <- optim(par=c(ps[4:5]),f=logplik.b.ns, k1=k1, k2=k2, bcovs=bcovt,
                          ys=newys, nbrs=nbss,control=list(maxit=500), hessian=T)$par
    gests.uni.ns <- optim(par=ps[6:8],f=logplikuni.g.ns, k1=k1, k2=k2, gcovs=gcovt,
                          zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    uni.ns[cnt,] <- c(0,0,0,bests.uni.ns[1:2],gests.uni.ns[c(1,2,3)])
    bscore[cnt,4,]<- brierscore(uni.ns[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=uni.ns [cnt,4:5], gcovt, gcovpars=uni.ns[cnt,6:7], uni.ns [cnt,8])
    gscore[cnt,4,]<- GaussScore(uni.ns[cnt,1:3],newys, newzs, nbss,
                                bcovt, bcovpars=uni.ns [cnt,4:5], gcovt, gcovpars=uni.ns [cnt,6:7], uni.ns [cnt,8])
    cnt<-cnt+1
    print(paste(cnt-1, "univariate non-spatial"))
    if(cnt==(M+1)){break}
  } 
      ests<-cbind(estsfull,ests.nocov, uni.sp, uni.ns, ests.bns)

  colnames(gscore)<- colnames(bscore)<-c("Full","No Covariates", "Uni", "UniNS", "BivSp")
  return(list(ests = ests, CRPSg = gscore, CRPSb = bscore))
}
##################################################
##################################################



