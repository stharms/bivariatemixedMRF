library(tidyverse);library(viridisLite); library(scoringRules)
source("ReducedBivariateLLPs.R")
source("ScoringRules.R")
source("Sstatistic.R");
source("BinaryUpperBounds.R")
source("ModelwCovs.R")

k1=30;k2=30;kk=k1*k2
NBHD4=rectanggrid.4nbrs(k1,k2)
locs=1:kk

#mean parameters
sig2 = 1
bmeans <- c(-1,0.5); gmeans <- c(1,0.1); spatpars <- c(0.2,0.3,1) 
truths=c(spatpars,bmeans,gmeans,sig2)

#generate covariates
#not run
##################################3
xmult <- 0.5*seq(1,k1,by=1)
ymult <- 0.0005*seq(1,k2,by=1)^2
s <-matrix(1,nrow=k1,ncol=k2)
s1 <- t(s*xmult)*ymult
covb<- spatial.genfieldbivar.covs(ys,zs,NBHD4,M=150,
                                  c(0,0.9,5),c(-0.5,1), c(1,-0.5), sigmasq=sig2, bcovs=c(s1), gcovs=c(s1)); colMeans(covg);
gcov <- rgamma(n=k1*k2, shape=3, scale=4)
######################################################
bcovt <- read.csv(file='ycovariates.csv', header=T)[,2]
gcovt <- read.csv(file='zcovariates.csv', header=T)[,2]
image(matrix(bcovt,nrow=30), col=magma(20))
image(matrix(gcovt,nrow=30), col=magma(20))

#generate the datasets and estimate the parameters for each parameter set
spatpars <- c(0.2,0.3,1) 
weakests <- bootstrapbivar(spatpars, bmeans, gmeans, sig2, ys, zs, k1=30, k2=30, nbss, B=300, M=1000, S=50)
write.csv(modests, file = "weakests3.csv")
spatpars <- c(1,0.3,1)
modests <- bootstrapbivar(spatpars, bmeans, gmeans, sig2, ys, zs, k1=30, k2=30, nbss, B=300, M=1000, S=30)
write.csv(modests, file = "modests3.csv")
spatpars <- c(0.5,0.9,3.5)
strongests <- bootstrapbivar(spatpars, bmeans, gmeans, sig2, ys, zs, k1=30, k2=30, nbss, B=300, M=1000, S=30)
write.csv(strongests, file = "strongests3.csv")
####################################################################
####################################################################
###simulation function##############################################
#####################################################################
iout=1
bootstrapbivar <- function(spatials,bmeanpars, gmeanpars, sig2, startys, startzs, k1, k2, nbrs, B, M,S){
  ##Compute bootstrap standard errors 
  ##spatials/bmeanpars/gmeanpars - vector of estimated parameters
  ##startys - arbitrary vector of 0's and 1's
  ##startzs - arbitrary vector of gaussian rv
  ##B - number of iterations until burn-in
  ##M - number of bootstrap samples
  #generate dataset
  burnin <- spatial.genfieldbivar.covs(startys,startzs,nbss,M=B,
                                        spatials,bmeanpars, gmeanpars, sigmasq=sig2, bcovs=bcovt, gcovs=gcovt)
  newys <- burnin[,1]
  newzs <- burnin[,2]
  locs=1:(k1*k2);n=length(locs)
  ests <- matrix(nrow=M, ncol=8*5)
  cntzeros <- NULL
  cnt <- 1
  repeat{
    #generate dataset
    dataset <- spatial.genfieldbivar.covs(newys,newzs,nbss,M=S,
                                          spatials,bmeanpars, gmeanpars, sigmasq=sig2, bcovs=bcovt, gcovs=gcovt)
    newys <- dataset[,1]
    newzs <- dataset[,2]
    #estimate all of the models
    ps <- initialparams.cov(newys,newzs,nbss,bcovt,gcovt)
    estsfull <- optim(par=ps,f=logplikbivar.covs, k1=k1, k2=k2,ys=newys,
                      zs=newzs, bcovs=bcovt,gcovs=gcovt, nbrs=nbss,hessian=T,control=list(maxit=1000))
    ests.nocov <- optim(par=initialparams.cov(newys,newzs,nbss),f=logplikbivar.covs, k1=k1, k2=k2,ys=newys,
                        zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    ests.ns <- optim(par=c(ps[1], ps[4:8]),f=logplikbivar.ns, k1=k1, k2=k2, gcovs=gcovt, bcovs=bcovt,
                     zs=newzs, ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    #new covariates for the univariate models
    bcovg = cbind(bcovt, newzs-mean(newzs))
    gcovb = cbind(gcovt, newys-mean(newys))
    ps <- initialparams.cov(newys,newzs,nbss,bcovg,gcovb)
    bests.uni.s <- optim(par=c(ps[3:6]),f=logplikuni.b, k1=k1, k2=k2, bcovs=bcovt,
                         ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    gests.uni.s <- optim(par=c(ps[2],ps[7:10]),f=logplikuni.g, k1=k1, k2=k2, gcovs=gcovt,
                         zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    uni.sp <- c(gests.uni.s[4],gests.uni.s[1],bests.uni.s[1:3],gests.uni.s[c(2,3,5)])
    
    bests.uni.ns <- optim(par=c(ps[3:6]),f=logplik.b.ns, k1=k1, k2=k2, bcovs=bcovt,
                          ys=newys, nbrs=nbss,control=list(maxit=500), hessian=T)$par
    gests.uni.ns <- optim(par=ps[7:10],f=logplikuni.g.ns, k1=k1, k2=k2, gcovs=gcovt,
                          zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
    uni.ns <- c(gests.uni.ns[3],0,0,bests.uni.ns[1:2],gests.uni.ns[c(1,2,4)])
    
    ests[cnt,] <- c(estsfull$par, c(ests.nocov[1:4],0,ests.nocov[5],0,ests.nocov[6]),
                    uni.sp, uni.ns, c(ests.ns[1],0,0,ests.ns[2:6]))

    if(estsfull$conv!=0){cntzeros <- c(1, cntzeros)}
    datasetys <- newys
    datasetzs <- newzs
    print(paste(cnt, iout))
    cnt <- cnt+1
    if(cnt==(M+1)){break}
  }

  return(ests)
}


##################################################
##################################################

#Creating Images and Boxplots of the results


library(extrafont)
loadfonts(device = "win")
par(family = "CM Roman Greek")

bcovariates <- read.csv(file='ycovariates.csv', header=T)[,2]
bcovmat <- matrix(bcovariates, nrow=30,ncol=30)
image(matrix(bcovariates, nrow=30,ncol=30), col=inferno(20), axes=F,
      useRaster=T,main = expression('Binary Covariate X'['y']))
pdf("ycovsim.pdf")
dev.off()

gcovariates <- read.csv(file='zcovariates.csv', header=T)[,2]
image(matrix(gcovariates, nrow=30,ncol=30), col=viridis(20), axes=F, main = expression('Gaussian Covariate X'['z']))
pdf("zcovsim.pdf")
dev.off()

weakests <- read.csv(file='weakests3.csv', header=T)[,-1]
colnames(weakests)<- paste(rep(c("full","nocovs","uni.sp","uni.ns","biv.ns"),each=8),
                           rep(c("rho", "eta_z", "eta_y", "b0.y","b1.y","b0.z","b1.z","sigma2"),times=5))
weakmeans <- round(apply(weakests, 2, mean),4)
truths=c(.2,.3,1,-1,0.5,1,0.1,1)
sweeperpars <- rep(truths, times=5)
weakbias <- round(apply(sweep(weakests, 2, sweeperpars, FUN="-"),2,mean),4)
weakpctbias <- round(100*weakbias/sweeperpars,5)
weakses <- round(apply(weakests, 2, sd),3)
weaklows <- round(apply(weakests, 2, quantile, probs=.05),4)
weakhighs<- round(apply(weakests, 2, quantile, probs=.95),4)
cbind(paste(weakmeans, "(",weaklows,",",weakhighs,")"))
#data.frame(matrix(paste(weakmeans, "(",weaklows,",",weakhighs,")"), nrow=8,ncol=5))
#weaktable<-data.frame(matrix(paste(weakmeans, " (",weaklows," , ",weakhighs,")", sep=""), nrow=5,ncol=8,byrow=T))
#weaktablemeans <- data.frame(matrix(paste(weakmeans), nrow=5,ncol=8,byrow=T))
#weaktableCIs <- data.frame(matrix(paste(" (",weaklows," , ",weakhighs,")", sep=""), nrow=5,ncol=8,byrow=T))
#xtable(weaktablemeans)
#xtable(weaktableCIs)
cbind(paste(weakbias, "(",weakpctbias,"%)", "(",weakses,")"))
data.frame(matrix(paste0(weakbias, "(",weakpctbias,"%)", "(",weakses,")"), nrow=8,ncol=5))
weaktable<-data.frame(matrix(paste0(weakbias, " (",weakpctbias,"%)", "(",weakses,")"), nrow=5,ncol=8,byrow=T))
weaktablebias <- data.frame(matrix(paste0(weakbias, " [",weakpctbias,"%]"), nrow=5,ncol=8,byrow=T))
weaktableSEs <- data.frame(matrix(paste0("(",weakses,")"), nrow=5,ncol=8,byrow=T))
xtable(weaktablebias)
xtable(weaktableSEs)

modests <- read.csv(file='modests3.csv', header=T)[,-1]
colnames(modests)<- paste(rep(c("full","nocovs","uni.sp","uni.ns","biv.ns"),each=8),
                           rep(c("rho", "eta_z", "eta_y", "b0.y","b1.y","b0.z","b1.z","sigma2"),times=5))
modmeans <- round(apply(modests, 2, mean),4)
modlows <- round(apply(modests, 2, quantile, probs=.05),4)
modhighs<- round(apply(modests, 2, quantile, probs=.95),4)
modtable<-data.frame(matrix(paste(modmeans, " (",modlows," , ",modhighs,")", sep=""), nrow=5,ncol=8,byrow=T))
truths = c(1,0.3,1,-1,0.5,1,0.1,1)
sweeperpars <- rep(truths, times=5)
modbias <- round(apply(sweep(modests, 2, sweeperpars, FUN="-"),2,mean),3)
modpctbias <- round(100*modbias/sweeperpars,2)
modses <- round(apply(modests, 2, sd),3)
xtable(modtable)
#modtablemeans <- data.frame(matrix(paste(modmeans), nrow=5,ncol=8,byrow=T))
#modtableCIs <- data.frame(matrix(paste(" (",modlows," , ",modhighs,")", sep=""), nrow=5,ncol=8,byrow=T))
#xtable(modtablemeans)
#xtable(modtableCIs)
#xtable(modtable)
modtable<-data.frame(matrix(paste0(modbias, " [",modpctbias,"%]", "(",modses,")"), nrow=5,ncol=8,byrow=T))
modtablebias <- data.frame(matrix(paste0(modbias, " [",modpctbias,"%]"), nrow=5,ncol=8,byrow=T))
modtableSEs <- data.frame(matrix(paste0("(",modses,")"), nrow=5,ncol=8,byrow=T))
xtable(modtablebias)
xtable(modtableSEs)

strongests <- read.csv(file='strongests3.csv', header=T)[,-1]
colnames(strongests)<- paste(rep(c("full","nocovs","uni.sp","uni.ns","biv.ns"),each=8),
                          rep(c("rho", "eta_z", "eta_y", "b0.y","b1.y","b0.z","b1.z","sigma2"),times=5))
truths = c(0.5,0.9,3.5,-1,0.5,1,0.1,1)
sweeperpars <- rep(truths, times=5)
strongbias <- round(apply(sweep(strongests, 2, sweeperpars, FUN="-"),2,mean),3)
strongpctbias <- round(100*strongbias/sweeperpars,2)
strongses <- round(apply(strongests, 2, sd),3)
xtable(strongtable)
strongmeans <- round(apply(strongests, 2, mean),4)
stronglows <- round(apply(strongests, 2, quantile, probs=.05),4)
stronghighs<- round(apply(strongests, 2, quantile, probs=.95),4)
#strongtable<-data.frame(matrix(paste(strongmeans, " (",stronglows," , ",stronghighs,")",sep=""), nrow=5,ncol=8,byrow=T))
#xtable(strongtable)
#strongtablemeans <- data.frame(matrix(paste(strongmeans), nrow=5,ncol=8,byrow=T))
#strongtableCIs <- data.frame(matrix(paste(" (",stronglows," , ",stronghighs,")", sep=""), nrow=5,ncol=8,byrow=T))
#xtable(strongtablemeans)
#xtable(strongtableCIs)

strongtable<-data.frame(matrix(paste0(strongbias, " [",strongpctbias,"%]", "(",strongses,")"), nrow=5,ncol=8,byrow=T))
strongtablebias <- data.frame(matrix(paste0(strongbias, " [",strongpctbias,"%]"), nrow=5,ncol=8,byrow=T))
strongtableSEs <- data.frame(matrix(paste0("(",strongses,")"), nrow=5,ncol=8,byrow=T))
xtable(strongtablebias)
xtable(strongtableSEs)
####################################
#plots

plotform <- function(parmat){
  colnames(parmat)<-rep(c("rho", "eta_z", "eta_y", "b0.y","b1.y","b0.z","b1.z","sigma2"),times=5)
  #outmat<-matrix(ncol=3); colnames(outmat)<-c("parameter","value", "model")
  modelno <- rep((1:5), each=8000)
  spreadmat <- stack(parmat)
  parnames <- rep(rep(c("rho", "eta_z", "eta_y", "b0.y","b1.y","b0.z","b1.z","sigma2"),each=1000),times=5)
  spreadmat <- cbind(modelno,spreadmat, parnames)[,-3]
  colnames(spreadmat)<- c("modelno", "values","ind")
  return(spreadmat)
}

library(gridExtra)
weakplot <- plotform(weakests); modplot<-plotform(modests); strongplot <- plotform(strongests)
weakplot <- weakplot %>% mutate(values = case_when((ind%in%c("eta_z","eta_y") & modelno >= 4)~NA_real_,TRUE~values))%>% 
  mutate(values = case_when((ind %in% c("b1.y","b1.z") & modelno ==2)~NA_real_,TRUE~values)) 
modplot <- modplot %>% mutate(values = case_when((ind%in%c("eta_z","eta_y") & modelno >= 4)~NA_real_,TRUE~values)) %>% 
  mutate(values = case_when((ind %in% c("b1.y","b1.z") & modelno ==2)~NA_real_,TRUE~values)) 
strongplot <- strongplot %>% mutate(values = case_when((ind%in%c("eta_z","eta_y") & modelno >= 4)~NA_real_,TRUE~values))%>% 
  mutate(values = case_when((ind %in% c("b1.y","b1.z") & modelno ==2)~NA_real_,TRUE~values)) 
summary(strongplot[16001:17000,])

plot.r.w<- ggplot(filter(weakplot,ind=="rho"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.2, linetype="dashed", size=1.5) + labs(x="", y=expression(rho)) + ggtitle("Weak")+
  ylim(0,.75) + theme(axis.title.y = element_text(size=24, angle=0, vjust=0.5) , legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.ez.w<- ggplot(filter(weakplot,ind=="eta_z"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.3, linetype="dashed", size=1.5)+ labs(x="", y=expression(eta[z])) + ggtitle("Weak")+ylim(0,0.6)+
  theme(axis.title.y = element_text(size=24, angle=0, vjust=0.5), legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.ey.w<- ggplot(filter(weakplot,ind=="eta_y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5)+ labs(x="", y=expression(eta[y])) + ggtitle("Weak")+ylim(-0.5,2.5)+
  theme(axis.title.y = element_text(size=24, angle=0, vjust=0.5), legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b0y.w<- ggplot(filter(weakplot,ind=="b0.y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=-1, linetype="dashed", size=1.5)+ labs(x="", y=expression(beta["0y"])) + ggtitle("Weak")+ylim(-1.4,-0.5)+
  theme(axis.title.y = element_text(size=24, angle=0, vjust=0.5), legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b1y.w<- ggplot(filter(weakplot,ind=="b1.y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.5, linetype="dashed", size=1.5)+ labs(x="", y=expression(beta["1y"])) + ggtitle("Weak")+ylim(0.25,0.7)+
  theme(axis.title.y = element_text(size=24, angle=0, vjust=0.5), legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b0z.w<- ggplot(filter(weakplot,ind=="b0.z"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5)+ labs(x="", y=expression(beta["0z"])) + ggtitle("Weak")+ylim(0.7,1.3)+
  theme(axis.title.y = element_text(size=24, angle=0, vjust=0.5), legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b1z.w<- ggplot(filter(weakplot,ind=="b1.z"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.1, linetype="dashed", size=1.5)+ labs(x="", y=expression(beta["1z"])) + ggtitle("Weak")+ ylim(0.07,0.13)+
  theme(axis.title.y = element_text(size=24, angle=0, vjust=0.5), legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.sig2.w<- ggplot(filter(weakplot,ind=="sigma2"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5)+ labs(x="", y=expression(sigma^2)) + ggtitle("Weak")+ylim(0,1.5)+
  theme(axis.title.y = element_text(size=24, angle=0, vjust=0.5), legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
gw <-grid.arrange(plot.r.w, plot.ez.w,plot.ey.w,plot.b0y.w,plot.b1y.w,plot.b0z.w,plot.b1z.w,plot.sig2.w,
             nrow=8,ncol=1)
gws <-grid.arrange(plot.r.w, plot.ez.w,plot.ey.w,nrow=3,ncol=1)
gwc <-grid.arrange(plot.b0y.w,plot.b1y.w,plot.b0z.w,plot.b1z.w,plot.sig2.w,nrow=5,ncol=1)
gwc2 <-grid.arrange(plot.b0y.w,plot.b1y.w,plot.b0z.w,plot.b1z.w,nrow=4,ncol=1)

plot.r.m<- ggplot(filter(modplot,ind=="rho"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5) + labs(x="", y="") + ggtitle("Moderate")+ ylim(0,2)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.ez.m<- ggplot(filter(modplot,ind=="eta_z"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.3, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Moderate")+ylim(0,0.6)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.ey.m<- ggplot(filter(modplot,ind=="eta_y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Moderate")+ylim(-0.5,3)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b0y.m<- ggplot(filter(modplot,ind=="b0.y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=-1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Moderate")+ylim(-1.5,-0.3)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b1y.m<- ggplot(filter(modplot,ind=="b1.y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.5, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Moderate")+ylim(0.2,0.7)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b0z.m<- ggplot(filter(modplot,ind=="b0.z"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Moderate") +ylim(0.7,1.4)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b1z.m<- ggplot(filter(modplot,ind=="b1.z",as.factor(modelno)!=2), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Moderate")+ ylim(0.07,0.13)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.sig2.m<- ggplot(filter(modplot,ind=="sigma2"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Moderate")+ylim(0,1.5)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
gm <-grid.arrange(plot.r.m, plot.ez.m,plot.ey.m,plot.b0y.m,plot.b1y.m,plot.b0z.m,plot.b1z.m,plot.sig2.m,
                  nrow=8,ncol=1)
gms <-grid.arrange(plot.r.m, plot.ez.m,plot.ey.m,nrow=3,ncol=1)
gmc <-grid.arrange(plot.b0y.m,plot.b1y.m,plot.b0z.m,plot.b1z.m,plot.sig2.m,nrow=5,ncol=1)
gmc2 <-grid.arrange(plot.b0y.m,plot.b1y.m,plot.b0z.m,plot.b1z.m,nrow=4,ncol=1)

plot.r.s<- ggplot(filter(strongplot,ind=="rho"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.5, linetype="dashed", size=1.5) + labs(x="", y="") + ggtitle("Strong")+ ylim(-1,2)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.ez.s<- ggplot(filter(strongplot,ind=="eta_z"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.9, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Strong") + ylim(.5,1.1)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.ey.s<- ggplot(filter(strongplot,ind=="eta_y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=3.5, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Strong")+ ylim(0.8,7)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b0y.s<- ggplot(filter(strongplot,ind=="b0.y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=-1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Strong")+ ylim(1,3.5)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b1y.s<- ggplot(filter(strongplot,ind=="b1.y"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.5, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Strong")+ ylim(0,1)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b0z.s<- ggplot(filter(strongplot,ind=="b0.z"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Strong")+ ylim(0,6)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.b1z.s<- ggplot(filter(strongplot,ind=="b1.z"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=0.1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Strong")+ ylim(0.07,0.13)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))
plot.sig2.s<- ggplot(filter(strongplot,ind=="sigma2"), aes(x=as.factor(modelno)))+geom_boxplot(aes(y=values, fill=(as.factor(modelno))))+
  geom_hline(yintercept=1, linetype="dashed", size=1.5)+ labs(x="", y="") + ggtitle("Strong")+ ylim(0,2.5)+
  theme(legend.position="none") + scale_fill_manual(values = c("blue","white","white","white","white"))

gs<-grid.arrange(plot.r.s, plot.ez.s,plot.ey.s,plot.b0y.s,plot.b1y.s,plot.b0z.s,plot.b1z.s,plot.sig2.s,
             nrow=8,ncol=1)
gss <-grid.arrange(plot.r.s, plot.ez.s,plot.ey.s,nrow=3,ncol=1)
gsc <-grid.arrange(plot.b0y.s,plot.b1y.s,plot.b0z.s,plot.b1z.s,plot.sig2.s,nrow=5,ncol=1)
gsc2 <-grid.arrange(plot.b0y.s,plot.b1y.s,plot.b0z.s,plot.b1z.s,nrow=4,ncol=1)

grid.arrange(gws,gms,gss,ncol=3)
grid.arrange(gwc,gmc,gsc,ncol=3)
grid.arrange(gwc2,gmc2,gsc2,ncol=3)
gr<-grid.arrange(plot.r.w,plot.r.m,plot.r.s,ncol=3)
gez<-grid.arrange(plot.ez.w,plot.ez.m,plot.ez.s,ncol=3)
gey<-grid.arrange(plot.ey.w,plot.ey.m,plot.ey.s,ncol=3)
grid.arrange(gr,gez,gey,nrow=3)
pdf("spatialests.pdf")

gb0y<-grid.arrange(plot.b0y.w,plot.b0y.m,plot.b0y.s,ncol=3)
gb1y<-grid.arrange(plot.b1y.w,plot.b1y.m,plot.b1y.s,ncol=3)
gb0z<-grid.arrange(plot.b0z.w,plot.b0z.m,plot.b0z.s,ncol=3)
gb1z<-grid.arrange(plot.b1z.w,plot.b1z.m,plot.b1z.s,ncol=3)
grid.arrange(gb0y,gb1y,gb0z,gb1z,nrow=4)
grid.arrange(gb0y,gb1y,nrow=2)
pdf("ymeanests.pdf")

grid.arrange(gb0z,gb1z,nrow=2)
pdf("zmeanests.pdf")
dev.off()
