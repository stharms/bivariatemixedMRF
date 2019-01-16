library(tidyverse);library(viridisLite); library(scoringRules); library(ggplot2)
source("ReducedBivariateLLPs.R")
source("Sstatistic.R");
source("BinaryUpperBounds.R")
source("ModelwCovs.R")
#####################################################
#################################################

#Scaling Test
sig2 = 1
begvals <- c(1, 0.5, 2, .5, 1 ,1)

testparams <- begvals
testsig2s <- seq(from=1, to=100,by=5)
nsimeach = 50; nsim1 <- 2; nsim2 <- length(testrhos); nsim3 <- length(testeta2s)

iout <- 0; jout <- 0; kout <- 0; lout <- 0;
nbd <- rectanggrid.4nbrs(30,30)
k1=30;k2=30
bout <-0
ymean <- array(dim=c(nsimeach, length(testsig2s), 2)) ; zmean<-array(dim=c(nsimeach, length(testsig2s), 2))
SB <- array(dim=c(nsimeach, length(testsig2s), 2)); SG<-array(dim=c(nsimeach, length(testsig2s), 2))
repeat{
iout<-iout+1
jout<-0
repeat{
      jout <- jout + 1
      testparams[6] <- testsig2s[jout]
      testparams[2] <- begvals[2]
      testparams[1]<-begvals[1]
      if(iout==2){testparams[1]<-begvals[1]*sqrt(testparams[6])}
      print(paste("sig2s: ", testparams[6]))
      kout <- 0
      repeat{
        kout <- kout+1
        print(paste(testparams[6],".",jout, " iteration: ", kout,".",iout))
        yss <- rbinom(prob = exp(testparams[4])/(1+exp(testparams[4])), size = 1, n = length(nbd[,1])) 
        zss <- rnorm(length(nbd[,1]), mean = testparams[5], sd = sqrt(testparams[6]))
        newf <- spatial.genfieldbivar.covs(yss,zss,nbd,M=300,
                                          testparams[1:3],testparams[4],testparams[5], sigmasq=testparams[6], bcovs=NULL, gcovs=NULL)
        ipd <- initialparams.cov(newf[,1], newf[,2], nbd)
        #epd <- try(optim(par=ipd,f=f=logplikbivar.covs, k1=k1, k2=k2,ys=newf[,1],
                           #zs=newf[,2], nbrs=nbd,hessian=T,control=list(maxit=500)))
        SB[kout,jout,iout] <- unname(Dbinary(newf[,1],nbd, k1=k1,k2=k2))[1]
        SG[kout, jout,iout] <- unname(Dgaussian(newf[,2],nbd,nbins=k2,sigma2=var(newf[,2])))[2]*testparams[6]
        ymean[kout,jout,iout] <- mean(newf[,1])
        zmean[kout,jout,iout] <- mean(newf[,2])
        if(kout==nsimeach){break}
      }
      #print(paste(jout))
      if(jout==length(testsig2s)){break}
}
  if(iout==2){break}
}


###################
sig2 = 10000
begvals <- c(100, 0, 1, 0.5, 1 ,sig2)

testparams <- begvals
testdeltas <- seq(0,4,by=.2)
testeta2s <- seq(0,2,by=.2); 
nsimeach = 50; nsim1 <- 2; nsim2 <- length(testrhos); nsim3 <- length(testeta2s)
nsimeach*nsim1*nsim2*nsim3


rhoconst <- sqrt((1-begvals[2])*(ubest(begvals[4])-begvals[3])); rhoconst

iout <- 0; jout <- 0; kout <- 0; lout <- 0;
nbd <- rectanggrid.4nbrs(30,30)
k1=30;k2=30
bout <-0
ymeand <- array(dim=c(nsimeach, length(testeta2s), 2)) ; zmeand<-array(dim=c(nsimeach, length(testeta2s), 2))
SBd <- array(dim=c(nsimeach, length(testeta2s), 2)); SGd<-array(dim=c(nsimeach, length(testeta2s), 2))
repeat{
  iout<-iout+1
  jout<-0
  repeat{
    jout <- jout + 1
    testparams[4] <- testeta2s[jout]
    testparams[1]<-begvals[1]
    if(iout==2){testparams[1]<-begvals[1]*sqrt((1-testparams[2])*(ubest(testparams[4])-testparams[3]))/rhoconst}
    print(paste("delta: ", testparams[4]))
    kout <- 0
    repeat{
      kout <- kout+1
      print(paste(testparams[4],".",jout, " iteration: ", kout,".",iout))
      yss <- rbinom(prob = exp(testparams[4])/(1+exp(testparams[4])), size = 1, n = length(nbd[,1])) 
      zss <- rnorm(length(nbd[,1]), mean = testparams[5], sd = sqrt(testparams[6]))
      newf <- spatial.genfieldbivar.covs(yss,zss,nbd,M=300,
                                         testparams[1:3],testparams[4],testparams[5], sigmasq=testparams[6], bcovs=NULL, gcovs=NULL)
      ipd <- initialparams.cov(newf[,1], newf[,2], nbd)

      SBd[kout,jout,iout] <- unname(Dbinary(newf[,1],nbd, k1=k1,k2=k2))[1]
      SGd[kout, jout,iout] <- unname(Dgaussian(newf[,2],nbd,nbins=k2,sigma2=var(newf[,2])))[2]
      ymeand[kout,jout,iout] <- mean(newf[,1])
      zmeand[kout,jout,iout] <- mean(newf[,2])
      if(kout==nsimeach){break}
    }
    if(jout==length(testeta2s)){break}
  }
  if(iout==2){break}
}
#####################
#Transmission
sig2 = 1
begvals <- c(1, 0.5, 1.5, 0.5, 1 ,1)

testparams <- begvals
spatials1<- c(0.7,0); spatials2<-c(0,2.5)
testrhos <- seq(0,1.2, by=0.1)
nsimeach = 50; nsim1 <- 2; nsim2 <- length(testrhos);
nsimeach*nsim1*nsim2


rhoconst <- 1/sqrt((1-begvals[6]*begvals[2])*(ubest(begvals[4])-begvals[3])); rhoconst

iout <- 0; jout <- 0; kout <- 0; lout <- 0;
nbd <- rectanggrid.4nbrs(30,30)
k1=30;k2=30
bout <-0
ymeanT <- array(dim=c(nsimeach, length(testrhos), 2)) ; zmeanT<-array(dim=c(nsimeach, length(testrhos), 2))
SBT <- array(dim=c(nsimeach, length(testrhos), 2)); SGT<-array(dim=c(nsimeach, length(testrhos), 2))
repeat{
  iout<-iout+1
  if(iout==1){testparams[c(2,3)]<-spatials1}
  if(iout==2){testparams[c(2,3)]<-spatials2}
  jout<-0
  repeat{
    jout <- jout + 1
    testparams[1] <- testrhos[jout]
    print(paste("rhos: ", testparams[1]))
    kout <- 0
    repeat{
      kout <- kout+1
      print(paste(testparams[1],".",jout, " iteration: ", kout,".",iout))
      yss <- rbinom(prob = exp(testparams[4])/(1+exp(testparams[4])), size = 1, n = length(nbd[,1])) 
      zss <- rnorm(length(nbd[,1]), mean = testparams[5], sd = sqrt(testparams[6]))
      newf <- spatial.genfieldbivar.covs(yss,zss,nbd,M=300,
                                         testparams[1:3],testparams[4],testparams[5], sigmasq=testparams[6], bcovs=NULL, gcovs=NULL)
      ipd <- initialparams.cov(newf[,1], newf[,2], nbd)
      #epd <- try(optim(par=ipd,f=f=logplikbivar.covs, k1=k1, k2=k2,ys=newf[,1],
      #zs=newf[,2], nbrs=nbd,hessian=T,control=list(maxit=500)))
      SBT[kout,jout,iout] <- unname(Dbinary(newf[,1],nbd, k1=k1,k2=k2))[1]
      SGT[kout, jout,iout] <- unname(Dgaussian(newf[,2],nbd,nbins=k2,sigma2=var(newf[,2])))[2]
      ymeanT[kout,jout,iout] <- mean(newf[,1])
      zmeanT[kout,jout,iout] <- mean(newf[,2])
      if(kout==nsimeach){break}
    }
    #print(paste(jout))
    if(jout==length(testrhos)){break}
  }
  if(iout==2){break}
}

#####################
ymeans <- apply(ymean, c(2,3), median)
zmeans <- apply(zmean, c(2,3), mean)
SYs <- apply(SB, c(2,3), mean)
SZs <- apply(SG, c(2,3), mean)
sigtests.1 <- data.frame(cbind(testsig2s,ymeans[,1],zmeans[,1],SYs[,1],SZs[,1]), scale = "rho/sigmasq")
sigtests.2 <- data.frame(cbind(testsig2s,ymeans[,2],zmeans[,2],SYs[,2],SZs[,2]), scale = "rho/sqrt(sigmasq)")
colnames(sigtests.2)<-colnames(sigtests.1)<- c("sig2","ymean","zmean","sY","sZ", "scale") 
sigtests <- rbind(sigtests.1,sigtests.2)

plot(data=sigtests.1,SYs[,1]~testsig2s, type="l")+lines(data=sigtests.2, zmean~sig2, col="red")
g <- ggplot(data=sigtests)
g + geom_smooth(mapping = aes(y=ymean , x=sig2), se = F) + geom_smooth(data=sigtests.2, mapping = aes(y=ymean , x=sig2 ),colour="black", se = F)
g + geom_smooth(mapping = aes(y=sY , x=sig2), se = F) + geom_smooth(data=sigtests.2, mapping = aes(y=sY,x=sig2),colour="black", se = F)
g+ geom_smooth(mapping = aes(y=sY , x=sig2, linetype=factor(scale)),size=1.5, se = F, colour="black") + #ylim(0,2.2)+
  labs(x=expression(sigma^2), y="Binary S-value")+
  scale_linetype_discrete(name="Scale", labels = c(expression(rho),expression("     "*sqrt(sigma^2)*rho)))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=16, face = "bold"), 
        legend.text=element_text(size=16), legend.key.size = unit(1, "cm"),
        axis.title.y = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text = element_text(size=16))
pdf("scaletestSY.pdf")

g+ geom_smooth(mapping = aes(y=sY , x=sig2, linetype=factor(scale)), se = F, colour="black") +
  labs(x=expression(sigma^2), y="Gaussian S-value")+
  scale_linetype_discrete(name="Scale", labels = c(expression(rho),expression("     "*sqrt(sigma^2)*rho)))+
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12), legend.key.size = unit(1, "cm"))

g+ geom_smooth(mapping = aes(y=ymean , x=sig2, linetype=factor(scale)), se = F, colour="black") +
  labs(x=expression(sigma^2), y="Binary Mean")+
  scale_linetype_discrete(name="Scale", labels = c(expression(rho),expression("     "*sqrt(sigma^2)*rho)))+ ylim(0.5,0.65)+
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=16, face = "bold"), 
        legend.text=element_text(size=16), legend.key.size = unit(1, "cm"),
        axis.title.y = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text = element_text(size=16))+
  geom_hline(yintercept = exp(0.5)/(1+exp(0.5)), linetype=6)
pdf("scaletestYmean.pdf")

#######################
ymeansT <- apply(ymeanT, c(2,3), mean)
zmeansT <- apply(zmeanT, c(2,3), mean)
SYsT <- apply(SBT, c(2,3), mean)
SZsT <- apply(SGT, c(2,3), mean)
Ttests.1 <- data.frame(cbind(seq(0,1.2, by=0.1),ymeansT[,1],zmeansT[,1],SYsT[,1],SZsT[,1]), set = "(0.7 , 0)")
Ttests.2 <- data.frame(cbind(seq(0,1.2, by=0.1),ymeansT[,2],zmeansT[,2],SYsT[,2],SZsT[,2]), set = "(0 , 2.5)")
colnames(Ttests.2)<-colnames(Ttests.1)<- c("rho","ymean","zmean","sY","sZ", "set") ;
Ttests <- rbind(Ttests.1,Ttests.2)

plot(data=Ttests.1,sZ~rho, type="l")+lines(data=Ttests.2, sZ~rho, col="red")
h <- ggplot(data=Ttests.1)
h + geom_smooth(mapping = aes(y=ymean , x=rho), se = F) + geom_smooth(data=Ttests.2, mapping = aes(y=ymean , x=rho ),colour="black", se = F)
h + geom_smooth(mapping = aes(y=sY , x=rho), se = F) + geom_smooth(data=Ttests.2, mapping = aes(y=sY , x=rho ),colour="black", se = F)+
  labs(x=expression(rho), y="Binary S-value")
h + geom_smooth(mapping = aes(y=sZ , x=rho), se = F) + geom_smooth(data=Ttests.2, mapping = aes(y=sZ , x=rho ),colour="black", se = F)+
  labs(x=expression(rho), y="Gaussian S-value")
h <- ggplot(data=Ttests)
h+ geom_smooth(mapping = aes(y=sY , x=rho, linetype=factor(set)), size = 1.5, se = F, colour="black") +
  labs(x=expression(rho), y="Binary S-value")+labs(linetype = expression("("*eta[z]*","*eta[y]*")"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=16, face = "bold"), 
        legend.text=element_text(size=16), legend.key.size = unit(1, "cm"),
        axis.title.y = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text = element_text(size=16))
pdf("transmissionYs.pdf")
h+ geom_smooth(mapping = aes(y=sZ , x=rho, linetype=factor(set)), size =1.5, se = F, colour="black") +
  labs(x=expression(rho), y="Gaussian S-value")+labs(linetype = expression("("*eta[z]*","*eta[y]*")"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=16, face = "bold"), 
        legend.text=element_text(size=16), legend.key.size = unit(1, "cm"),
        axis.title.y = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text = element_text(size=16))
pdf("transmissionZs.pdf")
#######################
ymeansd <- apply(ymeand, c(2,3), mean)
zmeansd <- apply(zmeand, c(2,3), mean)
SYsd <- apply(SBd, c(2,3), mean)
SZsd <- apply(SGd, c(2,3), mean)
etests.1 <- data.frame(cbind(testeta2s,ymeansd[,1],zmeansd[,1],SYsd[,1],SZsd[,1]), set = "unscaled")
etests.2 <- data.frame(cbind(testeta2s,ymeansd[,2],zmeansd[,2],SYsd[,2],SZsd[,2]), set = "scaled")
colnames(etests.2)<-colnames(etests.1)<- c("etaz","ymean","zmean","sY","sZ", "set") ;
e.tests <- rbind(etests.1,etests.2)

plot(data=etests.1,SZsd[,2]~testeta2s, type="l")+lines(data=Ttests.2, sZ~rho, col="red")
h <- ggplot(data=etests.1)
h + geom_smooth(mapping = aes(y=ymean , x=etaz), se = F) + geom_smooth(data=etests.2, mapping = aes(y=ymean , x=etaz ),colour="black", se = F)
h + geom_smooth(mapping = aes(y=sY , x=etaz), se = F) + geom_smooth(data=etests.2, mapping = aes(y=sY , x=etaz ),colour="black", se = F)+
  labs(x=expression(eta_y), y="Binary S-value")
h + geom_smooth(mapping = aes(y=sZ , x=etaz), se = F) + geom_smooth(data=etests.2, mapping = aes(y=sZ , x=etaz ),colour="black", se = F)+
  labs(x=expression(eta_y), y="Gaussian S-value")+ylim(0,1)
h <- ggplot(data=e.tests)
h+ geom_smooth(mapping = aes(y=sY , x=etaz, linetype=factor(set)), se = F, colour="black") +
  labs(x=expression(etaz), y="Binary S-value")+labs(linetype = expression("("*eta[z]*","*eta[y]*")")) + geom_abline(intercept=0,slope=1)

h+ geom_smooth(mapping = aes(y=ymean , x=etaz, linetype=factor(set)), se = F, colour="black", size=2) + ylim(0,0.62)+
  labs(x=expression(eta[z]), y="Binary Mean")+
  scale_linetype_discrete(name=expression(rho*" Transformation"), labels = c("unscaled", "scaled"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=16, face = "bold"), 
        legend.text=element_text(size=16), legend.key.size = unit(1, "cm"),
        axis.title.y = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text = element_text(size=16))

h+ geom_smooth(mapping = aes(y=sY , x=etaz, linetype=factor(set)), se = F, colour="black", size=2)+   geom_hline(yintercept = 1, linetype=11) +
  labs(x=expression(eta[z]), y="Binary S-Value")+
  scale_linetype_discrete(name=expression(rho*" Transformation"), labels = c("unscaled", "scaled", "actual"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=16, face = "bold"), 
        legend.text=element_text(size=16), legend.key.size = unit(1, "cm"),
        axis.title.y = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text = element_text(size=16))

################################################
Mist <- c(seq(-3,-0, by=.001),seq(0.001, 3, by =.01))
gams<- sapply(Mist, FUN=ubest)
bub <- ggplot(data=data.frame(Mist,gams))
bub + geom_line(mapping = aes(x=Mist, y = gams), colour = "black") +
  labs(x=expression(delta), y=expression("Upper Bound for  "*eta[y]))+theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=10), legend.key.size = unit(1, "cm"),
        axis.title.y = element_text(size=26),
        axis.title.x = element_text(size=32),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))
