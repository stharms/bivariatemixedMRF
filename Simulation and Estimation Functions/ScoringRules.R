library(scoringRules)
########################
brierscore <- function(spatialpars,ys, zs, nbrs, bcovs=NULL, bcovpars, gcovs=NULL, gcovpars, sigmasq){
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs)
  mu <- gcov%*%gcovpars[1:ncol(gcov)]; delta = bcov%*%bcovpars[1:ncol(bcov)];
  bet <- sigmasq
  rho <- spatialpars[1]
  eta1 <-  spatialpars[2]
  eta2 <-  spatialpars[3]
  margy <- exp(delta)/(1+exp(delta))
  
  #use only interior locations
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  
  #Calculate Brier Score for Constant mean model = p(1-p)
  brierconst <- var(ys[intlocs])
  #Calculate Brier Score for Each Location
  n<-length(intlocs)
  bs<-c()
  proby <- mean(ys[intlocs])
  crpsc <- c()
  crpsm <- c()
  cnt<-0
  repeat{
    cnt<-cnt+1
    tloc<-intlocs[cnt];tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    
    crpsc[cnt] <- crps_binom(y=ys[tloc], size=1, prob=proby)
    
    theta <- (delta[tloc] + (eta2/4)*sum(ys[tnbs]-margy[tnbs]) + (rho)*(zs[tloc]-mu[tloc]))/bet
    py <- exp(theta)/(1+exp(theta))
    
    bs[cnt] <- (ys[tloc]-py)^2
    
    crpsm[cnt] <- crps_binom(y=ys[tloc], size=1, prob=py)
    if(cnt==n) break
  }
  briermod <- mean(bs)
  crpsmod <- mean(crpsm)
  crpscons <- mean(crpsc)
  
  #skill score
  skillb <- 100* (1-briermod/brierconst)
  skillc <- 100* (1- crpsmod/(crpscons))
  return(c(skillc,crpsmod, crpscons))
}
##################################################################################


GaussScore <- function(spatialpars,ys, zs, nbrs, bcovs=NULL, bcovpars, gcovs=NULL, gcovpars, sigmasq){
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs)
  mu <- gcov%*%gcovpars[1:ncol(gcov)]; delta = bcov%*%bcovpars[1:ncol(bcov)];
  bet <- sigmasq
  rho <- spatialpars[1]
  eta1 <-  spatialpars[2]
  eta2 <-  spatialpars[3]
  margy <- exp(delta)/(1+exp(delta))
  
  #use only interior locations
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  
  #Calculate CRPS Score for constant mean and model for Each Location
  n<-length(intlocs)
  mz <- mean(zs[intlocs]); sdz <- sd(zs[intlocs])
  crpsc<- c()
  crpsm<-c()
  cnt<-0
  
  repeat{
    cnt<-cnt+1
    tloc<-intlocs[cnt];tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    
    crpsc[cnt] <- crps_norm(y=zs[tloc], mean=mz, sd=sdz)
  
    pz <- ((mu[tloc]) + (eta1/4)*sum(zs[tnbs]-mu[tnbs]) + rho*(ys[tloc]-(margy[tloc])))
    crpsm[cnt] <- crps_norm(y=zs[tloc], mean=pz, sd=sqrt(bet))
    if(cnt==n) break
  }
  scoremod <- mean(crpsm)
  scoreconst <- mean(crpsc)
  #skill score
  skillz <- 100* ((1-scoremod/scoreconst))
  return(c(skillz, scoremod,scoreconst))
}
#######################################


