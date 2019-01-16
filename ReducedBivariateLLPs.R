alpha.biv.ns <- function(fixedpars, ys, zs, nbrs, bcovs=NULL,gcovs=NULL){
    n=length(ys)
    bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs)
    rho <- fixedpars[1]
    zcnt <- 2+ncol(bcov)
    bmeanpar <- fixedpars[2:(zcnt-1)]; 
    gmeanpar <- fixedpars[zcnt:(zcnt+ncol(gcov)-1)]
    bet <- fixedpars[length(fixedpars)]
    
    mu <- gcov%*%gmeanpar[1:ncol(gcov)]; delta = bcov%*%bmeanpar[1:ncol(bcov)];
    margy <- exp(delta)/(1+exp(delta))
  
  cnt<-0
  alphas <- NULL
  As <- NULL;Bs <- NULL;Cs <- NULL;
  repeat{
    cnt <- cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    alphas[cnt] <- ((mu[tloc])+rho*(ys[tloc]-margy[tloc]))
    if(cnt==n){break}
  }
  return(alphas)
}
#compute canonicals for binary variables (non-spatial)
theta.biv.ns <- function(fixedpars, ys, zs, nbrs, bcovs=NULL,gcovs=NULL){
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs)
  rho <- fixedpars[1]
  zcnt <- 2+ncol(bcov)
  bmeanpar <- fixedpars[2:(zcnt-1)]; 
  gmeanpar <- fixedpars[zcnt:(zcnt+ncol(gcov)-1)]
  bet <- fixedpars[length(fixedpars)]
  
  mu <- gcov%*%gmeanpar[1:ncol(gcov)]; delta = bcov%*%bmeanpar[1:ncol(bcov)];
  margy <- exp(delta)/(1+exp(delta))
  
  cnt<-0
  thetas <- NULL
  repeat{
    cnt <- cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    thetas[cnt] <- delta[tloc] + rho*(zs[tloc]-mu[tloc])/bet
    if(cnt==n){break}
  }
  return(thetas)
}
##############################################

###############################

logplikbivar.ns<- function(k1,k2, allpars, ys, zs, nbrs, bcovs=NULL,gcovs=NULL){
  y <- ys
  z <- zs
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs)
  rho <- allpars[1]
  zcnt <- 2+ncol(bcov)
  bmeanpar <- allpars[2:(zcnt-1)]; 
  gmeanpar <- allpars[zcnt:(zcnt+ncol(gcov)-1)]
  bet <- allpars[length(allpars)]
  
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  
  thetys <- theta.biv.ns(allpars, y, z, nbrs,bcovs,gcovs)
  bthetys <- log(1+exp(thetys))
  Alphas <- alpha.biv.ns(allpars, y ,z, nbrs,bcovs,gcovs)

  logpliklihood <- (-1/2)*length(intlocs)*log((2*bet*pi)) - sum((1/(2*bet))*(zs[intlocs]-Alphas[intlocs])^2) + 
    sum(ys[intlocs]*thetys[intlocs]-bthetys[intlocs])

  return(-1*logpliklihood)
}

################################################
theta.uni.s <- function(fixedpars, ys, nbrs,bcovs=NULL){
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs)
  eta2 <- fixedpars[1]
  bmeanpar <- fixedpars[2:length(fixedpars)]; 
  delta = bcov%*%bmeanpar[1:ncol(bcov)];
  margy <- exp(delta)/(1+exp(delta))
  
  cnt<-0
  thetas <- NULL
  repeat{
    cnt <- cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    thetas[cnt] <- delta[tloc] + (eta2/length(ys[tnbs]))*sum(ys[tnbs]-margy[tnbs])
    if(cnt==n){break}
  }
  return(thetas)
}

logplikuni.b<- function(k1,k2, allpars, ys, nbrs, bcovs=NULL){
  y <- ys
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs)
  eta2 <- allpars[1]
  bmeanpar <- allpars[2:length(allpars)]; 
  delta = bcov%*%bmeanpar[1:ncol(bcov)];
  margy <- exp(delta)/(1+exp(delta))
  
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  
  thetys <- theta.uni.s(allpars, y, nbrs, bcovs)
  bthetys <- log(1+exp(thetys))
  
  logpliklihood <-  sum(ys[intlocs]*thetys[intlocs]-bthetys[intlocs])
  return(-1*logpliklihood)
}

#optim(par=c(0.5,2),f=logplikuni.b, k1=k1, k2=k2,
#ys=newys, nbrs=nbss,hessian=T,control=list(maxit=500))$par

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
theta.uni.ns <-  function(fixedpars, ys, nbrs,bcovs=NULL){
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs)
  bmeanpar <- fixedpars; 
  delta = bcov%*%bmeanpar[1:ncol(bcov)];
  margy <- exp(delta)/(1+exp(delta))
  
  cnt<-0
  thetas <- NULL
  repeat{
    cnt <- cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    thetas[cnt] <- delta[tloc]
    if(cnt==n){break}
  }
  return(thetas)
}

logplik.b.ns<- function(k1,k2, allpars, ys, nbrs, bcovs=NULL){
  y <- ys
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs)
  bmeanpar <- allpars; 
  delta = bcov%*%bmeanpar[1:ncol(bcov)];
  margy <- exp(delta)/(1+exp(delta))
  
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  
  thetys <- theta.uni.ns(allpars, ys, nbrs, bcovs)
  
  bthetys <- log(1+exp(thetys))
  
  logpliklihood <-  sum(ys[intlocs]*thetys[intlocs]-bthetys[intlocs])
  return(-1*logpliklihood)
}

#optim(par=c(0.5),f=logplik.b.ns, k1=k1, k2=k2,
#ys=newys, nbrs=nbss,control=list(maxit=500), method="Brent", lower=-10, upper=10)$par

#alphas for independent gaussian spatial
alphaug <- function(fixedpars, zs, nbrs,gcovs=NULL){
  n<-length(zs)
  gcov <- cbind(rep(1,times=n),gcovs)
  eta1 <- fixedpars[1]
  gmeanpar <- fixedpars[2:(length(fixedpars)-1)]
  bet <- fixedpars[length(fixedpars)]
  mu <- gcov%*%gmeanpar[1:ncol(gcov)];
  cnt<-0
  alphas <- NULL
  As <- NULL;Bs <- NULL;Cs <- NULL;
  repeat{
    cnt <- cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    alphas[cnt] <- ((mu[tloc])+(eta1/length(zs[tnbs]))*sum(zs[tnbs]-mu[tnbs]))
    if(cnt==n){break}
  }
  return(alphas)
}

#log pseudolikelihood for independent gaussian spatial
logplikuni.g<- function(k1,k2, allpars, zs, nbrs, gcovs=NULL){
  z <- zs
  n<-length(zs)
  gcov <- cbind(rep(1,times=n),gcovs)
  eta1 <- allpars[1]
  gmeanpar <- allpars[2:(length(allpars)-1)]
  bet <- allpars[length(allpars)]
  
  #interior locations
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  
  Alphas <- alphaug(allpars, z, nbrs, gcovs)
  betazs <- bet*(Alphas^2)/2
  czs <- (-(zs)^2)/(2*bet)
  
  logpliklihood <- (-1/2)*length(intlocs)*log((2*bet*pi)) - sum((1/(2*bet))*(zs[intlocs]-Alphas[intlocs])^2)
  return(-1*logpliklihood)
}


#alphas for independent gaussian non-spatial
alphaug.ns <- function(fixedpars, zs, nbrs, gcovs=NULL){
  n<-length(zs)
  gcov <- cbind(rep(1,times=n),gcovs)
  gmeanpar <- fixedpars[1:(length(fixedpars)-1)]
  bet <- fixedpars[length(fixedpars)]
  mu <- gcov%*%gmeanpar[1:ncol(gcov)];
  
  cnt<-0
  alphas <- NULL
  As <- NULL;Bs <- NULL;Cs <- NULL;
  repeat{
    cnt <- cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    alphas[cnt] <- mu[tloc]
    if(cnt==n){break}
  }
  return(alphas)
}

#log pseudolikelihood for independent gaussian non-spatial
logplikuni.g.ns<- function(k1,k2, allpars, zs, nbrs,gcovs=NULL){
  z <- zs
  n<-length(zs)
  gcov <- cbind(rep(1,times=n),gcovs)
  gmeanpar <- allpars[1:(length(allpars)-1)]
  bet <- allpars[length(allpars)]
  
  #interior locations
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  
  Alphas <- alphaug.ns(allpars,z,nbrs,gcovs)
  betazs <- bet*(Alphas^2)/2
  czs <- (-(zs)^2)/(2*bet)
  
  logpliklihood <- (-1/2)*length(intlocs)*log((2*bet*pi)) - sum((1/(2*bet))*(zs[intlocs]-Alphas[intlocs])^2)
  return(-1*logpliklihood)
}
#optim(par=c(100,100),f=logplikuni.g.ns, k1=k1, k2=k2,
#zs=newzs, nbrs=nbss,hessian=T,control=list(maxit=500))$par
