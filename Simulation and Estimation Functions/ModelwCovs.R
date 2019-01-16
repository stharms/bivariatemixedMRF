#4 nearest neighbors function
rectanggrid.4nbrs<-function(k1,k2){
  us<-1:k1
  vs<-1:k2
  nbmat<-matrix(0,k1*k2,4)
  ucnt<-0
  repeat{
    ucnt<-ucnt+1
    tu<-us[ucnt]
    vcnt<-0
    repeat{
      vcnt<-vcnt+1
      tv<-vs[vcnt]
      nbrs<-NULL
      tsite<-(tu-1)*k2+tv
      for(ui in (tu-1):(tu+1)){
        if((ui>0) && (ui<=k1) && (tv>0) && (tv<=k2))
          nbrs<-c(nbrs,(ui-1)*k2 +tv) 
        else nbrs<-c(nbrs,0)}
      for(vi in (tv-1):(tv+1)){
        if((tu>0) && (tu<=k1) && (vi>0) && (vi<=k2))
          nbrs<-c(nbrs,(tu-1)*k2 +vi) 
        else nbrs<-c(nbrs,0)}
      nbmat[tsite,]<-nbrs[nbrs!=tsite]
      if(vcnt==k2) break
    }
    if(ucnt==k1) break
  }
  return(nbmat)
}

#Simulation function via node-by-node Gibbs sampling scheme
simbivar.covs<-function(spatialpars,curys, curzs, nbrs, bcovs=NULL, bcovpars, gcovs=NULL, gcovpars, sigmasq){

  n=length(curys)
  #initialize covariates. add a column of 1s for the intercept.
  bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs);
  #parameters used for simulation
  mu <- gcov%*%gcovpars[1:ncol(gcov)]; delta = bcov%*%bcovpars[1:ncol(bcov)];
  bet <- sigmasq
  rho <- spatialpars[1]
  eta1 <-  spatialpars[2]
  eta2 <-  spatialpars[3]
  margy <- exp(delta)/(1+exp(delta))
  
  cnt<-0
  #iterate through the binary nodes
  repeat{
    cnt<-cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    theta <- (delta[tloc] + (eta2/(length(curys[tnbs])))*sum(curys[tnbs]-margy[tnbs]) + (rho)/bet*(curzs[cnt]-mu[tloc]))
    py <- exp(theta)/(1+exp(theta))
    ny <- ifelse(runif(1)<py, 1, 0)		
    curys[cnt]<-ny
    if(cnt==n) break
  }
  cnt<-0
  #iterate through the gaussian nodes
  repeat{
    cnt<-cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    nz <- rnorm(n=1, mean = ((mu[tloc]) + (eta1/(length(tnbs)))*sum(curzs[tnbs]-mu[tnbs]) + rho*(curys[tloc]-(margy[tloc]))), sd = sqrt(bet))
    curzs[cnt]<- nz
    if(cnt==n) break
  }
  #output the data
  data.frame(curys, curzs)
}

###################################################################################
######################################################################################
#generate lattice data by iterating the simulation function M times
spatial.genfieldbivar.covs<-function(startys,startzs,nbrs,M,spatialpars, bmeanpars, gmeanpars, sigmasq, bcovs=NULL,gcovs=NULL){
  curys<-startys
  curzs<-startzs
  cnt<-0
  repeat{
    cnt<-cnt+1
    newobs<-simbivar.covs(spatialpars,curys,curzs,nbrs, bcovs,bmeanpars,gcovs, gmeanpars, sigmasq)
    curys<-newobs[,1]
    curzs<-newobs[,2]
    if(cnt==M) break
  }
  data.frame(curys, curzs)
}

###################################################################################
##################################################################################
####################################
####################################
####################################
#####################################
##Log-pseudolikelihood functions#####
#####################################
#compute conditional canonicals for gaussian variables
alpha.covs <- function(fixedpars, ys, zs, nbrs,bcovs,gcovs){
  n<-length(zs)
  bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs)
  rho <- fixedpars[1]
  eta1 <- fixedpars[2]
  eta2 <- fixedpars[3]
  zcnt <- 4+ncol(bcov)
  bmeanpar <- fixedpars[4:(zcnt-1)]; 
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
    alphas[cnt] <- ((mu[tloc])+(eta1/length(zs[tnbs]))*sum(zs[tnbs]-mu[tnbs])+ rho*(ys[tloc]-margy[tloc]))
    if(cnt==n){break}
  }
  return(alphas)
}
#compute canonical parameters for binary variables
theta.covs <- function(fixedpars, ys, zs, nbrs,bcovs,gcovs){
  n=length(ys)
  bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs)
  rho <- fixedpars[1]
  eta1 <- fixedpars[2]
  eta2 <- fixedpars[3]
  zcnt <- 4+ncol(bcov)
  bmeanpar <- fixedpars[4:(zcnt-1)]; 
  gmeanpar <- fixedpars[zcnt:(zcnt+ncol(gcov)-1)]
  bet <- fixedpars[length(fixedpars)]
  
  mu <- gcov%*%gmeanpar[1:ncol(gcov)]; delta = bcov%*%bmeanpar[1:ncol(bcov)];
  margy <- exp(delta)/(1+exp(delta))
  
  cnt<-0
  thetas <- NULL
  repeat{
    cnt <- cnt+1
    tloc<-cnt;tnbs<-nbrs[tloc,];tnbs<-tnbs[tnbs>0]
    thetas[cnt] <- delta[tloc] + (eta2/length(ys[tnbs]))*sum(ys[tnbs]-margy[tnbs]) + rho*(zs[tloc]-mu[tloc])/bet
    if(cnt==n){break}
  }
  return(thetas)
}
##############################################
###############################
#function to estimate the log pseudo-likelihood
logplikbivar.covs<- function(k1,k2, allpars, ys, zs, nbrs, bcovs=NULL, gcovs=NULL){
  y <- ys
  z <- zs
  n=length(ys)
  #return(length(bcovs))
  bcov <- cbind(rep(1,times=n),bcovs); gcov <- cbind(rep(1,times=n),gcovs)
  rho <- allpars[1]
  eta1 <- allpars[2]
  eta2 <- allpars[3]
  zcnt <- 4+ncol(bcov)
  bmeanpar <- allpars[4:(zcnt-1)]; 
  gmeanpar <- allpars[zcnt:(zcnt+ncol(gcov)-1)]
  bet <- allpars[length(allpars)]
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  thetys <- theta.covs(allpars, y, z, nbrs,bcovs,gcovs)
  bthetys <- log(1+exp(thetys))
  Alphas <- alpha.covs(allpars, y ,z, nbrs,bcovs,gcovs)
  logpliklihood <- (-1/2)*length(intlocs)*log((2*bet*pi)) - sum((1/(2*bet))*(zs[intlocs]-Alphas[intlocs])^2) + 
    sum(ys[intlocs]*thetys[intlocs]-bthetys[intlocs])
  
  return(-1*logpliklihood)
}


#A function to get some initial parameters for optimizing the log pseudo-likelihood
initialparams.cov <- function(newys, newzs, nbrs, bcovs=NULL, gcovs=NULL){
  #return(ncol(gcovs))
  
  nbsr <- apply(nbrs, 1, FUN = prod)
  nbrr <- nbrs
  nbrr[which(nbrr==0)]<-1
  neighbs <- cbind(newzs, newzs[nbrr[,1]],newzs[nbrr[,3]],newzs[nbrr[,4]],newzs[nbrr[,2]],gcovs)[-which(nbsr==0),]
  fitzs <- lm(neighbs[,1]~neighbs[,2]+neighbs[,3]+neighbs[,4]+neighbs[,5])
  zbs <- lm(newzs~1); ybs <- glm(newys~1, family="binomial")
  if(!is.null(gcovs)){zbs <- lm(newzs~gcovs);}
  if(!is.null(bcovs)){ybs <- glm(newys~bcovs, family="binomial")}
  beti<-anova(fitzs)$'Mean Sq'[length(anova(fitzs)$'Mean Sq')]
  if(is.null(gcovs)){beti <- anova(fitzs)$'Mean Sq'[length(anova(fitzs)$'Mean Sq')]}
  
  eta1i <- sum(fitzs$coefficients[2:5])
  neighbsy <- cbind(newys, newys[nbrr[,1]],newys[nbrr[,3]],newys[nbrr[,4]],newys[nbrr[,2]])[-which(nbsr==0),]
  fitys <- glm(neighbsy[,1]~neighbsy[,2]+neighbsy[,3]+neighbsy[,4]+neighbsy[,5], family = "binomial")
  eta2i<-sum(fitys$coefficients[2:5])
  gcov <- gcovs; if(is.null(gcovs)){gcov<-rep(1,times=length(newys))}
  rho <- unname(lm(c(newzs-mean(newzs))~c(newys-mean(newys))+gcov)$coefficients[2])
  return(c(rh=rho, e1=eta1i, e2=eta2i,  bmp = coef(ybs), gmp = coef(zbs), sigsq=beti ))
}


#logplikbivar.covs(k1,k2, allpars=c(spatpars,bmeans,gmeans,1), newty,newtz,nbss,bcovs=bcovt[,2],gcovs=gcovt[,2])

