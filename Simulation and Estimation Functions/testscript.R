nbss <- rectanggrid.4nbrs(k1=30,k2=30); k1=30;k2=30

sig2 = 1
bmeans <- c(0.5); gmeans <- c(1); spatpars <- c(1,.2,1)/c(1,1,1); 
bcovt <- NULL; gcovt<-NULL

ys <- rbinom(prob = exp(bmeans[1])/(1+exp(bmeans[1])), size = 1, n = k1*k2) 
zs <- rnorm(k1*k2, mean = gmeans, sd = sqrt(sig2))
tries <-spatial.genfieldbivar.covs(ys,zs,nbss,M=300,
                                   spatpars,bmeans, gmeans, sigmasq=sig2, bcovs=bcovt, gcovs=gcovt)
newtz <- tries[,2];newty <- tries[,1]
mean(newtz);sd(newtz);mean(newty); cor(newtz,newty)
ipcs <- initialparams.cov(newty,newtz,nbss,bcovt,gcovt)
ipcs;
SB <- unname(Dbinary(newty,nbss, k1=k1,k2=k2));SB
SG <- unname(Dgaussian(newtz,nbss,nbins=k2,sigma2=var(newtz)));SG
zc <-unname(newtz-mean(newtz)); yc <-unname(newty-mean(newty)); xmod <- summary(lm(zc~yc))
SX <- c(unname(xmod$coef[2])/unname(xmod$sigma), unname(xmod$sigma));SX[1]/SX[2];SX
glm(newty~bcovt+newtz, family="binomial")$coef; lm(newtz~gcovt+newty)$coef


rhoconst <- sqrt((1-spatpars[2])*(ubest(bmeans[1])-spatpars[3]));
spatpars[1]<-spatpars[1]*sqrt((1-spatpars[2])*(ubest(bmeans[1])-spatpars[3]))/rhoconst

#full w/ covariates
pps <- optim(par=ipcs,f=logplikbivar.covs, k1=k1, k2=k2,ys=newty,
             zs=newtz, bcovs=bcovt,gcovs=gcovt, nbrs=nbss,hessian=T,control=list(maxit=1000))$par; round(pps,5)

#full w/ no covariates
ppnc <- round(optim(par=initialparams.cov(newty,newtz,nbss),f=logplikbivar.covs, k1=k1, k2=k2,ys=newty,
                    zs=newtz, nbrs=nbss,hessian=T,control=list(maxit=500))$par,5);ppnc

#sqrt((1-ppnc[6]*ppnc[2])*(ubest(ppnc[4])-ppnc[3])) /sqrt(ppnc[6])
bcovt <- cbind(bcovt,newtz-mean(newtz)); gcovt <- cbind(gcovt,newty-mean(newty))
ipcs <- initialparams.cov(newty,newtz,nbss,bcovt,gcovt); ipcs
#binary univariate spatial
optim(par=c(ipcs[3:6]),f=logplikuni.b, k1=k1, k2=k2, bcovs=bcovt,
      ys=newty, nbrs=nbss,hessian=T,control=list(maxit=500))$par
#binary univariate nonspatial
optim(par=c(ipcs[4:5]),f=logplik.b.ns, k1=k1, k2=k2, bcovs=bcovt,
      ys=newty, nbrs=nbss,control=list(maxit=500), hessian=T)$par
glm(newty~bcovt, family="binomial")

#gaussian univariate spatial
optim(par=c(ipcs[2],ipcs[6:8]),f=logplikuni.g, k1=k1, k2=k2, gcovs=gcovt,
      zs=newtz, nbrs=nbss,hessian=T,control=list(maxit=500))$par
#gaussian univariate nonspatial
optim(par=ipcs[6:8],f=logplikuni.g.ns, k1=k1, k2=k2, gcovs=gcovt,
      zs=newtz, nbrs=nbss,hessian=T,control=list(maxit=500))$par
#bivariate nonspatial
optim(par=c(ipcs[1], ipcs[4:6]),f=logplikbivar.ns, k1=k1, k2=k2, gcovs=gcovt, bcovs=bcovt,
      zs=newtz, ys=newty, nbrs=nbss,hessian=T,control=list(maxit=500))$par
lm(newtz~gcovt[,1]+gcovt[,2])
glm(newty~bcovt[,1]+bcovt[,2],family="binomial")



image(matrix(newty, nrow = sqrt(length(nbss[,1]))),col = gray.colors(2))
image(matrix(newtz, nrow = sqrt(length(nbss[,1]))), col = inferno(200))
image(matrix(bcovt[,1], nrow = sqrt(length(nbss[,1]))), col = inferno(20))

covt<- cbind(seq(1,k1, by = 1)*seq(1,k1*k2, by =1)*0.01, covg[,1]) ; mean(covt[,1]); mean(covt[,2])
s <-matrix(1,nrow=k1,ncol=k2)
xmult <- 0.5*seq(1,k1,by=1)
ymult <- 0.0005*seq(1,k2,by=1)^2
s1 <- t(s*xmult)*ymult
s2 <- rev(t(s)*ymult)
image(matrix(rev(s1), nrow = sqrt(length(nbss[,1]))), col = inferno(10))
summary(c(s1));summary(c(s2))
covg<- spatial.genfieldbivar.covs(ys,zs,nbss,M=300,
                                  c(0,0.5,5),c(-0.5,1), c(1,-0.5), sigmasq=sig2, bcovs=c(s1), gcovs=c(s1)); colMeans(covg);

gcov2 <- rgamma(n=k1*k2, shape=3, scale=4)
bcov2 <- rexp(n=k1*k2, rate=.2)
image(matrix(covg[,2], nrow=30,ncol=30), col = inferno(50))
hist(bcov2)
bcovt <- cbind(covg[,2],bcov2)
bcovt <-c(covg[,2]); gcovt <- c(gcov2)
gcovt <- cbind(covg[,2],gcov2)


initialparams.cov(newty,newtz,nbss,bcovt[,1],gcovt[,1])
