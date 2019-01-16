#Functions to estimate the constant mean S-value from Kaiser and Caragea (2009)
library(tidyverse)
#################################
#binary field first
Dbinary <- function(ys, nbrs, k1, k2){
  kk1=k1
  kk2=k2
  #only use interior locations
    locs=1:(kk1*kk2)
    locmat<-matrix(locs,kk1,kk2,byrow=T)
    #remove border locations
    locmat<-locmat[-c(1,kk1),-c(1,kk2)]
    intlocs<-as.vector(locmat)
    intlocs<-intlocs[order(intlocs)]
  #get kappa estimator first
  ktilde <- mean(ys)
  #set of possible values for (1/m)*sum(nbrs)
  hls <- c(1,(1/4),(2/4),(3/4), 0)
  
  yhls <- c()
  cnt <- 0
  repeat{
    cnt <- cnt + 1
    tloc<-intlocs[cnt];tnbs<-nbrs[tloc,];#tnbs<-tnbs[tnbs>0]
    yhls[cnt] <- (1/4)*sum(ys[tnbs])
    if(cnt == length(intlocs)){break}
  }
  #Get D(hl,k) for each y
  Dhls <- data.frame(loc=intlocs,y=ys[intlocs],hl=yhls, Dhlk = yhls-ktilde)

  #Now get Cs for each
  #define sets
  Hl <- seq(0,length(unique(Dhls$hl))-1)
  #bin data
  sety <- 4*Dhls$hl
  Dhls <- data.frame(Dhls, sety = as.factor(sety))
  #Get Cs for each H
  Cs<- Dhls %>% group_by(g=sety) %>% summarise(HL = length(sety), sumY = sum(y)) %>% mutate(Cl = (1/HL)*sumY)
  Cobs <- c()
  #Apply Cs to each observation
  for(j in Hl){
    yjs <- which(sety==j)
    Cobs[yjs]<- Cs$Cl[j+1]
    j<-j+1
  }
  #get r values
  Dhls <- data.frame(Dhls, Cl = Cobs) %>% mutate(rClk= log(Cl/(1-Cl))-log(ktilde/(1-ktilde))) %>%
    filter(!is.na(rClk) & !is.infinite(rClk))
  #get S statistic
  Sreg <- lm(data=Dhls, formula=rClk~0+Dhlk)
  Sstat <- coef(Sreg)
  Svar <- vcov(Sreg)
  return(c(Sstat,Svar))
}

#####################################################################
#gaussian
Dgaussian <- function(zs, nbrs, nbins, sigma2=0){
  #only use interior locations
  locs=1:(k1*k2)
  locmat<-matrix(locs,k1,k2,byrow=T)
  #remove border locations
  locmat<-locmat[-c(1,k1),-c(1,k2)]
  intlocs<-as.vector(locmat)
  intlocs<-intlocs[order(intlocs)]
  #get kappa estimator first
  sig2z<-var(zs)
  if(sigma2>0){sig2z <- sigma2}
  ktilde <- mean(zs)
  #set of possible values for (1/m)*sum(nbrs)
  cbs <- quantile(zs, prob= c(seq(0,1, length=nbins)))
  zhls <- c()
  cnt <- 0
  repeat{
   cnt <- cnt + 1
  tloc<-intlocs[cnt];
  loc<-intlocs[cnt];tnbs<-nbrs[tloc,];
  zhls[cnt] <- (1/4)*sum(zs[tnbs])
  if(cnt == length(intlocs)){break}
  }
  cbs <- quantile(zhls, prob= c(seq(0,1, length=nbins)))
  zhls <- .bincode(zhls, breaks=cbs, include.lowest=T)
  #return(zhls)
  #Get D(hl,k) for each z
  hl <- c()
  j=1
  for(j in unique(zhls)){
    zjs <- which(zhls==j)
    hl[zjs]<- cbs[j]
    j<-j+1
  }
  Dhlk = hl-ktilde
  #return(c(length(zjs),length(hl)))
  Dhls <- data.frame(loc=intlocs,z=zs[intlocs],hl=hl, Dhlk)
  #Now get Cs for each
  #define sets
  #Hl <- seq(0,length(unique(Dhls$hl))-1)
  #bin data
  #sety <- 4*Dhls$hl
  Dhls <- data.frame(Dhls, setzs = as.factor(zhls))
  #Get Cs for each H
  Cs<- Dhls %>% group_by(g=setzs) %>% summarise(HL = length(g), sumz = sum(z)) %>% mutate(Cl = (1/HL)*sumz)
  Cobs <- c()
  #Apply Cs to each observation
  j=1
  for(j in unique(zhls)){
    zjs <- which(zhls==j)
    Cobs[zjs]<- Cs$Cl[j]
    j<-j+1
  }
  #get r values
  Dhlz <- data.frame(Dhls, Cl = Cobs) %>% mutate(rClk= Cl/sig2z - ktilde/sig2z)%>%
    filter(!is.na(rClk) & !is.infinite(rClk))
  #get S statistic
  Sreg <- lm(data=Dhlz, formula=rClk~0+Dhlk)
  Sstat <- coef(Sreg)*sig2z
  Svar <- vcov(Sreg)*sig2z
  return(c(Sstat, Svar))
}


##########################################
##########################################
##########################################
