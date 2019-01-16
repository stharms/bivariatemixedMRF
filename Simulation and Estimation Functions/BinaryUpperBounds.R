
Qyest<- function(AI, kap){
  Qy <- (exp(AI)-kap*(1+exp(AI)))/((1+exp(AI))*(AI-log(kap)+log(1-kap)))
  #Fi <- rho*Qy
  return(-Qy)
  #(-1/2)*length(intlocs)*log((2*bet*pi)) 
}

#This is the Q function
Qys <- function(AI, Mi){
  Qy <- (((exp(AI))/(1+exp(AI))) - (exp(Mi))/(1+exp(Mi)))/(AI-Mi)
  return(-Qy)
}
#Optimize it in 1-D
ubest<-function(delta){
  (optim(par=0.5,f=Qys, Mi = delta,
         hessian=T,control=list(maxit=500), method = "Brent", lower = -10, upper=10)$value*(-1))^(-1)
}

