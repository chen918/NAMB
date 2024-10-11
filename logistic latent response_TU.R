set.seed(99999)
if (!require(mvtnorm)) {
  install.packages("mvtnorm")
  library(mvtnorm)
}
if (!require(invgamma)) {
  install.packages("invgamma")
  library(invgamma)
}
if (!require(sna)) {
  install.packages("sna")
  library(sna)
}
if (!require(coda)) {
  install.packages("coda")
  library(coda)
}
if (!require(truncnorm)) {
  install.packages("truncnorm")
  library(truncnorm)
}
# True value of Rho
# rho=0.2/-0.2
rho=0

itl=1000
medianbeta_m = matrix(0,nrow=4,ncol=itl)
medianRho_v = c()
coverage = 0
ubv = c()
lbv = c()

k=1

while(k<=itl){
  
  ################
  #Generating the data
  node=200 
  ddensity=0.4
  X=cbind(1,matrix(rnorm(node*3),ncol=3))
  BBeta=matrix(c(0.5,2,1,0.5),nrow =4)
  Tao=matrix(rlogis(node, location = 0, scale = 1),ncol = 1)
  W=rgraph(node, m=1, tprob=ddensity, mode="graph", diag = F, replace=FALSE,
           tielist=NULL, return.as.edgelist= F )
  
  W[lower.tri(W)][which(W[lower.tri(W)]!=0)]=rgamma(n = length(which(W[lower.tri(W)]!=0)), shape = 0.025, scale = 2000)
  
  W[upper.tri(W)]<-t(W)[upper.tri(W)]
  
  W=t(apply(W, 1, function(x) (x/sum(x))))
  
  W[which(is.na(W))]=1/(node-1)
  
  diag(W)<-0
  
  DDelta= solve(diag(x=1,nrow=node)-rho*W)%*%(X%*%BBeta+Tao)
  
  Y=as.numeric(DDelta>0)
  
  # ordered eigenvalues of W
  Lamda = c(sort(eigen(W)$values, decreasing = T))
  
  Lamda = Re(Lamda)
  
  
  ######
  #MCMC
  num=20000
  
  Rho = rep(0, num)
  betamt = matrix(0, 4, num)
  detamt = matrix(0, node, num)
  
  Rho[1]=rho
  betamt[,1]=BBeta
  detamt[,1]=DDelta
  ####
  I = t(rep(1,node))
  V = solve(t(X)%*%t(solve(diag(1,node)-rho*W))%*%diag(c(exp(solve(diag(1,node)-rho*W)%*%X%*%BBeta)/(1+exp(solve(diag(1,node)-rho*W)%*%X%*%BBeta))^2))%*%solve(diag(1,node)-rho*W)%*%X)
  count = 0
  countBt = 0
  countDt = 0
  i = 2
  while(i<=num){
    Rho[i] = Rho[i-1]+rnorm(1, mean = 0, sd = 0.15)
    
    if (Rho[i]>1/max(Lamda) | Rho[i]< 1/min(Lamda)){
      next
    }
    
    I1 = diag(1,node)-Rho[i]*W
    I2 = diag(1,node)-Rho[i-1]*W
    XB1 = X%*%betamt[,i-1]
    
    cal1 <- function(g){
      y=exp(I2[g,]%*%detamt[,i-1]-I1[g,]%*%detamt[,i-1])*((1+exp(-I2[g,]%*%detamt[,i-1]+XB1[g]))/(1+exp(-I1[g,]%*%detamt[,i-1]+XB1[g])))^2
      return(y)
    }
    
    proportion = (((1/max(Lamda)-Rho[i-1])*(Rho[i-1]-1/min(Lamda)))/((1/max(Lamda)-Rho[i])*(Rho[i]-1/min(Lamda))))*abs(det(diag(1,node)-Rho[i]*W)/det(diag(1,node)-Rho[i-1]*W))*prod(unlist(lapply(seq(1,node),cal1)))* 
      dtruncnorm(Rho[i-1], a=1/min(Lamda), b=1/max(Lamda), mean=Rho[i],sd=0.15)/dtruncnorm(Rho[i], a=1/min(Lamda), b=1/max(Lamda), mean=Rho[i-1],sd=0.15)
    
    proportion = min(1, Re(proportion))
    
    if ((runif(1) > proportion)){
      Rho[i]=Rho[i-1]
      count = count+1
    }
    
    
    betamt[,i] = betamt[,i-1]+rmvnorm(1, mean = rep(0,4), sigma=0.6*V)
    
    I3 = diag(1,node)-Rho[i]*W
    XB2 = X%*%betamt[,i]
    XB3 = X%*%betamt[,i-1]
    
    cal2 <- function(g){
      y=exp(XB2[g]-XB3[g])*((1+exp(-I3[g,]%*%detamt[,i-1]+XB3[g]))/(1+exp(-I3[g,]%*%detamt[,i-1]+XB2[g])))^2
      return(y)
    }
    
    proportionbeta = prod(unlist(lapply(seq(1,node),cal2)))
    
    proportionbeta = min(1, Re(proportionbeta))
    
    if (runif(1) > proportionbeta){
      betamt[,i]=betamt[,i-1]
      countBt = countBt+1
    }
    
    
    for(j in 1:node){
      
      detamt[j,i] = detamt[j,i-1]+rnorm(1, mean = 0, sd = 0.3)
      
      if(Y[j]==1){
        if (detamt[j,i]<0){
          detamt[j,i]=detamt[j,i-1]
        }
      }else{
        if (detamt[j,i]>0){
          detamt[j,i]=detamt[j,i-1]
        }
      }
    }
    
    XB4 = X%*%betamt[,i]
    
    cal3 <- function(g){
      y=exp(I3[g,]%*%detamt[,i-1]-I3[g,]%*%detamt[,i])*((1+exp(-I3[g,]%*%detamt[,i-1]+XB4[g]))/(1+exp(-I3[g,]%*%detamt[,i]+XB4[g])))^2
      return(y)
    }
    
    proportiondeta = prod(unlist(lapply(seq(1,node),cal3)))
    
    proportiondeta = min(1, Re(proportiondeta))
    
    if (runif(1) > proportiondeta){
      detamt[,i]=detamt[,i-1]
      countDt = countDt+1
    }
    
    i = i+1
    
  }
  
  medianbeta = apply(betamt[,-(1:2000)],MARGIN=1, FUN=median)
  
  medianbeta_m[,k] = medianbeta
  medianRho_v[k] = median(Rho[-(1:2000)])
  
  Rate<-(num-count-1)/(num-1)
  RateB<-(num-countBt-1)/(num-1)
  RateD<-(num-countDt-1)/(num-1)
  
  print(Rate)
  print(RateB)
  print(RateD)
  
  cat("k = ", k,"\n")
  
  lb=quantile(Rho[-(1:2000)],0.025)
  ub=quantile(Rho[-(1:2000)],0.975)
  lbv = c(lbv,lb)
  ubv = c(ubv,ub)
  if(rho<=ub & rho>=lb){
    coverage=coverage+1
  }
  
  k = k + 1
  
}

apply(medianbeta_m,MARGIN=1, FUN=mean)

mean(medianRho_v)

# Bias
mean(medianRho_v)-rho

# MSE
sum((medianRho_v-rho)^2)/itl

# Coverage rate
coverage/itl

mean(lbv)

mean(ubv)

mean(ubv)-mean(lbv)

