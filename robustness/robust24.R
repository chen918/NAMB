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
  rho=0
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
  
  Rho[1]=rho
  betamt[,1]=BBeta
  ####
  I = t(rep(1,node))
  V = solve(t(X)%*%t(solve(diag(1,node)-rho*W))%*%diag(c(exp(solve(diag(1,node)-rho*W)%*%X%*%BBeta)/(1+exp(solve(diag(1,node)-rho*W)%*%X%*%BBeta))^2))%*%solve(diag(1,node)-rho*W)%*%X)
  count = 0
  countBt = 0
  i = 2
  while(i<=num){
    
    Rho[i] = Rho[i-1]+rnorm(1, mean = 0, sd = 0.25)
    
    if (Rho[i]>1/max(Lamda) | Rho[i]< 1/min(Lamda)){
      next
    }
    
    INV_i=solve(diag(1,node)-Rho[i]*W)
    INV_ii=solve(diag(1,node)-Rho[i-1]*W)
    
    cal <- function(g){
      y=((pnorm(INV_i[g,]%*%(X%*%betamt[,i-1]))/pnorm(INV_ii[g,]%*%(X%*%betamt[,i-1])))^Y[g])*(((1-pnorm(INV_i[g,]%*%(X%*%betamt[,i-1])))/(1-pnorm(INV_ii[g,]%*%(X%*%betamt[,i-1]))))^(1-Y[g]))
      return(y)
    }
    
    proportion = prod(unlist(lapply(seq(1,node),cal)))*
      dtruncnorm(Rho[i-1], a=1/min(Lamda), b=1/max(Lamda), mean = Rho[i], sd = 0.25)/dtruncnorm(Rho[i], a=1/min(Lamda), b=1/max(Lamda), mean = Rho[i-1], sd = 0.25)
    
    proportion = min(1, Re(proportion))
    
    if (runif(1) > proportion){
      Rho[i]=Rho[i-1]
      count = count+1
    }
    
    betamt[,i] = betamt[,i-1]+rmvnorm(1, mean = rep(0,4), sigma=0.75*V)
    
    INV_iii=solve(diag(1,node)-Rho[i]*W)
    
    cal2 <- function(g){
      y=((pnorm(INV_iii[g,]%*%(X%*%betamt[,i]))/pnorm(INV_iii[g,]%*%(X%*%betamt[,i-1])))^Y[g])*(((1-pnorm(INV_iii[g,]%*%(X%*%betamt[,i])))/(1-pnorm(INV_iii[g,]%*%(X%*%betamt[,i-1]))))^(1-Y[g]))
      return(y)
    }
    
    proportionbeta = prod(unlist(lapply(seq(1,node),cal2)))
    
    proportionbeta = min(1, Re(proportionbeta))
    
    if (runif(1) > proportionbeta){
      betamt[,i]=betamt[,i-1]
      countBt = countBt+1
    }
    
    i = i+1
  }
  
  medianbeta = apply(betamt[,-(1:2000)],MARGIN=1, FUN=median)
  
  medianbeta_m[,k] = medianbeta
  medianRho_v[k] = median(Rho[-(1:2000)])
  
  Rate<-(num-count-1)/(num-1)
  RateB<-(num-countBt-1)/(num-1)
  
  print(Rate)
  print(RateB)
  
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

mean(medianRho_v)-rho

sum((medianRho_v-rho)^2)/itl

coverage/itl

mean(lbv)

mean(ubv)

mean(ubv)-mean(lbv)

