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

if (!require(tmvtnorm)) {
  install.packages("tmvtnorm")
  library(tmvtnorm)
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
  detamt = matrix(0, node, num)
  
  Rho[1]=rho
  betamt[,1]=BBeta
  detamt[,1]=DDelta
  ####
  count = 0
  i = 2
  while(i<=num){
    Rho[i]=rtruncnorm(1, a=1/min(Lamda), b=1/max(Lamda), mean=(t(detamt[,i-1])%*%t(W)%*%(detamt[,i-1]-X%*%betamt[,i-1]))/(sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1]),
                      sd=sqrt(1/(sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1])))
    
    proportion = abs(det(diag(1,node)-Rho[i]*W)/det(diag(1,node)-Rho[i-1]*W))*exp(
      -1/2*t((diag(1,node)-Rho[i]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])%*%((diag(1,node)-Rho[i]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])+
        1/2*t((diag(1,node)-Rho[i-1]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])%*%((diag(1,node)-Rho[i-1]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])
      -1/2*((Rho[i-1]-(t(detamt[,i-1])%*%t(W)%*%(detamt[,i-1]-X%*%betamt[,i-1]))/(sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1]))^2/(1/(sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1])))+
        1/2*((Rho[i]-(t(detamt[,i-1])%*%t(W)%*%(detamt[,i-1]-X%*%betamt[,i-1]))/(sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1]))^2/(1/(sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1])))
    )
    
    proportion = min(1, Re(proportion))
    
    if (runif(1) > proportion){
      Rho[i]=Rho[i-1]
      count = count+1
    }
    
    betamt[,i] = rmvnorm(1, mean = solve(t(X)%*%X)%*%t(X)%*%(diag(1,node)-Rho[i]*W)%*%detamt[,i-1],
                         sigma = solve(t(X)%*%X))
    
    tmpmean=as.vector((solve(diag(1,node)-Rho[i]*W))%*%(X%*%betamt[,i]))
    pretrix = t(diag(1,node)-Rho[i]*W)%*%(diag(1,node)-Rho[i]*W)
    
    lw=c()
    up=c()
    
    for(j in 1:node){
      
      if(Y[j]==1){
        lw=c(lw,0)
        up=c(up,Inf)
      }else{
        lw=c(lw,-Inf)
        up=c(up,0)
      }
      
    }
    
    detamt[,i]=rtmvnorm(1, mean = tmpmean, H = pretrix, lower=lw, upper=up, algorithm="gibbs", start.value=detamt[,i-1])
    i = i+1
    
  }
  
  medianbeta = apply(betamt[,-(1:2000)],MARGIN=1, FUN=median)
  
  medianbeta_m[,k] = medianbeta
  medianRho_v[k] = median(Rho[-(1:2000)])
  
  Rate<-(num-count-1)/(num-1)
  
  print(Rate)
  
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


