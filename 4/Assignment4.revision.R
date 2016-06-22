###################################################
## dependencies
library("matlab")
###################################################
## 1.Gibbs Sampling for biexponential distribution
gibbs_biexp <- function(niter,x0,y0,M){
  x = rep(x0,niter);y = rep(y0,niter)
  rexpM <- function(lambda,M){
    repeat{
      x = rexp(1,lambda)
      if(x < M) return(x)
    }
  }
  for(i in 2:niter){
    x[i] = rexpM(y[i-1],M)
    y[i] = rexpM(x[i],M)
  }
  result = cbind.data.frame(x,y)
}
#  assume B = 10, simulate x,y
plot(gibbs_biexp(niter = 1000,1,1,10))
#  Ex, Exy as a function of B
B = linspace(1,20,50);
f11 <- function(M)
  mean(gibbs_biexp(niter = 3000,1,1,M)[,1])
f12 <- function(M) {
  A = gibbs_biexp(niter = 3000,1,1,M)
  mean(A[,1]*A[,2])
}
Ex = sapply(B, f11)
Exy = sapply(B, f12)
plot(B,Ex)
plot(B,Exy)

###################################################
## 2.conditional expectation via simulation
#  x_1,2,3 ~ exp dist with mean 1 (satisfying iid). 
#  assume y = x_1 + 2*x_2 + 3*x_3, solve E1 = E(y|y>15); E2 = E(y|y<1)
f2 <- function(M,greater = TRUE){
  c = c(1,2,3);
  repeat{
    x = sum(rexp(3,1)*c)
    if(greater){if(x > M) return(x)}
    else if(x < M) return(x)
  }
}
E <- function(niter, bound, greater)
  mean(sapply(rep(bound,niter),f2,greater = greater))
E1 = E(10000,15,TRUE);
E2 = E(10000,1,FALSE);

###################################################
##3. f
triple_simulation <- function(niter, FUN){
    
}
f3 <- function(x,y,z){
  exp(-x-y-z-x*y-y*z-x*z)
}

###################################################
##4.
f4 <- function(niter, alpha = 2,beta = 3, lambda = 4){
  x = 
  for(i in 1:niter){
    
  }
  return(cbind.data.frame(x[floor(niter/10):niter],y[floor(niter/10):niter],n[floor(niter/10):niter]));
}

###################################################
## 5.Metropolis-Hasting Algorithm for 2-dim-GMM
#  algorithm: gmm
#  method:    Metropolis-Hasting
#  datatype:  x       data.frame(names = c("x1","x2"))
#             mu      list(names = c("a","b"))
#             sigma   list(names = c("a","b"))
#             sigma0  matrix
library(mvtnorm)
gmm_MH <- function(n, mu ,sigma, sigma0 = diag(2,2)){
  state = 1
  x = data.frame(x1=rep(-5,n),x2=rep(-5,n))
  for(i in 2:n){
    x[i,] = x[i-1,] + rmvnorm(1,c(0,0),sigma0)
    if(runif(1) > dmvnorm(x[i,],mu[[state]],sigma[[state]])/dmvnorm(x[i-1,],mu[[state]],sigma[[state]])){
      x[i,] = x[i-1,]
    }
    if(state == 1){if(runif(1) < 0.1) state = 2}
    else{if(runif(1) < 0.1) state = 1}
  }
  names(x)=c("x1",'x2')
  x
}  
mu = list(a=c(1,4),b=c(-2,-1))
sigma = list(a=matrix(c(1,0.3,0.3,2),2),b=matrix(c(3,0.4,0.4,1),2))
plot(gmm_MH(5000,mu,sigma))
