## Multi-dimentional Estimation-Maximization method for Mixed Gaussian Distribution of n terms(clusters):
library("mvtnorm")
EM = function(data, cluster_number) {
  #setting primitive values
  n = nrow(data)
  dim = ncol(data)
  mu = matrix(runif(cluster_number * dim),nrow = dim,ncol = cluster_number)##mean
  sigma = rep(c(diag(dim)), length = cluster_number * dim * dim);dim(sigma) = c(dim, dim, cluster_number)##covariance matrix
  tau = runif(cluster_number)##coefficients of mixed gaussian distribution
  
  step = 20
  for (s in 1:step) {
    #the E Step:
    T = matrix(nrow = cluster_number, ncol = n)##membership prob
    for (i in 1:n) {
      P = 0
      for (j in 1:cluster_number) {
        P = P + tau[j] * dmvnorm(data[i, ], mu[, j], sigma[, , j])
      }
      for (j in 1:cluster_number) {
        T[j, i] = tau[j] * dmvnorm(data[i, ], mu[, j], sigma[, , j]) / P
      }
    }
    #the M Step:
    tau = apply(T, 1, mean)
    sum = apply(T, 1, sum)
    for (j in 1:cluster_number) {
      P = rep(0, dim)
      Q = rep(0, dim * dim);dim(Q) = c(dim, dim)
      for (i in 1:n) {
        P = P + T[j, i] * data[i, ]
        Q = Q + T[j, i] * t(as.matrix(data[i, ] - mu[,j])) %*% as.matrix(data[i, ] - mu[,j])
      }
      mu[, j] = t(as.matrix(P / sum[j]))
      sigma[,,j] = as.matrix(Q / sum[j])
    }
  }
  return(list(mu=mu,sigma=sigma,tau=tau))
}
## Estimate coefficients in Data1.csv
data = read.table("~/Desktop/STAT/hw/1/Data1'.csv",sep = ",",header = TRUE)
EM(data,3)

