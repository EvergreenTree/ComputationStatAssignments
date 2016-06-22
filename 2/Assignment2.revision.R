###################################################
##1.observation probability
##  forward algorithm
#   description: given Markov Model parameters,
#     compute the prob of observation x with hidden statuses s at the final step
#   input:  x     observation
#           pi0   initial state prob
#           T     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#   output: P(x(1),...x(step),s(step))
forward <- function(x, pi0, T, E, step = length(x)){
  if (step > length(x))
    stop("step is greater than length of x")
  f = pi0 * E[,x[1]]
  for (i in 2:step)
    f = E[,x[i]]*(f%*%T)
  f
}
##  evaluation algorithm
#   input:  x     observation
#           pi0   initial state prob
#           T     Transmission matrix
#           E     Emission matrix
#   output: P(x)  probability of observation
evaluation <- function(x, pi0, T, E){
  sum(forward(x, pi0, T, E, length(x)))
}
##  apply algorithm above
x = c(1,2,2,2)
pi0 = c(0.5,0.5)
T = matrix(c(0.6,0.4,0.4,0.6),2)
E = matrix(c(0.5,0.8,0.5,0.2),2)
evaluation(x,pi0,T,E)

###################################################
##2.hidden state probability
##  backward algorithm
#   input:  x     observation
#           pi0   initial state prob
#           T     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#   output: P(x(i+1),...|s(i))
backward <- function(x,pi0,T,E,step = length(x)){
  n = length(x)
  if(step > n)
    return("errror:step is greater than length of x")
  if(step == n)
    return(c(1,1))
  f = t(T %*% E[,x[n]])
  i = n-1
  while(i >= step){
    f = f*t(T%*%E[,x[i]])
    i=i-1
  }
  f
}
##  forward-backward algorithm
#   input:  x     observation
#           pi0   initial state prob
#           T     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#   output: P(x,s(step))
forward_backward <-function(x,pi0,T,E,step = length(x)){
  forward(x,pi0,T,E,step)*backward(x,pi0,T,E,step)
}
##  apply algorithm,
lapply(1:4,forward_backward,x=x,pi0=pi0,T=T,E=E)

###################################################
##3.hidden-state probability
##  forward-backward algorithm
#   input:  x     observation
#           pi0   initial state prob
#           T     Transmission matrix
#           E     Emission matrix
#           s     hidden state
#   output: P(s|x)
decoding <- function(x,pi0,T,E,s){
  f = function(i){forward_backward(x,pi0,T,E,i)[s[i]]}
  exp(sum(log(sapply(1:length(x),f))))
}
##  apply algorithm,
s = c(1,1,2,2)
decoding(x,pi0,T,E,s)

###################################################
##4.decode optimal hidden state route
##  forward-backward algorithm
#   discription:  optimized state sequence given observed sequence
#   input:  x     observation
#           pi0   initial state prob
#           T     Transmission matrix
#           E     Emission matrix
#   output: argmax_s{P(s|x)}
viterbi <- function(x,pi0,T,E){
  L = length(x) ##length of observed sequence 
  S = dim(T)[1] ##S states in total
  s_ = rep(0,L) #initialize
  v_ = matrix(0,L+1,S)##v_k(i)=max_(s_1,...s_i-1)P(s_1,...s_i=k,x_1,...x_i)
  v_[1,1] = 1;
  for(i in 1:L){
    m = which.max(v_[i,]) ##max(v)=v[m]
    v_[i+1,] = v_[i,m] * T[m,] * E[,x[i]]
  }
  s_[L] = which.max(v_[L+1,])
  i=L
  while(i > 1){
    s_[i-1]=which.max(v_[i,] * T[,s_[i]])
    #print(v_[i,] * T[,s_[i]])
    i = i-1
  }
  s_
}
##  apply algorithm
viterbi(x,pi0,T,E)

###################################################
##4.learn optimal hidden state route
#   input:  x     observation
#   output: argmax_s{P(s|x)}
hmm_em = function(x,pi0,T,E){
  #to be updated
}

data = read.table("~/Desktop/STAT/hw/2/assign2.csv",sep = ",",header = TRUE)
data[1] = NULL
data = as.matrix(data)
data[data == "L"] = 1
data[data == "R"] = 2
data = matrix(as.numeric(data),ncol = 2)
hmm_em(data,pi0,T,E)
