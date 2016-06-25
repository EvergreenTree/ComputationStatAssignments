###################################################
##1.observation probability
##  forward algorithm
#   description: given Markov Model parameters,
#     compute the prob of observation x with hidden statuses s at the final step
#   input:  x     observation
#           pi0   initial state prob              ##question here
#           T     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#   output: P(x(1),...x(step),s(step))
forward <- function(x, pi0, T, E, step = length(x)){
  if (step > length(x))
    stop("step is greater than length of x")
  f = pi0 * E[,x[1]]
  if(step > 1)              # R language see 2:1 as c(2,1) instead of NULL
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
    stop("step is greater than length of x")
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
#           state s (if not specified,return a vector)
#   output: P(x,s(step))
forward_backward <-function(x,pi0,T,E,step = length(x),state = NULL){
  if(is.null(state))
    forward(x,pi0,T,E,step)*backward(x,pi0,T,E,step)
  else
    forward(x,pi0,T,E,step)[state]*backward(x,pi0,T,E,step)[state]
}
##  apply algorithm, compute P(s(step))
outcome = sapply(1:length(x),forward_backward,x=x,pi0=pi0,T=T,E=E)
outcome = outcome/(t(t(rep(1,2)))%*%apply(outcome,2,sum)) #normalize outcome by row
colnames(outcome) = c(1,2,3,4);rownames(outcome) = c("F","B")
outcome

###################################################
##3.hidden-state probability
##  decoding algorithm
#   input:  x     observation
#           pi0   initial state prob
#           T     Transmission matrix
#           E     Emission matrix
#           s     hidden state
#   output: P(s|x)
decoding <- function(x,pi0,T,E,s){
  exp(sum(log(sapply(1:length(x),function(i){
    forward_backward(x,pi0,T,E,i)[s[i]]
    }))))
}
##  apply algorithm,
s = c(1,1,2,2)
decoding(x,pi0,T,E,s)

###################################################
##4.decoding: optimal hidden state path (given model)
##  forward-backward algorithm
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
##4.learning: optimal hidden state path (without model)
##  forward-backward algorithm for double hidden states
#   input:  x     observation
#           pi0   initial state prob
#           T     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#           k,l   value of s(step-1), s(step)
#   output: P(x,s(step-1) = k, s(step) = l)
forward_backward_double_state <-function(x,pi0,T,E,step = length(x),k,l){
  forward(x,pi0,T,E,step-1)[k]*backward(x,pi0,T,E,step)[l]*T[k,l]*E[k,x[step]]
}
##  EM algorithm for learning HMM states
#   input:  n                         iteration step
#           x                         observation sequence, could be a list or a data frame
#           d = dim(table(x))         HMM dimension
#           pi0 = rep(1,2)            init value
#           T = diag(1,2)             init value
#           E = diag(1,2)             init value
#   output: argmax_s{P(s|x)}                 
hmm_em = function(n, x, d = dim(table(x)), pi0 = rep(0.5,d), T = diag(1,d), E = diag(1,d)){
  pi0 = forward_backward(x,pi0,T,E,1)
  L = length(x)
  for(i in 1:n){
    for(k in 1:d){
      for(l in 1:d){
        T[k,l] = sum(sapply(2:L, forward_backward_double_state,x=x,pi0=pi0,T=T,E=E,k=k,l=l))
        E[k,l] = sum(sapply(1:L, 
                            function(step)
                              forward_backward(x=x,pi0=pi0,T=T,E=E,step=step,state = k)*(l==x[step])
                            )
                     )#####question here
      }
    }
    T = T/(t(t(rep(1,d)))%*%apply(T,2,sum)) #normalize T by row
    E = E/(t(t(rep(1,d)))%*%apply(E,2,sum)) #normalize E by row
  }
  list(T = T,E = E)
}
dat = read.table("~/Desktop/STAT/hw/2/assign2.csv",sep = ",",header = TRUE)
dat[1] = NULL;dat = as.matrix(dat);dat[dat == "L"] = 1;dat[dat == "R"] = 2
dat = matrix(as.numeric(dat),ncol = 2)
dat[1:5,]
hmm_em(10,as.vector(dat))
#did not converge