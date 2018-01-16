###################################################
##1.observation probability
##  forward algorithm
#   description: given Markov Model parameters,
#     compute the prob of observation x with hidden statuses s at the final step
#   input:  x     observation
#           pi0   initial state prob              ##question here
#           P     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#   output: P(x(1),...x(step),s(step))
forward <- function(x, pi0, P, E, step = length(x)){
  if (step > length(x))
    stop("step is greater than length of x")
  f = pi0 * E[,x[1]]
  if(step > 1)              # R language see 2:1 as c(2,1) instead of NULL
    for (i in 2:step)
      f = E[,x[i]]*(f%*%P)
  f
}
##  evaluation algorithm
#   input:  x     observation
#           pi0   initial state prob
#           P     Transmission matrix
#           E     Emission matrix
#   output: P(x)  probability of observation
evaluation <- function(x, pi0, P, E){
  sum(forward(x, pi0, P, E, length(x)))
}
##  apply algorithm above
x = c(1,2,2,2)
pi0 = c(0.5,0.5)
P = matrix(c(0.6,0.4,0.4,0.6),2)
E = matrix(c(0.5,0.8,0.5,0.2),2)
evaluation(x,pi0,P,E)

###################################################
##2.hidden state probability
##  backward algorithm
#   input:  x     observation
#           pi0   initial state prob
#           P     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#   output: P(x(i+1),...|s(i))
backward <- function(x,pi0,P,E,step = length(x)){
  n = length(x)
  if(step > n)
    stop("step is greater than length of x")
  if(step == n)
    return(c(1,1))
  f = t(P %*% E[,x[n]])
  i = n-1
  while(i >= step){
    f = f*t(P%*%E[,x[i]])
    i=i-1
  }
  f
}
##  forward-backward algorithm
#   input:  x     observation
#           pi0   initial state prob
#           P     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#           state s (if not specified,return a vector)
#   output: P(x,s(step))
forward_backward <-function(x,pi0,P,E,step = length(x),state = NULL){
  if(is.null(state))
    forward(x,pi0,P,E,step)*backward(x,pi0,P,E,step)
  else
    forward(x,pi0,P,E,step)[state]*backward(x,pi0,P,E,step)[state]
}
##  apply algorithm, compute P(s(step))
outcome = sapply(1:length(x),forward_backward,x=x,pi0=pi0,P=P,E=E)
outcome = outcome/(t(t(rep(1,2)))%*%apply(outcome,2,sum)) #normalize outcome by row
colnames(outcome) = c(1,2,3,4);rownames(outcome) = c("F","B")
outcome

###################################################
##3.hidden-state probability
##  decoding algorithm
#   input:  x     observation
#           pi0   initial state prob
#           P     Transmission matrix
#           E     Emission matrix
#           s     hidden state
#   output: P(s|x)
decoding <- function(x,pi0,P,E,s){
  exp(sum(log(sapply(1:length(x),function(i){
    forward_backward(x,pi0,P,E,i)[s[i]]
    }))))
}
##  apply algorithm,
s = c(1,1,2,2)
decoding(x,pi0,P,E,s)

###################################################
##4.decoding: optimal hidden state path (given model)
##  forward-backward algorithm
#   input:  x     observation
#           pi0   initial state prob
#           P     Transmission matrix
#           E     Emission matrix
#   output: argmax_s{P(s|x)}
viterbi <- function(x,pi0,P,E){
  L = length(x) ##length of observed sequence 
  S = dim(P)[1] ##S states in total
  s_ = rep(0,L) #initialize
  v_ = matrix(0,L+1,S)##v_k(i)=max_(s_1,...s_i-1)P(s_1,...s_i=k,x_1,...x_i)
  v_[1,1] = 1;
  for(i in 1:L){
    m = which.max(v_[i,]) ##max(v)=v[m]
    v_[i+1,] = v_[i,m] * P[m,] * E[,x[i]]
  }
  s_[L] = which.max(v_[L+1,])
  i=L
  while(i > 1){
    s_[i-1]=which.max(v_[i,] * P[,s_[i]])
    #print(v_[i,] * P[,s_[i]])
    i = i-1
  }
  s_
}
##  apply algorithm
viterbi(x,pi0,P,E)

###################################################
##4.learning: optimal hidden state path (without model)
##  forward-backward algorithm for double hidden states
#   input:  x     observation
#           pi0   initial state prob
#           P     Transmission matrix
#           E     Emission matrix
#           step  position of markov chain to compute
#           k,l   value of s(step-1), s(step)
#   output: P(x,s(step-1) = k, s(step) = l)
forward_backward_double_state <-function(x,pi0,P,E,step = length(x),k,l){
  if(step <= 1) return(0)
  forward(x,pi0,P,E,step-1)[k]*backward(x,pi0,P,E,step)[l]*P[k,l]*E[l,x[step]]
}
##  EM algorithm for learning HMM states
#   input:  n                         iteration step
#           x                         observation sequence, could be a list or a data frame
#           d = dim(table(x))         HMM dimension
#           pi0 = rep(1,2)            init value
#           P = diag(1,2)             init value
#           E = diag(1,2)             init value
#   output: parameter, we can obtain argmax_s{P(s|x)} with algorithm above
hmm_em = function(n, x, d = dim(table(x)), pi0 = rep(0.5,d), P = diag(1,d), E = diag(1,d)){
  L = length(x);P0 = matrix(0,d,d);E0 = matrix(0,d,d)
  for(i in 1:n){
    pi0 = forward_backward(x,pi0,P,E,1)
    pi0 = pi0/sum(pi0)
    for(k in 1:d){
      for(l in 1:d){
        P0[k,l] = sum(sapply(2:L, forward_backward_double_state,x=x,pi0=pi0,P=P,E=E,k=k,l=l))
        E0[k,l] = sum(sapply(1:L, 
                             function(step)
                               forward_backward(x=x,pi0=pi0,P=P,E=E,step=step,state = k)*(x[step]==l)
                            )
                     )#####question here
        
      }
    }
    P = P0/(t(t(rep(1,d)))%*%apply(P,2,sum)) #normalize P by row
    E = E0/(t(t(rep(1,d)))%*%apply(E,2,sum)) #normalize E by row
  }
  list(P = P,E = E)
}
#rewrite with circulation
hmm_em_slow <- function(n, x, d = dim(table(x)), pi0 = rep(0.5,d), P = matrix(0.5,2,2), E = matrix(0.5,2,2)){
  L = length(x)
  for(i in 1:n){
    list(pi0 = pi0,P = P,E = E)
    pi0 = forward_backward(x,pi0,P,E,1)
    pi0 = pi0 / sum(pi0)
    P0 = matrix(0,d,d); E0 = matrix(0,d,d)
    for(k in 1:d){
      for(l in 1:d){
        for(j in 1:L){
          P0[k,l] = P0[k,l]+forward_backward_double_state(x=x,step=j,pi0=pi0,P=P,E=E,k=k,l=l); t(as.matrix(list(j = j,P = forward_backward_double_state(x=x,step=j,pi0=pi0,P=P,E=E,k=k,l=l)))
          E0[k,l] = E0[k,l]+forward_backward(x=x,pi0=pi0,P=P,E=E,step=j,state = k)*(x[j]==l)
        }
      }
    }
    P = P0/(t(t(rep(1,d)))%*%apply(P0,2,sum)) #normalize P by row
    E = E0/(t(t(rep(1,d)))%*%apply(E0,2,sum)) #normalize E by row
  }
  list(pi0 = pi0, P = P, E = E)
}

dat = read.table("~/Desktop/STAT/hw/2/assign2.csv",sep = ",",header = TRUE)
dat[1] = NULL; dat = as.matrix(dat);dat[dat == "L"] = 1; dat[dat == "R"] = 2
dat = matrix(as.numeric(dat),ncol = 2)
dat[1:5,]
#hmm_em(3,as.vector(dat))
#did not converge
parameter = hmm_em_slow(3,as.vector(dat))
#optimal hidden-state route
viterbi(dat,parameter$pi0,parameter$P,parameter$E)