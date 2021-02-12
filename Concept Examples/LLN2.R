#ABOUT THE LLN:

# THE LLN STATES THAT THE SAMPLE MEAN CONVERGES IN PROBABILITY TO THE EXPECTED VALUE AS THE SAMPLE SIZE INCREASES.
# THE LAW HOLDS FOR MOST DISTRIBUTIONS (EXPECT THE CAUCHY AND PARETO DISTRIBUTIONS)
# IT IS IMPORTANT AND ALLOWS US TO ESTIMATE PROBABILITIES FROM UNKNOWN MASS FUNCTIONS
# IT FORMS THE BACKBONE OF THE CENTRAL LIMIT THEOREM


#Example 1: Bernoulli/Binomial Experiment  ---- coin flip

#expectation is p

LLNbern<-function(S,p,n){
  #S should be the max sample
  #p is the probability of the outcome
  #n is the number of trials in the experiment (for bernoulli = 1)
  #NOTE: when n > 1 this function transfers to a binomial exper.
  
  avg=rep(0, S)
  expected = n*p
  outcomes=cumsum(rbinom(S, n, p))
  trials=rep(0,S)
  
  for (i in 1:S){
    trials[i]=i
    avg[i]=outcomes[i]/i
    
  }
  
  plot(trials, avg, type = "l", col="blue",main = "Bernoulli LLN")
  legend(x=0.4*S, y=max(avg)+0.01, col = c("red","blue"), c("expected value", "running sample mean"), lty = c("solid", "solid"))
  abline(h=expected, col = "red")
  
  return(avg)
  
}
p=0.5
n=1
S=1000

sample.means.bern=LLNbern(S,p,n)


#=================================================================================================
#-------------------------------------------------------------------------------------------------
#=================================================================================================



#Example 2: Binomial Experiment

#expectation is n*p

LLNbinom<-function(S,p,n){
  #S should be the max sample
  #p is the probability of each outcome
  #n is the number of trials in the experiment (for bernoulli = 1)
  
  avg=rep(0, S)
  expected = n*p
  outcomes=cumsum(rbinom(S, n, p))
  trials=rep(0,S)
  
  for (i in 1:S){
    trials[i]=i
    avg[i]=outcomes[i]/i
    
  }
  
  plot(trials, avg, type = "l", col="blue",main = "Binomial LLN")
  legend(x=0.4*S, y=max(avg), col = c("red","blue"), c("expected value", "running sample mean"), lty = c("solid", "solid"))
  abline(h=expected, col = "red")
  
  return(avg)
  
}
p=runif(1,0,1)
n=10
S=1000

sample.means.binom=LLNbinom(S,p,n)


#=================================================================================================
#-------------------------------------------------------------------------------------------------
#=================================================================================================


#Example 3: poisson experiment

#expectation is lambda

LLNpois<-function(S,lambda){
  #S should be the max sample
  #lambda is the rate parameter of the poisson
  
  avg=rep(0, S)
  expected = lambda
  outcomes=cumsum(rpois(S, lambda))
  trials=rep(0,S)
  
  for (i in 1:S){
    trials[i]=i
    avg[i]=outcomes[i]/i
    
  }
  
  plot(trials, avg, type = "l", col="blue",main = "Poisson LLN")
  legend(x=0.4*S, y=max(avg), col = c("red","blue"), c("expected value", "running sample mean"), lty = c("solid", "solid"))
  abline(h=expected, col = "red")
  
  return(avg)
  
}
lambda=abs(round(rnorm(1,3,1)))
S=1000

sample.means.pois=LLNpois(S,lambda)


