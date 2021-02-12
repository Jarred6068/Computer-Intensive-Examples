

#binomial Distribution LLN

LLNbinom<-function(S,p,n){
  #S should be a vector of sample sizes
  #p is the probability of the outcome
  #n is the number of trials in the experiment

  avg=rep(0, length(S))
  avg2=rep(0, length(S))
  expected = n*p
  expected2=rep(0, length(S))
  sampler=rep(0, length(S))
  
  
  for (i in 1:length(S)){
    
    #expected[i]=n*p
    avg[i]=sum(rbinom(S[i], n, p))/S[i]
    
  }
  
  for (j in 1:length(S)){
    sampler[j]=rbinom(1, S[j], p)
    avg2[j]=sum(sampler[j])/S[j]
    expected2=S[j]*p
  }
  
  plot(S, avg, type = "l", col="blue")
  #lines(x=n, y = expected, type = "l", col = "red")
  abline(h=expected, col = "red")
  
  plot(S, avg2, type="l", col = "green")
  
  return(expected2)
  lines(S,expected2, col="red")
  
  return(expected2)
  #return(avg2)
}