#Example of MCMC 

#simulate from a gamma(a,b) distribution with normal proposal

sampler=NULL
n=100000
alpha=2
beta=0.5
stdev=1
init.value=alpha/beta
sampler[1]=init.value

for(i in 2:n){
  
  xt=sampler[i-1]
  Yt=rnorm(1, mean = xt, sd=stdev)
  ro=min( (dgamma(Yt, shape = alpha, rate = beta)/dgamma(xt, shape = alpha, rate = beta)) * 
            (dnorm(xt, mean = Yt, sd=stdev)/dnorm(Yt, mean = xt, sd=stdev) ), 1)
  #print(ro)
  
  r=runif(1)
  if(r<=ro){
    
    sampler[i]=Yt
    
    
  }else{
    
    sampler[i]=xt
    
  }
  
}


plot(c(1:n), sampler, type = "l", 
     xlab = "iteration",
     ylab = "sample value",
     main = "MCMC Gamma with Normal Proposal")
abline(h=alpha/beta, col="red", lwd=3, lty="dashed")
acf(sampler, main="ACF sample")











#simulate from a beta(a,b) distribution with normal proposal

sampler=NULL
n=100000
alpha=7
beta=7
stdev=1
init.value=alpha/(alpha+beta)
sampler[1]=init.value

for(i in 2:n){
  
  xt=sampler[i-1]
  Yt=rnorm(1, mean = xt, sd=stdev)
  ro=min( (dbeta(Yt, shape1 = alpha, shape2 = beta)/dbeta(xt, shape1 = alpha, shape2 = beta)) * 
            (dnorm(xt, mean = Yt, sd=stdev)/dnorm(Yt, mean = xt, sd=stdev) ), 1)
  #print(ro)
  
  r=runif(1)
  if(r<=ro){
    
    sampler[i]=Yt
    
    
  }else{
    
    sampler[i]=xt
    
  }
  
}


plot(c(1:n), sampler, type = "l", 
     xlab = "iteration",
     ylab = "sample value",
     main = "MCMC Beta with Normal Proposal")
abline(h=alpha/(alpha+beta), col="red", lwd=3, lty="dashed")
acf(sampler, main="ACF sample")














