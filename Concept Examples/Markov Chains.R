

#a simple random walk example
X=NULL
X[1]=0
nsteps=10000

for(i in 2:nsteps){
  
  X[i]=X[i-1]+rnorm(1)
  
}

plot(X, type = "l", lwd=2)
abline(h=0, col="red", lty="dashed", lwd=2.5)



