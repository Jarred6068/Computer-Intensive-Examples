
#by newton raphson algorithm 

#initialize data
#sample size of the data we observed
n=1000
#the number of times to run NR
its=10
#the true value of the parameter
lambda.true=6
#generate the data we observe
x.vec=rexp(n, lambda.true)
#enter a starting guess
starting.guess=0.5
#allocate space
Lam.n=rep(0,its)
#give the storage the starting guess
Lam.n[1]=starting.guess

#first derivative of log-likelihood
d1.dl=function(lambda=NULL, ss=NULL, data=NULL){
  n/lambda-sum(x.vec)
}

#second derivative of log-likelihood
d2.dl=function(lambda=NULL, ss=NULL){
  -n/lambda^2
}

#the newton raphson algoritm
for(i in 1:its){

Lam.n[i+1]=Lam.n[i]-(d1.dl(lambda = Lam.n[i], ss=n, data = x.vec)/d2.dl(lambda = Lam.n[i], ss=n))
}

#we plot the final value retuned by newton raphson
plot(c(1:(its+1)), Lam.n, 
     type = "p", 
     xlab ="iteration", 
     ylab="Value of Lambda", 
     main="Values of Lambda at Each Iteration",
     pch=19,
     col="blue")
#plot the line of the theoretical MLE
abline(h=1/mean(x.vec), col="red")

NR.estimate=Lam.n[length(Lam.n)]
print(NR.estimate)

MLE.Estimate=1/mean(x.vec)
print(MLE.Estimate)
