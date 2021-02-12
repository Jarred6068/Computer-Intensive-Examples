#(1)--Monte-Carlo-Integration-Example-Problem-One:---

#We'll start with basic integration. Suppose we have an instance of a Normal distribution with a mean of 1 and a standard deviation of 10.
#Then we want to find the integral from {3<=X<=6} N(x;1,10) as visualized below

#firstly, use what you know about random number generation to obtain a sample for Monte Carlo Integration
#=============================================

n=10000
#Box-Muller Transformation
u1=runif(n,0,1)
u2=runif(n,0,1)

z1=sqrt(-2*log(u1))*cos(2*pi*u2)
z2=sqrt(-2*log(u1))*sin(2*pi*u2)

#check to make sure normal:
par(mfrow=c(2,1))
hist(z1)
hist(z2)
par(mfrow=c(1,1))

#linear transformation of z1 to normal needed for problem context
#recall: x~N(mu, sigma) and Y=aX+b, then Y~N(a*mu+b, a^2*sigma^2)
xx1=10*z1+1
MC.val=NULL

for(i in 1:length(xx1)){
  
  MC.val[i]=mean(xx1[1:i]>=3&xx1[1:i]<=6)
  
}
print(paste("Integral.Value=",MC.val[length(xx1)]))

plot(1:length(xx1), MC.val, type="l", col="red", lwd=2.5,
     main="integral for 3<=x<=6 of X~N(1,10)",
     xlab="Sample Size",
     ylab="value of Integral")







