
# let the RV x represent shortest the distance between the center of the needle
# and a line
# let the RV theta represent the acute angle formed by the needle and the closest line

# let both RV's be uniformly distributed such that 0<=X<=T/2 and 0<=theta<=L/2

# let L represent the length of the needle and T represent the distance between two lines
# additionally let L==T

# we now have that the needle crosses a line whenever x<=L/2
# because X and theta are independent the joint density is the product of their 
# individual densities. 

#===========================================================
#-----------Monte--Carlo--Method?---------------------------
#===========================================================
start.time=Sys.time()
#set the number of resamplings
iter=10000
#allocation of space...
prob=rep(0,iter)

#-----resample loop----
for(i in 1:iter){
#set/select random length between lines
T=runif(1,min=0,max=2)
L=T
#set sample size
n=iter/10
#generate sample of random variables
theta=runif(n,min = 0, max=pi/2)
X=runif(n, min = 0, max = T/2)

#allocate space
P=rep(1, n)
stat=rep(0, n)
#calculate the line cross criteria
stat=(L/2)*sin(theta)
#find the number of needles that cross
logi=(X<=stat)
P=P[logi]
good=sum(P)
#calculate probability
prob[i]=good/n

}
end.time=Sys.time()
#----end--loop---------
#estimate pi from samples
est.pi=(2*L)/(T*prob)
mean(est.pi)

run.avg=cumsum(est.pi)/c(1:iter)
#remove unnecessary vars
rm(n, i, logi, P)

#plotting---------------------------------------


plot(1:iter, est.pi, type = "l", col="blue", xlab = "sample", ylab = "estimated pi", main=("Buffon's Needle to Est. Pi"))
abline(h=mean(est.pi), col="green", lwd=2.5)
abline(h=pi, col="red", lwd=2.5)

plot(1:iter, run.avg, type = "l", col="blue", xlab = "sample", ylab = "estimated pi", main=("Buffon's Needle to Est. Pi"))
abline(h=pi, col="red", lwd=2.5)

Run.Time=end.time-start.time
Run.Time
