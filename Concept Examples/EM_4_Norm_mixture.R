


#generate random values
n=400
#generate random mean vectors
mu1=c(runif(1,0,5), runif(1,3,9))
mu2=c(runif(1,3,7), runif(1,0,5))
mu3=c(runif(1,0,5), runif(1,7,12))
mu4=c(runif(1,1,9), runif(1,4,5))

means.list=list(mu1, mu2, mu3, mu4)

prop1=0.3
prop2=0.2
prop3=0.1
prop4=0.4

props.vec=c(prop1, prop2, prop3, prop4)

library(Matrix)
#generate random sigmas and find the nearest positive definite form:
sigma1=as.matrix(nearPD(matrix(c( runif(1,0,1), runif(1,-3,5), runif(1,-1,1), runif(1,5,6) ), ncol = 2), keepDiag = TRUE)$mat)
sigma2=as.matrix(nearPD(matrix(c( runif(1,2,10), runif(1,-1,1), runif(1,-5,9), runif(1,1,4) ), ncol = 2), keepDiag = TRUE)$mat)
sigma3=as.matrix(nearPD(matrix(c( runif(1,5,7), runif(1,-7,3), runif(1,-3,3), runif(1,2,4) ), ncol = 2), keepDiag = TRUE)$mat)
sigma4=as.matrix(nearPD(matrix(c( runif(1,1,9), runif(1,2,5), runif(1,2,4), runif(1,3,4) ), ncol = 2), keepDiag = TRUE)$mat)

sigmas.list=list(sigma1=sigma1, sigma2=sigma2, sigma3=sigma3, sigma4=sigma4)

#==================function==============================
Gen.mix.data=function(prop, varcov.mat, Samp.size, mu){
  
  n=Samp.size*prop
  z1=rnorm(n)
  z2=rnorm(n)
  
  Z=cbind(z1,z2)
  
  sig.chol=chol(varcov.mat)
  
  mixi.obs=t(t(sig.chol)%*%t(Z))+mu
  return(mixi.obs)
  
}
#========================================================

mix1.obs=Gen.mix.data(prop1, sigma1, n, mu1)
mix2.obs=Gen.mix.data(prop2, sigma2, n, mu2)
mix3.obs=Gen.mix.data(prop3, sigma3, n, mu3)
mix4.obs=Gen.mix.data(prop4, sigma4, n, mu4)

obs.list=list(mix1=mix1.obs, mix2=mix2.obs, mix3=mix3.obs, mix4=mix4.obs)

mix.mat=rbind(mix1.obs, mix2.obs, mix3.obs, mix4.obs)

#===================function=================================
calc.mix=function(observations=NULL, props=NULL, sigmas, means){
  
  library(mvtnorm)
  dat.obs=list()
  
  for(i in 1:4){dat.obs[[i]]=props[i]*dmvnorm(observations, mean = means[[i]], sigma = sigmas[[i]])}
  
  mixture.obs=Reduce("+", dat.obs)
  return(mixture.obs)
}
#============================================================

dens.mixmat=calc.mix(observations = mix.mat, props=props.vec, sigmas = sigmas.list, means=means.list)

library(MASS)


dens=kde2d(mix.mat[,1], mix.mat[,2], n=dim(mix.mat)[1])$z
image(dens)
contour(dens, add = TRUE)


library(knitr)
kable(as.data.frame(means.list), digits = 4, col.names = c("Mu.1", "Mu.2","Mu.3","Mu.4"),
      caption="True Mean vectors for jth mixture component")
kable(as.data.frame(sigmas.list), digits = 4,
      caption="True Variance-Covariance matrices for jth mixture component")
names(props.vec)=c("mix1","mix2","mix3","mix4")
kable(as.data.frame(props.vec), digits = 4, col.names = c("Proportion"),
      caption="True Mixing Proportions")



#functions needed to compute EM

#==================function=========================================================
MLE.mean=function(mix.data=NULL, T.ij.matrix=NULL, n=NULL){
  mu1k = matrix(0, nrow = n, ncol = 2)
  mu2k = matrix(0, nrow = n, ncol = 2)
  mu3k = matrix(0, nrow = n, ncol = 2)
  mu4k = matrix(0, nrow = n, ncol = 2)
  
  for(i in 1:n){
    
    mu1k[i,]=T.ij.matrix[i,1]%*%mix.data[i,]
    mu2k[i,]=T.ij.matrix[i,2]%*%mix.data[i,]
    mu3k[i,]=T.ij.matrix[i,3]%*%mix.data[i,]
    mu4k[i,]=T.ij.matrix[i,4]%*%mix.data[i,]
    
  }
  mu1k=colSums(mu1k)/sum(T.ij.matrix[,1])
  mu2k=colSums(mu2k)/sum(T.ij.matrix[,2])
  mu3k=colSums(mu3k)/sum(T.ij.matrix[,3])
  mu4k=colSums(mu4k)/sum(T.ij.matrix[,4])
  
  return(list(mu1=mu1k, mu2=mu2k, mu3=mu3k, mu4=mu4k))
}
#===================================================================================


#==================function=========================================================
MLE.sigma=function(mix.data = NULL, T.ij.matrix=NULL, updated.mu=NULL, n=NULL){
  
  sigma1k=list()
  sigma2k=list()
  sigma3k=list()
  sigma4k=list()
  
  
  for(i in 1:n){
    
    sigma1k[[i]]=T.ij.matrix[i,1]*((mix.data[i,]-updated.mu[[1]])%*%t(mix.data[i,]-updated.mu[[1]]))
    sigma2k[[i]]=T.ij.matrix[i,2]*((mix.data[i,]-updated.mu[[2]])%*%t(mix.data[i,]-updated.mu[[2]]))
    sigma3k[[i]]=T.ij.matrix[i,3]*((mix.data[i,]-updated.mu[[3]])%*%t(mix.data[i,]-updated.mu[[3]]))
    sigma4k[[i]]=T.ij.matrix[i,4]*((mix.data[i,]-updated.mu[[4]])%*%t(mix.data[i,]-updated.mu[[4]]))
  }
  
  sigma1k=Reduce("+", sigma1k)/sum(T.ij.matrix[,1])
  sigma2k=Reduce("+", sigma2k)/sum(T.ij.matrix[,2])
  sigma3k=Reduce("+", sigma3k)/sum(T.ij.matrix[,3])
  sigma4k=Reduce("+", sigma4k)/sum(T.ij.matrix[,4])
  
  return(list(sig1=sigma1k, sig2=sigma2k, sig3=sigma3k, sig4=sigma4k))
  
}


#=====================================================================================

#==============================function=====================================================
EM.4mix=function(iterations=NULL, start.props=NULL, start.mu=NULL, start.sigma=NULL){
  
  # E-step initialization
  library(mvtnorm)
  pi.est.vec=start.props
  mu.est.list=start.mu
  sigma.est.list=start.sigma
  T.mat=matrix(0, nrow = 400, ncol = 4)
  stor=NULL
  stor2=NULL
  
  prop.estimates=as.data.frame(matrix(0, nrow = iterations, ncol = 4))
  means.estimates=list()
  sigmas.estimates=list()
  
  for(k in 1:its){
    
    
    #calculate expectation
    for(i in 1:n){
      stor=NULL
      for(j in 1:4){
        T.ij=(dmvnorm(mix.mat[i,], mean = mu.est.list[[j]], sigma = sigma.est.list[[j]])*pi.est.vec[j])/
          calc.mix(mix.mat[i,], pi.est.vec, sigmas = sigma.est.list, means = mu.est.list)
        T.mat[i,j]=T.ij
        #print(T.ij)
        logs=log(pi.est.vec[j])-(1/2)*log(2*pi)-(1/2)*log(det(sigma.est.list[[j]]))
        #print(logs)
        kernel=-(1/2)*t(mix.mat[i,]-mu.est.list[[j]])%*%sigma.est.list[[j]]%*%(mix.mat[i,]-mu.est.list[[j]]) 
        #print(kernel)
        stor[j]=T.ij%*%(logs*kernel)
        
        stor2[i]=sum(stor)
        
      }
      
      #print(stor)
      
    }   
    
    #---------------------The-M-Step---------------------
    #update parameters 
    #mean vector MLE
    updated.mu=MLE.mean(mix.data = mix.mat, T.ij.matrix = T.mat, n= 400)
    mu.est.list=updated.mu
    means.estimates[[k]]=as.data.frame(updated.mu)
    #sigma.MLE
    updated.sigma=MLE.sigma(mix.data = mix.mat, T.ij.matrix = T.mat, updated.mu = updated.mu, n=400 )
    sigma.est.list=updated.sigma
    sigmas.estimates[[k]]=as.data.frame(updated.sigma)
    #proportion MLE
    pi.est.vec=colMeans(T.mat)
    prop.estimates[k,]=colMeans(T.mat)
  }  
  
  return(list(means=means.estimates, sigmas=sigmas.estimates, props=prop.estimates))
  
}  

#===================================================================================================  






#---------------------------run--the--EM---------------------------------------

its=500

mu0=list(c(1,1),c(0,0),c(2,2),c(3,3))
sigma0=list(matrix(c(1,0,0,1),ncol=2),matrix(c(1,0,0,1),ncol=2),matrix(c(1,0,0,1),ncol=2),matrix(c(1,0,0,1),ncol=2))
pi0=c(0.25,0.25,0.25,0.25)

start.time=Sys.time()
output=EM.4mix(iterations = its, start.props = pi0, start.mu = mu0, start.sigma = sigma0)
end.time=Sys.time()

Time.to.compute=end.time-start.time
print(paste("Time To Compute EM:", round(Time.to.compute,4),"min;", "Iterations:", its , sep = " "))
library(knitr)

kable(as.data.frame(output$means[[its]]), digits = 4, caption="Mean Estimates at final iteration of EM")
kable(as.data.frame(output$sigmas[[its]]), digits = 4, caption="Variance-Covariance Estimates at final iteration of EM")
kable(as.data.frame(output$props[its,]), digits = 4, caption="Mixing Proportion Estimates at final iteration of EM")

meanii.vec=function(out=NULL,its=NULL,mix.idx=NULL){
  mean.EM.estimate=matrix(0, nrow = its, ncol=2)
  for(i in 1:its){mean.EM.estimate[i,]=out$means[[i]][,mix.idx]}
  
  return(mean.EM.estimate)
}

means1.mat=meanii.vec(output, its, mix.idx=1)
means2.mat=meanii.vec(output, its, mix.idx=2)
means3.mat=meanii.vec(output, its, mix.idx=3)
means4.mat=meanii.vec(output, its, mix.idx=4)

x.axis=c(1:its)
par(mfrow=c(2,2))
plot(x.axis, means1.mat[,1], type = "l", lwd=2.5, col="blue2", ylim=c(min(means1.mat), max(means1.mat)),
     main=paste("Estimate mean vec1"))
lines(x.axis, means1.mat[,2], lwd=2.5, col="lightblue2")
plot(x.axis, means2.mat[,1], type = "l", lwd=2.5, col="blue2", ylim=c(min(means2.mat), max(means2.mat)),
     main=paste("Estimate mean vec2"))
lines(x.axis, means2.mat[,2], lwd=2.5, col="lightblue2")
plot(x.axis, means3.mat[,1], type = "l", lwd=2.5, col="blue2", ylim=c(min(means3.mat), max(means3.mat)),
     main=paste("Estimate mean vec3"))
lines(x.axis, means3.mat[,2], lwd=2.5, col="lightblue2")
plot(x.axis, means4.mat[,1], type = "l", lwd=2.5, col="blue2", ylim=c(min(means4.mat), max(means4.mat)),
     main=paste("Estimate mean vec4"))
lines(x.axis, means4.mat[,2], lwd=2.5, col="lightblue2")
par(mfrow=c(1,1))
