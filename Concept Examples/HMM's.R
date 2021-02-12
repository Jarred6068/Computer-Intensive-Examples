

#generate observations from the hidden markov model described in question 1
#===============================function====================================================
sim.proteinHMM=function(T.mat=NULL, emission.mat=NULL, state.probs=NULL, sample.size=NULL){
  
  hid.seq=NULL
  emi.seq=NULL
  
  r=runif(1,0,1)
  
  hid.seq[1]=ifelse(r>state.probs[1], 1, 0)
  
  #-----------generate-hidden----------------
  for(i in 2:(sample.size)){
    r=runif(1,0,1)
    
    
    if(hid.seq[i-1]==1){
      
      tprobs=T.mat[2,]
      hid.seq[i]=ifelse(r<tprobs[1], 0, 1)
      
    }else{
      
      tprobs=T.mat[1,]
      hid.seq[i]=ifelse(r<tprobs[1], 0, 1)
      
    }
  }
  #-----------generate-observed--------------
  for(j in 1:(sample.size)){
    
    r=runif(1,0,1)
    
    if(hid.seq[j]==1){
      
      eprobs=emission.mat[2,]
      if(r<=eprobs[1]){emi.seq[j]="L"}
      if(r>eprobs[1] & r<=sum(eprobs[1:2])){emi.seq[j]="M"}
      if(r>sum(eprobs[1:2])){emi.seq[j]="H"}
      
    }else{
      
      eprobs=emission.mat[1,]
      if(r<=eprobs[1]){emi.seq[j]="L"}
      if(r>eprobs[1] & r<=sum(eprobs[1:2])){emi.seq[j]="M"}
      if(r>sum(eprobs[1:2])){emi.seq[j]="H"}
      
    }
    #------------------------------------------
    
  }
  
  #print(length(hid.seq))
  #print(length(emi.seq))
  HMM.sample=cbind.data.frame(hid.seq, emi.seq)
  colnames(HMM.sample)=c("Hidden", "Emission")
  return(HMM.sample)
  
}

#====================================================================


T.ij=matrix(c(0.7,0.3,0.4,0.6), nrow = 2, ncol = 2, byrow = T)
colnames(T.ij)=c("zero","one")

Pi=c((4/7),(3/7))

Gamma=matrix(c(0.7,0.2,0.1, 0.2, 0.4, 0.4), nrow = 2, ncol = 3, byrow=T)
colnames(Gamma)=c("L","M","H")

set.seed(239)

dat=sim.proteinHMM(T.ij, Gamma, Pi, 10)
dat2=sim.proteinHMM(T.ij, Gamma, Pi, 10000)
dat2$Emission=as.factor(dat2$Emission)
library(knitr)
kable(t(dat), caption="10 simulated values from Protein/expression HMM")

hid.state.obs=c(sum(dat2$Hidden)/length(dat2$Hidden), 1-(sum(dat2$Hidden)/length(dat2$Hidden)))
names(hid.state.obs)=c("P(X=1)=", "P(X=0)=")

kable(hid.state.obs, caption = "observed probability of the hidden states for a sample of 10,000 simulated observations")






library(HMM)
#using the HMM package to verify results
hmm1=initHMM(c("0","1"), c("L","M","H"), as.vector(Pi), as.matrix(T.ij), as.matrix(Gamma))













#Forward Pass Algorithm

Forward.pass=function(obs.seq=NULL, T.mat=NULL, E.mat=NULL, state.probs=NULL){
  #Emission matrix (E.mat) and Transition matrix (T.mat) must have columns labeled according
  # to the respective target states.
  pi.p=state.probs
  steps=length(obs.seq)
  
  if(obs.seq[1]=="L"){a.0i=diag(pi.p%*%t(E.mat[,1]))}
  if(obs.seq[1]=="M"){a.0i=diag(pi.p%*%t(E.mat[,2]))}
  if(obs.seq[1]=="H"){a.0i=diag(pi.p%*%t(E.mat[,3]))}
  
  a.tj=matrix(0, nrow = steps, ncol = 2)
  a.tj[1,]=a.0i
  
  for(t in 2:steps){
    
    if(obs.seq[t]=="L"){
      
      a.tj[t,]=diag(E.mat[,1])%*%t(T.mat)%*%a.tj[t-1,]
      
    }
    
    if(obs.seq[t]=="M"){
      
      a.tj[t,]=diag(E.mat[,2])%*%t(T.mat)%*%a.tj[t-1,]
      
    }
    
    if(obs.seq[t]=="H"){
      
      a.tj[t,]=diag(E.mat[,3])%*%t(T.mat)%*%a.tj[t-1,]
      
    }
    
    
  }
  
  
  
  return(a.tj)
  
  
}




seq.probs=Forward.pass(obs.seq=dat$Emission, T.mat=T.ij, E.mat=Gamma, state.probs=Pi)

output=cbind(seq.probs, rowSums(seq.probs))
colnames(output)=c("P(Y=y|Xt=0)", "P(Y=yt|Xt=1)", "P(Y=y|Theta)" )

library(knitr)
kable(output, caption = "Forward Algorithm results for simulated sequence of 10 observations")


expect=forward(hmm1, as.vector(as.character(dat$Emission)))

kable(exp(expect), caption = "Expected probabilities using the HMM package", digits = 6)






























