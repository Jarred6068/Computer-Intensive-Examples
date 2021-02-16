

#Simulation and estimation of a hidden Markov model. Consider a certain protein that moves along the
#genome, binds to the genome at a few sites and falls off the genome. When it binds to the genome, it
#may help activate the gene where it binds such that the gene may have an elevated expression compared
#with when the protein is not bound to this gene. The protein may potentially bind to 10 genes that are
#located linearly along the genome. We observe the gene expression levels but not the protein binding
#tatus.
#The movement of the protein along the genome can be described by a Markov chain with two states: 0
#(off the genome) and 1 (bound to the genome). The initial distribution is


#pi = c(4/7, 3/7)

#The transition probability matrix is
#                       0    1
#                    0  0.7  0.3
#                    1  0.4  0.6

#A gene can have three levels of expression: low (L), medium (M) and high (H). The emission probability
#matrix going from the protein binding state to gene expression is

#                         L     M     H
#                     0   0.7   0.2   0.1
#                     1   0.2   0.4   0.4

#generate observations from the hidden markov model described:
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
print(dat)
print(dat2)

hid.state.obs=c(sum(dat2$Hidden)/length(dat2$Hidden), 1-(sum(dat2$Hidden)/length(dat2$Hidden)))
names(hid.state.obs)=c("P(X=1)=", "P(X=0)=")
print(hit.state.obs)






#=================================================================================================

library(HMM)
#using the HMM package to verify results
hmm1=initHMM(c("0","1"), c("L","M","H"), as.vector(Pi), as.matrix(T.ij), as.matrix(Gamma))


#=================================================================================================

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

#print(output)


expect=forward(hmm1, as.vector(as.character(dat$Emission)))

print(expect)
#================================================================================================



#Viterbi Algorithm

ViterbiA=function(obs.seq=NULL, T.mat=NULL, E.mat=NULL, state.probs=NULL){
  
  
  pi.p=state.probs
  T.ij=as.matrix(T.mat)
  Gamma=as.matrix(E.mat)
  hid=NULL
  steps=length(obs.seq)
  
  if(obs.seq[1]=="L"){g.0i=diag(pi.p%*%t(E.mat[,1]))}
  if(obs.seq[1]=="M"){g.0i=diag(pi.p%*%t(E.mat[,2]))}
  if(obs.seq[1]=="H"){g.0i=diag(pi.p%*%t(E.mat[,3]))}
  
  #print(g.0i)
  g.ti=matrix(0, nrow=steps, ncol = 2)
  colnames(g.ti)=c("zero","one")
  g.ti[1,]=g.0i
  
  
  for(t in 2:steps){
    
    if(obs.seq[t]=="L"){
      
      p.zero=max(diag(diag(g.ti[t-1,])%*%diag(T.ij[,1])*Gamma[1,1]))
      p.one=max(diag(diag(g.ti[t-1,])%*%diag(T.ij[,2])*Gamma[2,1]))
      
      g.ti[t,]=c(p.zero, p.one)
      
    }
    
    if(obs.seq[t]=="M"){
      
      p.zero=max(diag(diag(g.ti[t-1,])%*%diag(T.ij[,1])*Gamma[1,2]))
      p.one=max(diag(diag(g.ti[t-1,])%*%diag(T.ij[,2])*Gamma[2,2]))
      
      g.ti[t,]=c(p.zero, p.one)
      
    }
    
    if(obs.seq[t]=="H"){
      
      p.zero=max(diag(diag(g.ti[t-1,])%*%diag(T.ij[,1])*Gamma[1,3]))
      p.one=max(diag(diag(g.ti[t-1,])%*%diag(T.ij[,2])*Gamma[2,3]))
      
      g.ti[t,]=c(p.zero, p.one)
      
    }
    
    
    
  }
  
  
  for(j in 1:steps){
    hid[j]=match(max(g.ti[j,]), g.ti[j,])
  }
  
  hid=ifelse(hid==1, 0, 1)
  
  return(list(it.vals=g.ti, ML.hidden=hid))
}


ML.S=ViterbiA(obs.seq=dat$Emission, T.mat=T.ij, E.mat=Gamma, state.probs=Pi)


expect.seq=viterbi(hmm1, as.vector(as.character(dat$Emission)))
expect.seq=rbind(ML.S$ML.hidden, expect.seq)
row.names(expect.seq)=c("My Version", "Pkg: HMM")
colnames(expect.seq)=as.character(c(1:10))
print(expect.seq)

ps=ML.S$it.vals
colnames(ps)=c("P(Xt=0,Y=y | Theta)", "P(Xt=1,Y=y | Theta)")
print(ps)
#=========================================================================================


Backward.pass=function(obs.seq=NULL, T.mat=NULL, E.mat=NULL, state.probs=NULL){
  #Emission matrix (E.mat) and Transition matrix (T.mat) must have columns labeled according
  # to the respective target states.
  pi.p=state.probs
  T.mat=as.matrix(T.mat)
  E.mat=as.matrix(E.mat)
  steps=length(obs.seq)
  obs.seq=rev(obs.seq)
  
  b.ni=c(1,1)
  
  b.tj=matrix(0, nrow = steps, ncol = 2)
  b.tj[1,]=b.ni
  
  for(t in 2:steps){
    
    if(obs.seq[t-1]=="L"){
      
      b.tj[t,]=T.mat%*%diag(E.mat[,1])%*%b.tj[t-1,]
      
    }
    
    if(obs.seq[t-1]=="M"){
      
      b.tj[t,]=T.mat%*%diag(E.mat[,2])%*%b.tj[t-1,]
      
    }
    
    if(obs.seq[t-1]=="H"){
      
      b.tj[t,]=T.mat%*%diag(E.mat[,3])%*%b.tj[t-1,]
      
    }
    
    
  }
  
  
  
  return(b.tj)
  
  
}

seq.probs=Backward.pass(obs.seq=dat$Emission, T.mat=T.ij, E.mat=Gamma, state.probs=Pi)
row.names(seq.probs)=as.character(c(10:1))
colnames(seq.probs)=c("P(Y=y|Xt=0)", "P(Y=yt|Xt=1)")
print(seq.probs)
print(exp(backward(hmm1, dat$Emission)))


#===============================================================================================

#Baum-Welch Algorithm


BW.algorithm=function(obs.seq=NULL, T.mat=NULL, E.mat=NULL, state.probs=NULL){
  #initialize for space and calculations
  pi.p=state.probs
  T.mat=as.matrix(T.mat)
  E.mat=as.matrix(E.mat)
  steps=length(obs.seq)
  deltat=matrix(0, nrow = steps, ncol = 2)
  Xit.ij=list()
  #storage and while loop initialize
  transitions=list()
  emissions=list()
  states=list()
  
  #begin loop
  test.con=1
  loops=0
  while(test.con>0.0001){
    #calculate forward and backward passes
    at=Forward.pass(obs.seq, T.mat, E.mat, state.probs)
    bt=Backward.pass(obs.seq, T.mat, E.mat, state.probs)
    bt[,1]=rev(bt[,1])
    bt[,2]=rev(bt[,2])
    
    #E-step
    for(i in 1:steps){
      at.bt=diag(diag(at[i,])%*%diag(bt[i,]))
      
      #calculate Delta
      deltat[i,]=at.bt/c(at[i,]%*%bt[i,])
      
      for(t in 2:steps){
        
        #calculate Xi
        if(obs.seq[t]=="L"){
          
          combos=at[t-1,]%*%t(bt[t,])
          row1=diag(combos[1,])%*%diag(T.mat[1,])
          row2=diag(combos[2,])%*%diag(T.mat[2,])
          mat1=rbind(diag(row1), diag(row2))
          num1=mat1[,1]*E.mat[1,1]
          num2=mat1[,2]*E.mat[2,1]
          mat2=cbind(num1, num2)
          Xit.ij[[t-1]]=mat2/as.numeric(at[t-1,]%*%bt[t-1,])
          
        }
        
        if(obs.seq[t]=="M"){
          
          combos=at[t-1,]%*%t(bt[t,])
          row1=diag(combos[1,])%*%diag(T.mat[1,])
          row2=diag(combos[2,])%*%diag(T.mat[2,])
          mat1=rbind(diag(row1), diag(row2))
          num1=mat1[,1]*E.mat[1,2]
          num2=mat1[,2]*E.mat[2,2]
          mat2=cbind(num1, num2)
          Xit.ij[[t-1]]=mat2/as.numeric(at[t-1,]%*%bt[t-1,])
          
        }
        
        if(obs.seq[t]=="H"){
          
          combos=at[t-1,]%*%t(bt[t,])
          row1=diag(combos[1,])%*%diag(T.mat[1,])
          row2=diag(combos[2,])%*%diag(T.mat[2,])
          mat1=rbind(diag(row1), diag(row2))
          num1=mat1[,1]*E.mat[1,3]
          num2=mat1[,2]*E.mat[2,3]
          mat2=cbind(num1, num2)
          Xit.ij[[t-1]]=mat2/as.numeric(at[t-1,]%*%bt[t-1,])
          
        }
        
        
        
        
      }
    }#end for loops
    
    #update Parameters
    
    #transitions
    #print(Reduce("+", Xit.ij))
    T.mat=Reduce("+", Xit.ij)
    T.mat=t(T.mat)/colSums(deltat[-1,])
    #print(T.mat)
    #state probs
    Pi=deltat[1,]
    #Emissions
    total=colSums(deltat)
    gammat=cbind.data.frame(deltat,obs.seq)
    row1=colSums(as.matrix(subset(gammat, gammat[,3]=="L")[,-3]))
    row2=colSums(as.matrix(subset(gammat, gammat[,3]=="M")[,-3]))
    row3=colSums(as.matrix(subset(gammat, gammat[,3]=="H")[,-3]))
    
    #update while loop condition
    tester=cbind(row1/total, row2/total, row3/total)
    test.con=abs(tester-E.mat)
    E.mat=tester
    loops=loops+1
    # print(paste("Iteration",loops, sep = "-"))
    
    transitions[[loops]]=T.mat   
    emissions[[loops]]=E.mat
    states[[loops]]=Pi
  }
  
  calc.info=list(deltat, Xit.ij)
  
  model.parameters=list(T.mat=transitions, E.mat=emissions, state.probs=states)
  return( list(model.parameters=model.parameters, calc.info=calc.info, iterations=loops))
  
}




dat3=sim.proteinHMM(T.ij, Gamma, Pi, 50)
output=BW.algorithm(dat3$Emission, T.ij, Gamma, Pi)
T.new=output$model.parameters$T.mat[[output$iterations]]
Gam.new=output$model.parameters$E.mat[[output$iterations]]

colnames(T.new)=c(0,1)
row.names(T.new)=c(0,1)
print(T.new)
colnames(Gam.new)=c("L","M","H")
row.names(Gam.new)=c(0,1)
print(Gam.new)



#===========================================================================================================


