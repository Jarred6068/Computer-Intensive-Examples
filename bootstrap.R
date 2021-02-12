
boots=function(x, resamples=NULL, is.matrix=FALSE, alpha.level=0.05){
  #SYNTAX:
  #========================================================
  #x--can be a vector or dataframe of data for resampling
  #resamples--the number of bootstrap resamplings desired
  #is.matrix--logical indicating if data in x is a matrix 
  #or not--default:FALSE
  #alpha.level--confidence level of interval estimate
  #========================================================
  #space allocation
  idx=NULL
  indicies=NULL
  rowind=NULL
  colind=NULL
  sub.samp=NULL
  upper.limit=NULL
  lower.limit=NULL
  est.mean=NULL
  SE=NULL
  mean.vec=NULL
  conf.interval=NULL
#-------------bootstrap of matrix---------
  if(is.matrix==TRUE){
    
    #get dimensions of matrix
    idx=dim(x)
    #vectorize row indicies
    rowind=c(1:idx[1])
    
    #additional matrix-specific allocation
    mean.vec=as.data.frame(matrix(rep(0, idx[2]*resamples), ncol = idx[2], nrow = resamples))
    conf.interval=as.data.frame(matrix(rep(0,idx[2]*3), ncol = 3, nrow = idx[2]))
    
    
    #resample dataset
    for(i in 1:resamples){
      
      #random sample w/replacement of row indicies
      indicies=sample(rowind, idx[1], replace = TRUE)
      #index matrix x for indicies
      sub.samp=x[indicies, ]
      #obtain column means of resampled matrix
      mean.vec[i, ]=colMeans(sub.samp)
        
    }
    
 
    #confidence interval construction by column
    for(j in 1:idx[2]){
      
      #calculate standard error of means vector
      SE[j]=sd(mean.vec[,j])
      #compute 100*(1-alpha)% conf. interval
      est.mean[j]=mean(x[, j])
      upper.limit=est.mean[j]+qnorm(1-alpha.level/2, 0, 1)*SE[j]
      lower.limit=est.mean[j]+qnorm(alpha.level/2, 0, 1)*SE[j]
      conf.interval[j, ]=c(lower.limit, est.mean[j], upper.limit)
      
      
    }
    
    #--naming--
    names(SE)=colnames(x)
    colnames(conf.interval)=c("lower.limit", "estiamted.mean","upper.limit")
    row.names(conf.interval)=colnames(x)
    
    
#------------bootstrap of vector----------
  }else{

    #get length of vector
    idx=length(x)
    #vectorize indicies
    rowind=c(1:idx)
    #additional allocation
    mean.vec=NULL
    conf.interval=NULL
    
    for(k in 1:resamples){
      
      #random sample w/replacement of indicies
      indicies=sample(rowind, idx, replace = TRUE)
      #index vector x for indicies
      sub.samp=x[indicies]
      #obtain mean of resampled vector
      mean.vec[k]=mean(sub.samp)
      
    }
    
    
    #calculate standard error of means vector
    SE=sd(mean.vec)
    #compute 100*(1-alpha)% conf. interval
    est.mean=mean(x)
    upper.limit=est.mean+qnorm(1-alpha.level/2, 0, 1)*SE
    lower.limit=est.mean+qnorm(alpha.level/2, 0, 1)*SE
    conf.interval=c(lower.limit,est.mean,upper.limit)
    names(conf.interval)=c("lower.limit", "estiamted.mean","upper.limit")
    
    
  }
  
  
  
  return(list(c("-----Confidence Interval-----"), conf.interval, c("----Estimated Standard Error----"), SE))
  
  
}





