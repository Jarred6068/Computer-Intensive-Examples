bsclass=function(data, classifier="lda", LOOCV=FALSE ,resamples=NULL, CVfolds=NULL, alpha.level=0.05, k=NULL,
                 scale.it=TRUE, kernel.type='linear', degrees=3, gamma.p=NULL, coef0.p=0, cost.p=1){
  #==========================================+SYNTAX+===========================================
  #____________________________________GENERAL-PARAMETERS_______________________________________
  
  #data -- a data frame with the first column containing the classes
  #Classifier -- a string specifying one of 'lda' , 'qda' , 'svm' , 'knn' 
  #LOOCV -- logical--if true conducts LOOCV in class of mfold (CVfolds then negligable paramter)
  #resamples -- the number of bootstrap resamples to take on the estimates
  #CVfolds -- the parameter k for k-fold cross validation
  #alpha.level -- the confidence level for the approximate CI around the est.
  
  #___________________________________FOR-K-NEAREST-NEIGHBOR____________________________________
  
  #k -- if using knn; it is the paramter for number of neighbors to consider
  
  #__________________________________FOR-SUPPORT-VECTOR-MACHINE_________________________________
  
  #scale.it -- logical indicating if the data should be scaled prior to training
  #kernel.type -- type of kernel function to be used: 'linear', 'polynomial, 'radial basis', 'sigmoid'
  #degrees -- parameter needed for kernel of type 'polynomial' (default: 3)
  #gamma.p -- parameter needed for all kernels except 'linear' (default: 1/(data dimension))
  #coef0.p -- parameter needed for kernels of type 'polynomial' and 'sigmoid' (default: 0)
  #cost.p -- cost of constraints violation (default: 1): general svm parameter
  
  #=============================================================================================
  library(MASS)
  library(caret)
  library(class)
  library(e1071)
  
  idx=dim(data[,-1])
  rowind=c(1:idx[1])
  
  
  
  
  #==========================determine=if=LOOCV===============================
  if(LOOCV==FALSE){
    #---space--allocation---
    indices=NULL
    trainer=NULL
    tester=NULL
    group.train=NULL
    group.test=NULL
    apparent.error=rep(0,CVfolds)
    actual.error=rep(0,CVfolds)
    timespec=rep(0,CVfolds)
    SE=NULL
    mean.error=rep(0,resamples)
    upper.limit=NULL
    lower.limit=NULL
    indexes=NULL
    #-----------------------
    #partition the data
    indices=createFolds(rowind, CVfolds, list = TRUE, returnTrain = FALSE)
    #-----------------------
    
    #--------------------------LDA----------------------------------- 
    
    if(classifier=="lda"){
      
      
      
      for( i in 1:CVfolds){
        start.time=Sys.time()
        #select training/testing sets
        trainer=data[-indices[[i]], ]             
        tester=data[indices[[i]], ]
        #extract set labels
        group.train=as.factor(trainer[,1])
        group.test=as.factor(tester[,1])
        #model training
        lda.fit=lda(trainer[,-1], group.train)
        #predict training set
        lda.train=predict(lda.fit, trainer[,-1])
        preds=table(actual=group.train, predicted = lda.train$class)
        #apparent error rate
        apparent.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        #calibration
        #predict test set
        lda.calib=predict(lda.fit, tester[,-1])
        preds=table(actual=group.test, predicted = lda.calib$class)
        #calculate actual error rate
        actual.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        end.time=Sys.time()
        timespec[i]=end.time-start.time
      }
      
      
      #bootstrap the estimates to get a mean and SE
      
      for ( j in 1:resamples){
        
        indexes=sample(c(1:length(actual.error)), CVfolds, replace = TRUE)
        info=actual.error[indexes]
        mean.error[j]=mean(info)
        
        
      }
      
      #calculate sample statistics and produce CI around Error
      SE=sd(mean.error)
      names(SE)=c("Est. SE")
      upper.limit=mean(actual.error)+qnorm(1-alpha.level/2, 0, 1)*SE
      lower.limit=mean(actual.error)+qnorm(alpha.level/2, 0, 1)*SE
      conf.interval=c(upper.limit, mean(actual.error), lower.limit)
      names(conf.interval)=c("upper.limit", "mean", "lower.limit")
      
    }
    
    
    #--------------------------QDA----------------------------------- 
    
    if(classifier=="qda"){
      
      
      
      for( i in 1:CVfolds){
        start.time=Sys.time()
        #select training/testing sets
        trainer=data[-indices[[i]], ]
        tester=data[indices[[i]], ]
        #extract set labels
        group.train=as.factor(trainer[,1])
        group.test=as.factor(tester[,1])
        #model training
        qda.fit=qda(trainer[,-1], group.train)
        #predict training set
        qda.train=predict(qda.fit, trainer[,-1])
        preds=table(actual=group.train, predicted = qda.train$class)
        #apparent error rate
        apparent.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        #calibration
        #predict test set
        qda.calib=predict(qda.fit, tester[,-1])
        preds=table(actual=group.test, predicted = qda.calib$class)
        #calculate actual error rate
        actual.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        end.time=Sys.time()
        timespec[i]=end.time-start.time
      }
      
      
      
      #bootstrap the estimates to get a mean and SE
      
      for ( j in 1:resamples){
        
        indexes=sample(c(1:length(actual.error)), CVfolds, replace = TRUE)
        info=actual.error[indexes]
        mean.error[j]=mean(info)
        
      }
      
      #calculate sample statistics and produce CI around Error
      SE=sd(mean.error)
      names(SE)=c("Est. SE")
      upper.limit=mean(actual.error)+qnorm(1-alpha.level/2, 0, 1)*SE
      lower.limit=mean(actual.error)+qnorm(alpha.level/2, 0, 1)*SE
      conf.interval=c(upper.limit, mean(actual.error), lower.limit)
      names(conf.interval)=c("upper.limit", "mean", "lower.limit")
      
      
    }
    
    
    #--------------------------SVM----------------------------------- 
    
    
    if(classifier=="svm"){
      
      
      
      
      for( i in 1:CVfolds){
        start.time=Sys.time()
        #select training/testing sets
        trainer=data[-indices[[i]], ]
        tester=data[indices[[i]], ]
        #extract set labels
        group.train=as.factor(trainer[,1])
        group.test=as.factor(tester[,1])
        #model training
        gamma.p=ifelse(is.null(gamma.p)==TRUE, 1/ncol(trainer), gamma.p )
        svm.fit=svm(trainer[,-1], group.train, scale = scale.it, kernel = kernel.type, 
                    degree = degrees, gamma = gamma.p, coef0 = coef0.p, cost = cost.p)
        #predict training set
        preds=table(actual=group.train, predicted = svm.fit$fitted)
        #apparent error rate
        apparent.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        #calibration
        #predict test set
        svm.calib=predict(svm.fit, tester[,-1])
        preds=table(actual=group.test, predicted = svm.calib)
        #calculate actual error rate
        actual.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        end.time=Sys.time()
        timespec[i]=end.time-start.time
      }
      
      
      
      
      #bootstrap the estimates to get a mean and SE
      
      for ( j in 1:resamples){
        
        indexes=sample(c(1:length(actual.error)), CVfolds, replace = TRUE)
        info=actual.error[indexes]
        mean.error[j]=mean(info)
        
      }
      
      #calculate sample statistics and produce CI around Error
      SE=sd(mean.error)
      names(SE)=c("Est. SE")
      upper.limit=mean(actual.error)+qnorm(1-alpha.level/2, 0, 1)*SE
      lower.limit=mean(actual.error)+qnorm(alpha.level/2, 0, 1)*SE
      conf.interval=c(upper.limit, mean(actual.error), lower.limit)
      names(conf.interval)=c("upper.limit", "mean", "lower.limit")
      
    }
    
    
    #--------------------------KNN----------------------------------- 
    
    if(classifier=="knn"){
      
      
      
      for( i in 1:CVfolds){
        start.time=Sys.time()
        #select training/testing sets
        trainer=data[-indices[[i]], ]
        tester=data[indices[[i]], ]
        #extract set labels
        group.train=as.factor(trainer[,1])
        group.test=as.factor(tester[,1])
        #repredict training set
        knn.fit=knn(trainer[,-1], trainer[,-1], group.train, k)
        preds=table(actual=group.train, predicted = knn.fit)
        #apparent error rate
        apparent.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        #calibration
        #predict test set
        knn.fit=knn(trainer[,-1], tester[,-1], group.train, k)
        preds=table(actual=group.test, predicted = knn.fit)
        #calculate actual error rate
        actual.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        end.time=Sys.time()
        timespec[i]=end.time-start.time
      }
      
      
      
      #bootstrap the estimates to get a mean and SE
      
      for ( j in 1:resamples){
        
        indexes=sample(c(1:length(actual.error)), CVfolds, replace = TRUE)
        info=actual.error[indexes]
        mean.error[j]=mean(info)
        
      }
      
      #calculate sample statistics and produce CI around Error
      SE=sd(mean.error)
      names(SE)=c("Est. SE")
      upper.limit=mean(actual.error)+qnorm(1-alpha.level/2, 0, 1)*SE
      lower.limit=mean(actual.error)+qnorm(alpha.level/2, 0, 1)*SE
      conf.interval=c(upper.limit, mean(actual.error), lower.limit)
      names(conf.interval)=c("upper.limit", "mean", "lower.limit")
      
      
    }
    
    
    return(list(output.info=cbind.data.frame(apparent.error, actual.error, timespec),
                c("---bootstrapped--error--rates---"),mean.error, SE, conf.interval))
  } 
  #========================conduct=leav=one=out=cross=validation===================================
  if(LOOCV==TRUE){
    #---space--allocation---
    indices=NULL
    trainer=NULL
    tester=NULL
    group.train=NULL
    group.test=NULL
    apparent.error=rep(0,idx[1])
    actual.error=rep(0,idx[1])
    timespec=rep(0,idx[1])
    SE=NULL
    mean.error=rep(0,resamples)
    upper.limit=NULL
    lower.limit=NULL
    indexes=NULL
    #-----------------------
    #partition the data
    indices=createFolds(rowind, idx[1], list = TRUE, returnTrain = FALSE)
    #-----------------------
    
    
    if(classifier=="lda"){
      
      
      
      for( i in 1:idx[1]){
        start.time=Sys.time()
        #select training/testing sets
        trainer=data[-indices[[i]], ]             
        tester=data[indices[[i]], ]
        #extract set labels
        group.train=as.factor(trainer[,1])
        group.test=as.factor(tester[,1])
        #model training
        lda.fit=lda(trainer[,-1], group.train)
        #predict training set
        lda.train=predict(lda.fit, trainer[,-1])
        preds=table(actual=group.train, predicted = lda.train$class)
        #apparent error rate
        apparent.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        #calibration
        #predict test observation
        lda.calib=predict(lda.fit, tester[,-1])
        #calculate actual error rate
        actual.error[i]=ifelse(lda.calib$class==group.test, 1, 0)
        end.time=Sys.time()
        timespec[i]=end.time-start.time
      }
      
      
      #bootstrap the estimates to get a mean and SE
      
      for ( j in 1:resamples){
        
        indexes=sample(c(1:length(actual.error)), idx[1], replace = TRUE)
        info=actual.error[indexes]
        mean.error[j]=1-mean(info)
        
        
      }
      
      #calculate sample statistics and produce CI around Error
      SE=sd(mean.error)
      names(SE)=c("Est. SE")
      upper.limit=(1-mean(actual.error))+qnorm(1-alpha.level/2, 0, 1)*SE
      lower.limit=(1-mean(actual.error))+qnorm(alpha.level/2, 0, 1)*SE
      conf.interval=c(upper.limit, 1-mean(actual.error), lower.limit)
      names(conf.interval)=c("upper.limit", "mean", "lower.limit")
      
    }
    
    
    #--------------------------QDA----------------------------------- 
    
    if(classifier=="qda"){
      
      
      
      for( i in 1:idx[1]){
        start.time=Sys.time()
        #select training/testing sets
        trainer=data[-indices[[i]], ]             
        tester=data[indices[[i]], ]
        #extract set labels
        group.train=as.factor(trainer[,1])
        group.test=as.factor(tester[,1])
        #model training
        qda.fit=qda(trainer[,-1], group.train)
        #predict test set
        qda.train=predict(qda.fit, trainer[,-1])
        preds=table(actual=group.train, predicted = qda.train$class)
        #apparent error rate
        apparent.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        #calibration
        #predict test observation
        qda.calib=predict(qda.fit, tester[,-1])
        #calculate actual error rate
        actual.error[i]=ifelse(qda.calib$class==group.test, 1, 0)
        end.time=Sys.time()
        timespec[i]=end.time-start.time
      }
      
      
      
      #bootstrap the estimates to get a mean and SE
      
      for ( j in 1:resamples){
        
        indexes=sample(c(1:length(actual.error)), idx[1], replace = TRUE)
        info=actual.error[indexes]
        mean.error[j]=1-mean(info)
        
      }
      
      #calculate sample statistics and produce CI around Error
      SE=sd(mean.error)
      names(SE)=c("Est. SE")
      upper.limit=(1-mean(actual.error))+qnorm(1-alpha.level/2, 0, 1)*SE
      lower.limit=(1-mean(actual.error))+qnorm(alpha.level/2, 0, 1)*SE
      conf.interval=c(upper.limit, 1-mean(actual.error), lower.limit)
      names(conf.interval)=c("upper.limit", "mean", "lower.limit")
      
      
    }
    
    
    #--------------------------SVM----------------------------------- 
    
    
    if(classifier=="svm"){
      
      
      
      for( i in 1:idx[1]){
        start.time=Sys.time()
        #select training/testing sets
        trainer=data[-indices[[i]], ]             
        tester=data[indices[[i]], ]
        #extract set labels
        group.train=as.factor(trainer[,1])
        group.test=as.factor(tester[,1])
        #model training
        gamma.p=ifelse(is.null(gamma.p)==TRUE, 1/ncol(trainer), gamma.p )
        svm.fit=svm(trainer[,-1], group.train, scale = scale.it, kernel = kernel.type, 
                    degree = degrees, gamma = gamma.p, coef0 = coef0.p, cost = cost.p)
        #predict test set
        preds=table(actual=group.train, predicted = svm.fit$fitted)
        #apparent error rate
        apparent.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        #calibration
        svm.calib=predict(svm.fit, tester[,-1])
        #calculate actual error rate
        actual.error[i]=ifelse(svm.calib==group.test, 1, 0)
        end.time=Sys.time()
        timespec[i]=end.time-start.time
      }
      
      
      
      
      #bootstrap the estimates to get a mean and SE
      
      for ( j in 1:resamples){
        
        indexes=sample(c(1:length(actual.error)), idx[1], replace = TRUE)
        info=actual.error[indexes]
        mean.error[j]=1-mean(info)
        
      }
      
      #calculate sample statistics and produce CI around Error
      SE=sd(mean.error)
      names(SE)=c("Est. SE")
      upper.limit=(1-mean(actual.error))+qnorm(1-alpha.level/2, 0, 1)*SE
      lower.limit=(1-mean(actual.error))+qnorm(alpha.level/2, 0, 1)*SE
      conf.interval=c(upper.limit, 1-mean(actual.error), lower.limit)
      names(conf.interval)=c("upper.limit", "mean", "lower.limit")
      
      
      
    }
    
    
    #--------------------------KNN----------------------------------- 
    
    if(classifier=="knn"){
      
      
      
      for( i in 1:idx[1]){
        start.time=Sys.time()
        #select training/testing sets
        trainer=data[-indices[[i]], ]             
        tester=data[indices[[i]], ]
        #extract set labels
        group.train=as.factor(trainer[,1])
        group.test=as.factor(tester[,1])
        knn.fit=knn(trainer[,-1], trainer[,-1], group.train, k)
        preds=table(actual=group.train, predicted = knn.fit)
        #apparent error rate
        apparent.error[i]=1-sum(diag(preds))/sum(colSums(preds))
        #calibration
        #predict test set
        knn.fit=knn(trainer[,-1], tester[,-1], group.train, k)
        preds=table(actual=group.test, predicted = knn.fit)
        #calculate actual error rate
        actual.error[i]=ifelse(knn.fit==group.test, 1, 0)
        end.time=Sys.time()
        timespec[i]=end.time-start.time
      }
      
      
      
      #bootstrap the estimates to get a mean and SE
      
      for ( j in 1:resamples){
        
        indexes=sample(c(1:length(actual.error)), idx[1], replace = TRUE)
        info=actual.error[indexes]
        mean.error[j]=1-mean(info)
        
      }
      
      #calculate sample statistics and produce CI around Error
      SE=sd(mean.error)
      names(SE)=c("Est. SE")
      upper.limit=(1-mean(actual.error))+qnorm(1-alpha.level/2, 0, 1)*SE
      lower.limit=(1-mean(actual.error))+qnorm(alpha.level/2, 0, 1)*SE
      conf.interval=c(upper.limit, 1-mean(actual.error), lower.limit)
      names(conf.interval)=c("upper.limit", "mean", "lower.limit")
      
      
    }
    
    
    return(list(output.info=cbind.data.frame(apparent.error, actual.error, timespec),
                c("---bootstrapped--error--rates---"),mean.error, SE, conf.interval))
    
    
  }
  
  
  
  
}
#-------------------------------------end-function----------------------------------------------
