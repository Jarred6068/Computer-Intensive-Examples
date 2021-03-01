#NOTICE of INTELLECTUAL PROPERTY: The following code was shared with the 
#STAT 565 class of 2020 as an example of the MCMC algorithm. The original Developer and owner of 
#the following code is Dr. Audrey Fu of the University of Idaho - Department of Mathematics and Statistical Science

#############################################################
# Implement an MH algorithm to sample from Beta(alpha, beta)
# using uniforms as the proposal distribution
#############################################################
mhBetaUnif <- function (alpha, beta, initial=0.5, N=1000, verbose=TRUE) {
    # Initialization
    theta <- rep (0, N)
    theta[1] <- initial
    
    # MH iteration
    for (i in 2:N) {
        # Generate a proposal
        theta.new <- runif (1, min=0, max=1)
        if (verbose) {
            print (paste ("i=", i, " current: ", theta[i-1], "; propose:", theta.new))
        }

        # Calculate the log-transformed acceptance probability
        # Calculate each component separately
        loglik.new <- calcLogLikBeta (theta=theta.new, alpha=alpha, beta=beta)
        log.prop.num <- dunif (theta[i], min=0, max=1, log=TRUE)
        loglik.curr <- calcLogLikBeta (theta=theta[i-1], alpha=alpha, beta=beta)
        log.prop.denom <- dunif (theta.new, min=0, max=1, log=TRUE)
        
        log.ratio <- loglik.new + log.prop.num - loglik.curr - log.prop.denom
        log.accept.prob <- min (log.ratio, 0)
        
        # Decide whether to accept the proposal
        log.u <- log (runif (1))
        if (log.u < log.accept.prob) {
            theta[i] <- theta.new
        } else {
            theta[i] <- theta[i-1]
        }
    }
    
    return (theta)
    
}

# Function to calculate the log likelihood
# of a beta random variable
calcLogLikBeta <- function (theta, alpha, beta){
    loglik <- - log(beta(alpha, beta)) + (alpha-1)*log(theta) + (beta-1)*log(1-theta)
    return (loglik)
}

# Function to post-process an MCMC sample
mhPostProcess <- function (post.sample, burn.in = 0.2, step.size=10) {
    # calculate the length of the MCMC sample
    N <- length (post.sample)
    # create the new index sequence
    seq.new <- 0:(N - floor (N*burn.in))
    # find the indices at an interval defined by step.size
    seq.choose <- which ((seq.new %% step.size) == 0)
    # return realizations at chosen indices
    return (post.sample[floor (N*burn.in)-1+seq.choose])
}


# test runs
# short run
beta.sample <- mhBetaUnif (alpha=7, beta=7)
beta.sample.post <- mhPostProcess (beta.sample)

# long run
beta.sample <- mhBetaUnif (N=1e6, alpha=7, beta=7, verbose=FALSE)
beta.sample.post <- mhPostProcess (beta.sample, step.size=200)

# diagnostic plots
# trace plot
plot (1:length (beta.sample), beta.sample, type = 'b', pch=16, cex=0.4, col="blue", xlab="index", ylab="MCMC sample", ylim=c(0,1))
plot (1:length (beta.sample.post), beta.sample.post, type = 'b', pch=16, cex=0.4, col="blue", xlab="index", ylab="MCMC sample", ylim=c(0,1))
# autocorrelation plot
acf (beta.sample)
acf (beta.sample.post)
# histogram
hist (beta.sample.post, freq=FALSE, col="yellow", xlab="MCMC sample", main="")


#############################################################
# Implement an MH algorithm to sample from Beta(alpha, beta)
# using normals as the proposal distribution
#############################################################
mhBetaNorm <- function (alpha, beta, initial=0.5, N=1000, prop.sd=0.05, verbose=TRUE) {
    # Initialization
    theta <- rep (0, N)
    theta[1] <- initial
    
    # MH iteration
    for (i in 2:N) {
        # Generate a proposal
        theta.new <- rnorm (1, mean=theta[i-1], sd=prop.sd)
        if (verbose) {
            print (paste ("i=", i, " current: ", theta[i-1], "; propose:", theta.new))
        }

        # Calculate the log-transformed acceptance probability
        # Calculate each component separately
        loglik.new <- calcLogLikBeta (theta=theta.new, alpha=alpha, beta=beta)
        log.prop.num <- dnorm (theta[i-1], mean=theta.new, sd=prop.sd, log=TRUE)
        loglik.curr <- calcLogLikBeta (theta=theta[i-1], alpha=alpha, beta=beta)
        log.prop.denom <- dnorm (theta.new, mean=theta[i-1], sd=prop.sd, log=TRUE)
        
        log.ratio <- loglik.new + log.prop.num - loglik.curr - log.prop.denom
        log.accept.prob <- min (log.ratio, 0)
        
        # Decide whether to accept the proposal
        log.u <- log (runif (1))
        if (log.u < log.accept.prob) {
            theta[i] <- theta.new
        } else {
            theta[i] <- theta[i-1]
        }
    }
    
    return (theta)
    
}

# test runs
# short run
beta.sample.norm <- mhBetaNorm (alpha=7, beta=7, prop.sd=0.02)
beta.sample.post.norm <- mhPostProcess (beta.sample.norm)

# long run
beta.sample.norm <- mhBetaNorm (N=1e6, alpha=7, beta=7, prop.sd=0.02, verbose=FALSE)
beta.sample.post.norm <- mhPostProcess (beta.sample.norm, step.size=200)

# diagnostic plots
# trace plot
plot (1:length (beta.sample.post.norm), beta.sample.post.norm, type = 'b', pch=16, cex=0.4, col="blue", xlab="index", ylab="MCMC sample", ylim=c(0,1))

# autocorrelation plot
acf (beta.sample.norm)
acf (beta.sample.post.norm)

# histogram
hist (beta.sample.post.norm, freq=FALSE, col="yellow", xlab="MCMC sample", main="")

###############################
# explore different proposals
###############################
beta.sample.norm <- mhBetaNorm (alpha=7, beta=7, prop.sd=0.005)
plot (1:length (beta.sample.norm), beta.sample.norm, type = 'b', pch=16, cex=0.4, col="blue", xlab="index", ylab="MCMC sample", ylim=c(0,1))

beta.sample.norm <- mhBetaNorm (alpha=7, beta=7, prop.sd=0.01)
plot (1:length (beta.sample.norm), beta.sample.norm, type = 'b', pch=16, cex=0.4, col="blue", xlab="index", ylab="MCMC sample", ylim=c(0,1))

beta.sample.norm <- mhBetaNorm (alpha=7, beta=7, prop.sd=0.05)
plot (1:length (beta.sample.norm), beta.sample.norm, type = 'b', pch=16, cex=0.4, col="blue", xlab="index", ylab="MCMC sample", ylim=c(0,1))

beta.sample.norm <- mhBetaNorm (alpha=7, beta=7, prop.sd=0.1)
plot (1:length (beta.sample.norm), beta.sample.norm, type = 'b', pch=16, cex=0.4, col="blue", xlab="index", ylab="MCMC sample", ylim=c(0,1))
