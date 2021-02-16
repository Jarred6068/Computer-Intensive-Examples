# Function to generate p-values under the null hypothesis
#
# Input
# N: the number of p-values to generate
# n: the sample size for each hypothesis test
# lambda_max: the max value for lambda
# plot_hist: logical argument for plotting the histogram
generatePValues <- function (N, n, lambda_max, plot_hist = TRUE) {
    # Initialize
    pvalues <- NULL
    lambdas <- NULL
    data <- matrix (NA, nrow=n, ncol=N)
    z_stats <- NULL
    
    # Iterations
    for (i in 1:N) {
        # In each iteration
        # Simulate lambda
        lambdas[i] <- runif (1, min=0, max=lambda_max)
        
        # Simulate n observations from exp(lambda)
        # Can call your own exp simulator here!
        data[,i] <- rexp (n, rate=lambdas[i])
        
        # Calculate sample mean and sample sd
        sample_mean <- mean (data[,i])
        sample_sd <- sd (data[,i])
        
        # Calculate z statistic
        z_stats[i] <- (sample_mean - 1/lambdas[i])/(sample_sd/sqrt(n))
        
        # Calculate p-value
        pvalues[i] <- pnorm (z_stats[i], mean=0, sd=1, lower.tail=FALSE)
    }
    
    if (plot_hist == TRUE) {
        # Generate a histogram of the p-values
        hist (pvalues, freq=FALSE)
    }
    
    return (list(lambdas=lambdas, data=data, z_stats=z_stats, pvalues=pvalues))
}

# test the function
simu <- generatePValues (N=100, n=100, lambda_max=10, plot_hist = TRUE)

# longer runs
start_time <- Sys.time()
simu <- generatePValues (N=1e4, n=100, lambda_max=10, plot_hist = TRUE)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
simu <- generatePValues (N=1e4, n=1e3, lambda_max=10, plot_hist = TRUE)
end_time <- Sys.time()
end_time - start_time
