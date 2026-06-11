## Changelog: 
# KB 0.0.1 2026-06-11: initial 



## Function definition
Mpsr1 <- function(tseries, nEp = 2){
  # Returns the multivariate potential scale reduction (mpsr) factor assessing trend
  # Argument: 
  #   tseries
  #     is either a single or multiple time series
  #   nEp
  #     is the number of epoches
  # Value: 
  #     is Mpsr1
  
  N <- nrow(tseries); d <- ncol(tseries)
  L <- N / nEp
  index <- rep(1:nEp, each = L) #index for epoch
  
  # Epochs as separate matrices
  epochs <- vector("list", nEp)
  for(j in 1:nEp){
    epochs[[j]] <- tseries[index == j, , drop = FALSE]
  }
  
  ## means per epoch 
  means.w <- matrix(0, nEp, ncol(tseries))
  for(j in 1:nEp){
    means.w[j,] <- colMeans(epochs[[j]], na.rm = TRUE) #ignores NA's
  }
  
  ## W: average within-epoch covariance 
  W <- matrix(0, ncol(tseries), ncol(tseries))
  for(j in 1:nEp){
    W <- W + cov(epochs[[j]], use = "pairwise.complete.obs") #only considers complete pairs
  }
  W <- W / nEp
  
  # Epoch row counts (individual L's)
  epoch_L <- sapply(epochs, function(m) sum(complete.cases(m)))
  
  ## B/L: covariance of epoch means
  BL <- cov(means.w)
  
  ## Largest eigenvalue of W^{-1}(B/L)
  M <- solve(W) %*% BL
  eigvals <- eigen(M)$values
  lambda1 <- max(Re(eigvals))
  
  ## PSR^P value:
  Mpsr1 <- sqrt((median(epoch_L) - 1) / median(epoch_L) + ((nEp + 1) / nEp) * lambda1)
  
  return(Mpsr1)
}

################################################################################################