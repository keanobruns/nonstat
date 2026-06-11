## Changelog: 
# KB 0.0.1 2026-06-11: initial 



## Function definition
Mpsr2 <- function(tseries, nEp = 2, span = 3) {
  # Returns the multivariate potential scale reduction (mpsr) factor assessing variance
  # Argument: 
  #   tseries
  #     is either a single or multiple time series
  #   nEp
  #     is the number of epoches
  #   span
  #     degree of smoothing
  # Value: 
  #     is Mpsr2
  
  N <- nrow(tseries); d <- ncol(tseries)
  L <- N / nEp
  
  # Detrend each variable with LOESS  
  time <- seq_len(N)
  
  fitted_mat <- vapply(
    seq_len(ncol(tseries)),
    function(j) {
      y <- tseries[, j] #draw column
      ok <- !is.na(y)   #indicator NA
      df <- data.frame(y = y[ok], time = time[ok]) #TS without NA's
      
      fit <- loess(y ~ time, data = df, span = span) 
      pred <- predict(fit, newdata = data.frame(time = time)) #predict on only existing values
      
      pred[!ok] <- NA_real_ #put NA's back in positions
      pred
    },
    numeric(N)
  )
  
  tseries_detrended <- tseries - fitted_mat
  # ---------------------------------------------------------------
  
  # Absolute deviations from within-epoch means
  index <- rep(1:nEp, each = L)  # epoch index along rows
  # Period means per column, aligned back to rows
  period_means <- sapply(
    seq_len(d),
    function(j) ave(tseries_detrended[, j], index, FUN = function(x) mean(x, na.rm = TRUE)) #applies ave() to each column/variable; 
    #ave() applies a function (mean) to a subset of a vector (column) indexed by index
  )
  per_dev <- abs(tseries_detrended - period_means)
  
  # Epochs as separate matrices
  epochs <- vector("list", nEp)
  for(j in 1:nEp){
    epochs[[j]] <- per_dev[index == j, , drop = FALSE]
  }
  
  ## means per epoch 
  means.w <- matrix(0, nEp, ncol(per_dev))
  for(j in 1:nEp){
    means.w[j,] <- colMeans(epochs[[j]], na.rm = TRUE) #ignores NA's
  }
  
  ## W: average within-epoch covariance 
  W <- matrix(0, ncol(per_dev), ncol(per_dev))
  for(j in 1:nEp){
    W <- W + cov(epochs[[j]], use = "pairwise.complete.obs") #only considers complete pairs
  }
  W <- W / nEp
  
  # Epoch row counts (individual L's)
  epoch_L <- sapply(epochs, function(m) sum(complete.cases(m)))
  
  # B/L: Covariance of epoch means 
  BL <- cov(means.w)
  
  # Largest eigenvalue of W^{-1}(B/L)
  M <- solve(W) %*% BL
  eigvals <- eigen(M)$values
  lambda1 <- max(Re(eigvals))
  
  # PSR^P value
  Mpsr2 <- sqrt((median(epoch_L) - 1) / median(epoch_L) + ((nEp + 1) / nEp) * lambda1)
  return(Mpsr2)
}

########################################################################################################