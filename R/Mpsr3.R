## Changelog: 
# KB 0.0.1 2026-06-11: initial 



## Function definition
Mpsr3 <- function(tseries, nEp = 2, span = 3, lag = 1) {
  # Returns the multivariate potential scale reduction (mpsr) factor assessing auto-correlation
  # Argument: 
  #   tseries
  #     is either a single or multiple time series
  #   nEp
  #     is the number of epoches
  #   span 
  #     degree of smoothing
  #   lag
  #     is specified lag
  # Value: 
  #     is Mpsr3
  
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
  # --------------------------------------------------------------
  
  #Standardization
  index <- rep(1:nEp, each = L)
  
  period_means <- sapply(
    seq_len(d),
    function(j) ave(tseries_detrended[, j], index, FUN = function(x) mean(x, na.rm = TRUE))
  )
  
  s_j <- sapply(
    seq_len(d),
    function(l) ave(tseries_detrended[, l], index, FUN = function(x) {
      ok <- !is.na(x)
      m <- mean(x[ok])
      sqrt(sum((x[ok] - m)^2) / sum(ok))
    })
  )
  
  per_dev <- (tseries_detrended - period_means) / s_j
  # ---------------------------------------------------------------
  
  
  #Lag-1 products
  prev <- sapply(seq_len(ncol(per_dev)), function(j)
    ave(per_dev[, j], index, FUN = function(x) c(rep(NA_real_, lag), head(x, -lag)))
  ) #creates a lagged version of the matrix with first row = NA
  
  per_dev <- per_dev * prev #multiplies the matrices
  per_dev[is.na(per_dev)] = NA_real_ #make sure all NA's are missing data
  #-----------------------------------------------------------------
  
  
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
  Mpsr3 <- sqrt((median(epoch_L) - 1) / median(epoch_L) + ((nEp + 1) / nEp) * lambda1)
  return(Mpsr3)
}

########################################################################################################