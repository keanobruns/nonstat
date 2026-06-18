## Changelog: 
# KB 0.0.1 2026-06-11: initial



## Documentation: 
#' @title Test for (multivariate) nonstationarity
#' @description Applies a nonvisual, diagnostic-based screening procedure to determine whether a time series (univariate or multivariate) violates the assumption of stationarity. Specifically, in the univariate case, the function evaluates (a) the presence of a trend, (b) changes in variance over time, and (c) changes in the auto-correlation over time. In the multivariate case contemporaneous and lagged correlation are also evaluated. These dimensions of nonstationarity are assessed using R-hat-type statistics adapted from Bayesian convergence diagnostics and Levene's test.
#' @param tseries a numerical vector, matrix or data frame
#' @param nEp the number of epochs (in which time series is cut for PSR calculation)
#' @param cut.mpsr1 threshold for the trend diagnostic, which assesses trending over time
#' @param cut.mpsr2 threshold for the changing variance diagnostic, which assesses whether the variances are changing over time
#' @param cut.mpsr3 threshold for the changing auto-correlation diagnostic, which assesses whether the auto-correlations are changing over time
#' @param cut.mpsr4 threshold for the changing correlation diagnostic, which assesses whether the contemporaneous correlation is changing over time
#' @param cut.mpsr5 threshold for the changing lagged correlation diagnostic, which assesses whether the lagged correlation is changing over time
#' @param span numerical value that is passed to the \code{loess} function
#' @param lag specification of the time lag for the auto-correlation and lagged correlation diagnostics
#' @return a logical scalar indicating whether the process has been diagnosed as non-stationary (\code{TRUE}) or stationary (\code{FALSE})
#' @example 
#' set.seed(356479)
#' x <- rnorm(100)
#' is.mult.nonstat(x)


## Function definition
is.mult.nonstat <- function(tseries, nEp = 2, cut.mpsr1 = 1.2, cut.mpsr2 = 1.1, cut.mpsr3 = 1.01, 
                            cut.mpsr4 = 1.01, cut.mpsr5 = 1.01, span = 3, lag = 1 ){
  
  
  # Convert data.frame to matrix
  if (is.data.frame(tseries)) {
    tseries <- as.matrix(tseries)
  }
  
  # check if tseries is univariate and convert accordingly
  if (is.null(dim(tseries))) {
    tseries <- matrix(tseries, ncol = 1)
    univariate <- TRUE
  } else {
    univariate <- ncol(tseries) == 1
  }
  
  if (univariate) {
    message("Input data is univariate: Mpsr4 and Mpsr5 cannot be calculated.")
  }
  
  
  # loop to remove elements until the length is both divisible by nEp and greater than nEp
  while (nrow(tseries) %% nEp != 0 && nrow(tseries) > 3*nEp) {
    tseries <- tseries[-nrow(tseries), ]  # Remove last element
  }
  
  if (nrow(tseries) %% nEp == 0 && nrow(tseries) >= 3*nEp ){
    
    Mpsr1 <- try( Mpsr1( tseries=tseries, nEp=nEp ) )
    Mpsr2 <- try( Mpsr2( tseries=tseries, nEp=nEp, span=span ) )
    Mpsr3 <- try( Mpsr3( tseries=tseries, nEp=nEp, span=span , lag=lag) )
    
    if (!univariate) {
      Mpsr4 <- try( Mpsr4( tseries=tseries, nEp=nEp, span=span ) )
      Mpsr5 <- try( Mpsr5( tseries=tseries, nEp=nEp, span=span , lag=lag) )
    } else {
      Mpsr4 <- NULL
      Mpsr5 <- NULL
    }
    
    if( inherits( Mpsr1, "try-error" ) ) Mpsr1 <- NULL
    if( inherits( Mpsr2, "try-error" ) ) Mpsr2 <- NULL
    if( inherits( Mpsr3, "try-error" ) ) Mpsr3 <- NULL
    if( inherits( Mpsr4, "try-error" ) ) Mpsr4 <- NULL
    if( inherits( Mpsr5, "try-error" ) ) Mpsr5 <- NULL
    
    
    if( !is.null( Mpsr1 ) && !is.na( Mpsr1 ) && !is.infinite( Mpsr1 ) ) mult.mean <- Mpsr1 > cut.Mpsr1 else {
      mult.mean <- NULL
      warning( "Unable to assess the multivariate nonstationarity of the mean." )
    }
    if( !is.null( Mpsr2 ) && !is.na( Mpsr2 ) && !is.infinite( Mpsr2 ) ) mult.var  <- Mpsr2 > cut.Mpsr2 else {
      mult.var  <- NULL
      warning( "Unable to assess the multivariate nonstationarity of the variance." )
    }
    if( !is.null( Mpsr3 ) && !is.na( Mpsr3 ) && !is.infinite( Mpsr3 ) ) mult.ar  <- Mpsr3 > cut.Mpsr3 else {
      mult.ar  <- NULL
      warning( "Unable to assess the multivariate nonstationarity of the autocorrelation." )
    }
    
    if (!univariate) {
      if( !is.null( Mpsr4 ) && !is.na( Mpsr4 ) && !is.infinite( Mpsr4 ) ) mult.cor  <- Mpsr4 > cut.Mpsr4 else {
        mult.cor  <- NULL
        warning( "Unable to assess the multivariate nonstationarity of the correlation." )
      }
      if( !is.null( Mpsr5 ) && !is.na( Mpsr5 ) && !is.infinite( Mpsr5 ) ) mult.lagcor  <- Mpsr5 > cut.Mpsr5 else {
        mult.lagcor  <- NULL
        warning( "Unable to assess the multivariate nonstationarity of the lagged correlation." )
      }
    } else {
      mult.cor <- NULL
      mult.lagcor <- NULL
    }
    
    if( !is.null( mult.mean ) && !is.null( mult.var ) && !is.null( mult.ar ) &&
        (univariate || (!is.null(mult.cor) && !is.null(mult.lagcor))) ){
      
      if (univariate) {
        ret <- mult.mean | mult.var | mult.ar
      } else {
        ret <- mult.mean | mult.var | mult.ar | mult.cor | mult.lagcor
      }
    } else {
      ret <- NULL
    }
    
    if( !is.null( ret ) ){
      attr(ret, "Mpsr1") <- Mpsr1
      attr(ret, "mult.nonstat.mean") <- mult.mean
      attr(ret, "Mpsr2") <- Mpsr2
      attr(ret, "mult.nonstat.var") <- mult.var
      attr(ret, "Mpsr3") <- Mpsr3
      attr(ret, "mult.nonstat.ar") <- mult.ar
      attr(ret, "Mpsr4") <- Mpsr4
      attr(ret, "mult.nonstat.cor") <- mult.cor
      attr(ret, "Mpsr5") <- Mpsr5
      attr(ret, "mult.nonstat.lagcor") <- mult.lagcor
    }
    
  } else {
    ret <- NULL
    warning("The length of 'tseries' should be at least 3 times 'nEp', and it must be divisible by 'nEp' without a remainder.")
  }
  
  return( ret )
}

