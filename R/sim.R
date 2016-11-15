#' Function poinar.sim
#'
#' Simulate from an INAR model
#'
#' @param n the length of outputs series. A strictly positive integer.
#' @param order.max the integer component p is the INAR order.
#' @param alpha a vector of INAR coefficients.
#' @param lambda the mean of the poisson distribution.
#' @param n.start the length of 'burn-in' period. If na, the default, a reasonable valve is computed.
#'
#'
#'@return A time-series object of class "ts".
#'
#'@seealso   \code{\link{poinar}}
#'
#'@references
#'  Du, J.G. and Li,Y. (1991).
#'  The integer-valued autorregressive (INAR(p)) model.
#'  \emph{Journal of time series analysis}. \bold{12}, 129--142.
#'
#'@examples
#'# A Poisson INAR simulation
#'ts.sim <- poinar.sim(n = 100, order.max = 2, alpha = c(0.1,0.4),lambda = 2, n.start=200)
#'ts.plot(ts.sim)
#'
#' @export
#' 
poinar.sim <- function(n, order.max, alpha,lambda, n.start=NA){
  length. <- n + n.start
  x <- rep(NA, times = length.)
  error <- rpois(length., lambda)
  for (i in 1:order.max) {
    x[i] <- error[i]
  }
  for (t in (order.max + 1):length.) {
    x[t] <- 0
    for (j in 1:order.max) {
      x[t] <- x[t] + rbinom(1, x[t - j], alpha[j])
    }
    x[t] <- x[t] + error[t]
  }
  ts(x[(n.start+1):length.],frequency = 1,start=1)
}



#' Function nginar.sim
#'
#' Simulate from an Inar model
#'
#' @param n A strictly positive integer.
#' @param alpha a vector of INAR coefficients.
#' @param mu a vector.
#' @param n.start the length of 'burn-in' period. If na, the default, a reasonable valve is computed.
#'
#'@return A time-series object of class "ts".
#'
#'@seealso   \code{\link{nginar}}
#'
#'@references
#'
#'Ristic, M.M., Bakouchb, H.S. and Nasti, A.S. (2009) 
#'A newgeometricfirst-orderinteger-valuedautoregressive (NGINAR(1)) process . 
#'\emph{Journal of Statistical Planning and Inference}, \bold{139}, 2218--2226.
#'
#'@examples
#'
#'# A new geometric INAR simulation
#'ts.sim <- nginar.sim(n = 100, alpha = 0.4, mu = 2)
#'ts.plot(ts.sim)
#'
#' @export
#' 
nginar.sim <- function(n, alpha,mu, n.start=150){
  length. <- n + n.start
  error.nginar <- function(length.,alpha,mu){
    epsilon <- rep(NA, times = n)
    a = (alpha*mu)/(mu - alpha)
    for(i in 1:length.){
      u <- runif(1,0,1)
      ifelse(u < a,epsilon[i] <- rgeom(1, 1 - (alpha/(1 + alpha))),epsilon[i] <- rgeom(1, 1 - (mu/(1+mu))))
    }
    return(epsilon)
  }
  
  x <- rep(NA, times = length.)
  epsilon <- error.nginar(length.,alpha,mu)
  x[1] = rgeom(1, 1-(mu/(1+mu)))
  for (t in 2:length.){
    if(x[t-1]==0){
      x[t] = epsilon[t]
    }
    else{
      x[t] = rnbinom(n = 1, size = x[t-1], prob = 1 - (alpha/(1+alpha))) + epsilon[t]
    }
  }
  ts(x[(n.start+1):length.],frequency = 1,start=1)
}
