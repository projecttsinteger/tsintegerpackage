#'Function expoinar
#'
#'Fit in inar model to a univariate time series by yule-walker method
#'
#'@usage
#'expoinar(x, order.max,series=NULL,family="gen.poisson")
#'
#'@param x a numeric vector or time series.
#'@param order.max a one dimensional integer vestor giving the order of the model to fit. This value corresponds the INAR order.
#'@param series name for the series.
#'
#'@return
#'A list of class "inar" with the following elements:
#'
#'\describe{
#'\item{coef:}{estimated INAR coefficients for the fitted model.}
#'\item{residual:}{the series of residuals.}
#'\item{fitted.values:}{the fitted series.}
#'\item{series:}{the name of the series x.}
#'}
#'@references
#'
#'Du, J.G. and Li,Y. (1991).
#'The integer-valued autorregressive (INAR(p)) model.
#'\emph{Journal of time series analysis}. \bold{12}, 129--142.
#'
#'Freeland R. K. (1998). 
#'\emph{Statistical analysis of discrete time series with applications to the analysis of workers compensation
#'claims data [unpublished doctoral dissertation]}. Vancouver (Canada): University of British Columbia.
#'
#'@examples
#'data(claims)
#'claims5 <- claims[,5]
#'mean(claims5)
#'var(claims5)
#'var(claims5)/mean(claims5)  # dispersion index
#'acf(claims5)
#'pacf(claims5)
#'poinar(claims5, 1)
#'
#'@export
#'
#'


#'Function poinar
#'
#'Fit in inar model to a univariate time series by yule-walker method
#'
#'@usage
#'poinar(x, order.max,series=NULL)
#'
#'@param x a numeric vector or time series.
#'@param order.max a one dimensional integer vestor giving the order of the model to fit. This value corresponds the INAR order.
#'@param series name for the series.
#'
#'@return
#'A list of class "inar" with the following elements:
#'
#'\describe{
#'\item{coef:}{estimated INAR coefficients for the fitted model.}
#'\item{residual:}{the series of residuals.}
#'\item{fitted.values:}{the fitted series.}
#'\item{series:}{the name of the series x.}
#'}
#'@references
#'
#'Du, J.G. and Li,Y. (1991).
#'The integer-valued autorregressive (INAR(p)) model.
#'\emph{Journal of time series analysis}. \bold{12}, 129--142.
#'
#'Freeland R. K. (1998). 
#'\emph{Statistical analysis of discrete time series with applications to the analysis of workers compensation
#'claims data [unpublished doctoral dissertation]}. Vancouver (Canada): University of British Columbia.
#'
#'@examples
#'data(claims)
#'claims5 <- claims[,5]
#'mean(claims5)
#'var(claims5)
#'var(claims5)/mean(claims5)  # dispersion index
#'acf(claims5)
#'pacf(claims5)
#'poinar(claims5, 1)
#'
#'@export

poinar <-
  function(x, order.max,series=NULL)
  {
    if (is.null(series))
      series <- deparse(substitute(x))
    xfreq <- frequency(x)
    n <- length(x)
    r0 <- acf(x, plot = FALSE)$acf[1]
    r <- acf(x, plot = FALSE)$acf[2:(order.max+1)]
    R <- diag(order.max)
    for(i in 1:order.max){
      for(j in 1:order.max){
        if(i!=j){
          R[i,j] <- r[abs(i-j)]
        }
      }
    }
    
    residual <-
      function(x,coef,p,lambda)
      {
        x<- x
        e <- NULL
        for(t in (p+1):length(x) )
        {
          e[t] <- x[t] - ( sum(coef*x[(t-1):(t-p)]) + lambda)
        }
        return(e)
      }
    
    coef <- round(solve(R, r), 4)
    xbar <- mean(x)	#mean of the serie
    mu.e <- xbar*(1-sum(coef)) #mean and variance of the error
    var.error <- r0 - sum(coef*r)
    mu.x <- mu.e/(1-sum(coef))
    sum.var <- sum( coef*(1-coef))
    Vp <- var.error + mu.x*sum.var
    fitted <- poinar.sim(length(x), order.max = length(coef), alpha = coef,lambda = mu.e, n.start=200)
    resid <- residual(x,coef,order.max,mu.e)
    rms <- sqrt(mean(resid^2,na.rm = TRUE))
    
    #AICc. <- n*log(Vp) + n*((1+order.max/n)/(1-(order.max+2)/n))
    #AIC. <-  n*log(Vp) + 2*order.max
    #BIC. <-  n*log(Vp) + (order.max/n)*log(n)
    
    inar <- list(order = order.max,
                 coef = coef,
                 mean.e = mu.e,
                 var = var.error,
                 rms = rms,
                 fitted.values <- fitted,
                 #        bic= BIC.,
                 #        aicc=AICc.,
                 #        aic = AIC.,
                 n.used = n,
                 order.max = order.max,
                 resid = resid,
                 method = "yule-walker",
                 series = series,
                 frequency = xfreq,
                 call = match.call())
    class(inar) <- "inar"
    return(inar)
  }



#' Function nginar
#'
#' Fit in inar model to a univariate time series by yule-walker method
#'
#' @param x a numeric vector or time series.
#' @param order.max the integer component p is the NGINAR order.
#' @param mu the mean of the geometric distribution.
#'
#'@return
#'A list of class "inar" with the following elements:
#'
#'\describe{
#'\item{coef:}{estimated NGINAR coefficients for the fitted model.}
#'\item{residual:}{the series of residuals.}
#'\item{fitted.values:}{the fitted series.}
#'\item{series:}{the name of the series x.}
#'}
#'
#'@references
#'
#'Ristic, M.M., Bakouchb, H.S. and Nasti, A.S. (2009) 
#'A newgeometricfirst-orderinteger-valuedautoregressive (NGINAR(1)) process . 
#'\emph{Journal of Statistical Planning and Inference}, \bold{139}, 2218--2226.
#'
#'@examples
#'
#'data(sexoffences)
#'mean(sexoffences)
#'var(sexoffences)
#'nginar(sexoffences)
#'
#'@export

nginar <- function(x,series=NULL)
{
  if (is.null(series))
    series <- deparse(substitute(x))
  xfreq <- frequency(x)
  n <- length(x)
  order.max = 1
  r0 <- acf(x, plot = FALSE)$acf[1]
  r <- acf(x, plot = FALSE)$acf[2:(order.max+1)]
  R <- diag(order.max)
  for(i in 1:order.max){
    for(j in 1:order.max){
      if(i!=j){
        R[i,j] <- r[abs(i-j)]
      }
    }
  }
  
  residual <- function(x,coef,mu)
  {
    x<- x
    p<- 1
    e <- NULL
    for(t in (p+1):length(x) )
    {
      e[t] <- x[t] - ( sum(coef*x[(t-1):(t-p)]) + mu)
    }
    return(e)
  }
  
  coef <- round(solve(R, r), 4)
  xbar <- mean(x)
  mu.e <- xbar*(1-sum(coef)) #mean of the error
  var.error <- r0 - sum(coef*r) #variance of the erro
  mu.x <- mu.e/(1-sum(coef))
  sum.var <- sum( coef*(1-coef))
  Vp <- var.error + mu.x*sum.var
  fitted <- nginar.sim(length(x), alpha = coef, mu = mu.e, n.start=200)
  resid <- residual(x,coef,mu.e)
  rms <- sqrt(mean(resid^2,na.rm = TRUE))
  
  # AICc. <- n*log(Vp) + n*((1+order.max/n)/(1-(order.max+2)/n))
  # AIC. <-  n*log(Vp) + 2*order.max
  #  BIC. <-  n*log(Vp) + (order.max/n)*log(n)
  
  res <- list(order = order.max,
              coef = coef,
              mean.e = mu.e,
              var = var.error,
              rms = rms,
              fitted.values <- fitted,
              #  bic= BIC.,
              #  aicc=AICc.,
              #   aic = AIC.,
              n.used = n,
              order.max = order.max,
              resid = resid,
              method = "yule-walker",
              series = series,
              frequency = xfreq,
              call = match.call())
  class(res) <- "inar"
  return(res)
}