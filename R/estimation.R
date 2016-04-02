#' Function nginar
#'
#' Fit in inar model to a univariate time series by yule-walker method
#'
#' @param x a numeric vector or time series.
#' @param order.max the integer component p is the INAR order.
#' @param mu the mean of the geometric distribution.
#'
#'@return
#'
#'@references
#'
#'Miroslav M Ristic, Hassan S Bakouch, Aleksandar S Nastic.
#'A new geometric first-order integer-valued autoregressive (NGINAR (1)) process . 
#'Journal of Statistical Planning and Inference.VOl. 139 (2009) 2218 -- 2226
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

    residuals <- function(x,coef,mu)
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
    resid <- residuals(x,coef,mu.e)
    rms <- sqrt(mean(resid^2,na.rm = TRUE))
    
    AICc. <- n*log(Vp) + n*((1+order.max/n)/(1-(order.max+2)/n))
    AIC. <-  n*log(Vp) + 2*order.max
    BIC. <-  n*log(Vp) + (order.max/n)*log(n)
    
    res <- list(order = order.max,
                 coef = coef,
                 mean.e = mu.e,
                 var = var.error,
                 rms = rms,
                 fitted.values <- fitted,
                 bic= BIC.,
                 aicc=AICc.,
                 aic = AIC.,
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