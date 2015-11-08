#' Function poinar
#'
#' Fit in inar model to a univariate time series by yule-walker method
#'
#' @param x a numeric vector or time series.
#' @param order.max a one dimensional integer vestor giving the order of the model to fit. This value corresponds the INAR order.
#' @param series name for the series.
#'
#'@return Resultados
#'
#'@references
#'
#'Du, J.G. and Li,Y.(1991):
#'The integer-valued autorregressive (INAR(p)) model.
#'\emph{Journal of time series analysis} \bold{12}, 129--142.
#'
#'Freeland RK. Statistical analysis of discrete time series with applications to the analysis of workers compensation
#'claims data [unpublished doctoral dissertation]. Vancouver (Canada): University of British Columbia; 1998.
#'
#'@examples
#'data(claims)
#'mean(claims[,5])
#'var(claims[,5])
#'var(claims[,5])/mean(claims[,5])  # dispersion index
#'acf(claims[,5])
#'pacf(claims[,5])
#'poinar(claims[,5], 1)
#'
#' @export
poinar <-
  function(x, order.max,series=NULL)
  {
    if (is.null(series))
      series <- deparse(substitute(x))
    xfreq <- frequency(x)
    n <- length(x)
    
    r <- acf(x, plot = FALSE)$acf[2:(order.max+1)]
    R <- diag(order.max)
    for(i in 1:order.max){
      for(j in 1:order.max){
        if(i!=j){
          R[i,j] <- r[abs(i-j)]
        }
      }
    }
    
    residuals <- function(x,coef,p,lambda)
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
    mu.e <- var.error <- xbar*(1-sum(coef)) #mean and variance of the error
    mu.x <- mu.e/(1-sum(coef))
    sum.var <- sum( coef*(1-coef))
    Vp <- var.error + mu.x*sum.var
    fitted <- poinar.sim(length(x), order.max = length(coef), alpha = coef,lambda = mu.e, n.start=200)
    resid <- residuals(x,coef,order.max,mu.e)
    r <- na.remove(resid,PACKAGE = )
    rms <- sqrt(mean(r^2))
    
    AICc. <- n*log(Vp) + n*((1+order.max/n)/(1-(order.max+2)/n))
    AIC. <-  -n*log(Vp) + 2*order.max
    BIC. <-  n*log(Vp) + (order.max/n)*log(n)
    
    inar <- list(order = order.max,
                coef = coef,
                lambda= mu.e,
                x.mean = xbar,
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
    class(inar) <- "inar"
    return(inar)
  }



coef.inar <-
  function(object, ...)
  {
    if(!inherits(object, "inar"))
      stop("method is only for arma objects")
    return(object$coef)
  }


residuals.inar <-
  function(object, ...)
  {
    if(!inherits(object, "inar"))
      stop("method is only for arma objects")
    return(object$resid)
  }

fitted.inar <-
  function(object, ...)
  {
    if(!inherits(object, "inar"))
      stop("method is only for arma objects")
    return(object$fitted.values)
  }

print.inar <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(!inherits(x, "inar"))
      stop("method is only for inar objects")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Coefficient(s):\n")
    print.default(format(coef(x), digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\nFit:\n")
    cat("RMS =", format(x$rms, digits = digits),
        ",  lambda estimated as = ", format(round(x$lambda, 2)),
        ",  AIC = ", format(round(x$aic, 2)), "\n", sep = "")
    cat("\n")
    invisible(x)
  }

summary.inar <-
  function(object, ...)
  {
    if(!inherits(object, "inar"))
      stop("method is only for arma objects")
    ans <- NULL
    ans$residuals <- na.remove(object$resid)
    ans$criterion <- c(object$aic,object$aicc,object$bic)
    ans$coef <- cbind(object$coef, c(object$aic,object$aicc,object$bic)  )
    dimnames(ans$coef) <- list(names(object$coef), c(" Estimate"," Inf. Criterion"))
    ans$call <- object$call
    ans$nn <- object$nn
    ans$var <- var(ans$residuals)
    ans$rms <- object$rms
    ans$p <- object$order.max
    class(ans) <- "summary.inar"
    return(ans)
  }

print.summary.arma <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(!inherits(x, "summary.inar"))
      stop("method is only for summary.inar objects")
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    cat("\nModel:\nINAR(",x$p,")\n", sep = "")
    cat("\nResiduals:\n")
    rq <- structure(quantile(x$residuals),
                    names = c("Min","1Q","Median","3Q","Max"))
    print(rq, digits = digits, ...)
    cat("\nCoefficient(s):\n")
    printCoefmat(x$coef, digits = digits, ...)
    cat("\nFit:\n")
    cat("sigma^2 estimated as ", format(x$var, digits = digits), 
        ",  RMS = ", format(round(x$rms, 2)) , "\n", sep = "")
    cat("\n")
    invisible(x)
  }


plot.inar <-
  function(x, ask = interactive(), ...)
  {
    if(!inherits(x, "inar"))
      stop("method is only for inar objects")
    op <- par()
    par(ask = ask, mfrow = c(2, 1))
    data <- eval.parent(parse(text = x$series))
    if(any(is.na(data))) stop(paste("NAs in", x$series))
    plot(data, main = x$series, ylab = "Series")
    plot(x$resid, main = "Residuals", ylab = "Series")
    acf(data, main = paste("ACF of", x$series))
    acf(x$resid, main = "ACF of Residuals",
        na.action = na.remove)
    pacf(data, main = paste("PACF of", x$series))
    pacf(x$resid, main = "PACF of Residuals",
         na.action = na.remove)
    par(ask = op$ask, mfrow = op$mfrow)
    invisible(x)
  }