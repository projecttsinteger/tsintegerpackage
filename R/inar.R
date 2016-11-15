#'Methods for Fitted INAR Models
#'
#'@description Methods for fitted INAR model objects.
#'
#'@usage 
#'## S3 method for class 'inar'
#'coef(object, ...)
#'## S3 method for class 'inar'
#' residuals(object, ...)
#'## S3 method for class 'inar'
#' fitted(object, ...)
#'## S3 method for class 'inar'
#' print(x, digits = max(3, getOption("digits") - 3), ...)
#' ## S3 method for class 'inar'
#' plot(x, ask = interactive(), ...)
#'
#' 
#'@param object, x an object of class "inar"; usually, a result of a call to inar.
#'@param digits see \code{\link{printCoefmat}}.
#'@param ark Should the plot method work interactively? See \code{\link{interactive}}.
#'
#'@param ...
#'further arguments passed to or from other methods.
#'
#'@return
#'For coef, a numeric vector; 
#'for residuals and fitted a univariate time series; 
#'for plot and print, the fitted INAR model object.
#'
#'@export

coef.inar <-
  function(object, ...)
  {
    if(!inherits(object, "inar"))
      stop("method is only for arma objects")
    return(object$coef)
  }

#'Methods for Fitted INAR Models
#'
#'@description Methods for fitted INAR model objects.
#'
#'@usage 
#'## S3 method for class 'inar'
#'coef(object, ...)
#'## S3 method for class 'inar'
#' residuals(object, ...)
#'## S3 method for class 'inar'
#' fitted(object, ...)
#'## S3 method for class 'inar'
#' print(x, digits = max(3, getOption("digits") - 3), ...)
#' ## S3 method for class 'inar'
#' plot(x, ask = interactive(), ...)
#'
#' 
#'@param object, x an object of class "inar"; usually, a result of a call to inar.
#'@param digits see \code{\link{printCoefmat}}.
#'@param ark Should the plot method work interactively? See \code{\link{interactive}}.
#'
#'@param ...
#'further arguments passed to or from other methods.
#'
#'@return
#'For coef, a numeric vector; 
#'for residuals and fitted a univariate time series; 
#'for plot and print, the fitted INAR model object.
#'
#'@export

residual.inar <- function(object, ...)
{
  if(!inherits(object, "inar"))
    stop("method is only for inar objects")
  return(object$resid)
}

#'Methods for Fitted INAR Models
#'
#'@description Methods for fitted INAR model objects.
#'
#'@usage 
#'## S3 method for class 'inar'
#'coef(object, ...)
#'## S3 method for class 'inar'
#' residuals(object, ...)
#'## S3 method for class 'inar'
#' fitted(object, ...)
#'## S3 method for class 'inar'
#' print(x, digits = max(3, getOption("digits") - 3), ...)
#' ## S3 method for class 'inar'
#' plot(x, ask = interactive(), ...)
#'
#' 
#'@param object, x an object of class "inar"; usually, a result of a call to inar.
#'@param digits see \code{\link{printCoefmat}}.
#'@param ark Should the plot method work interactively? See \code{\link{interactive}}.
#'
#'@param ...
#'further arguments passed to or from other methods.
#'
#'@return
#'For coef, a numeric vector; 
#'for residuals and fitted a univariate time series; 
#'for plot and print, the fitted INAR model object.
#'
#'@export

fitted.inar <-
  function(object, ...)
  {
    if(!inherits(object, "inar"))
      stop("method is only for arma objects")
    return(object$fitted.values)
  }

#'Methods for Fitted INAR Models
#'
#'@description Methods for fitted INAR model objects.
#'
#'@usage 
#'## S3 method for class 'inar'
#'coef(object, ...)
#'## S3 method for class 'inar'
#' residuals(object, ...)
#'## S3 method for class 'inar'
#' fitted(object, ...)
#'## S3 method for class 'inar'
#' print(x, digits = max(3, getOption("digits") - 3), ...)
#' ## S3 method for class 'inar'
#' plot(x, ask = interactive(), ...)
#'
#' 
#'@param object, x an object of class "inar"; usually, a result of a call to inar.
#'@param digits see \code{\link{printCoefmat}}.
#'@param ark Should the plot method work interactively? See \code{\link{interactive}}.
#'
#'@param ...
#'further arguments passed to or from other methods.
#'
#'@return
#'For coef, a numeric vector; 
#'for residuals and fitted a univariate time series; 
#'for plot and print, the fitted INAR model object.
#'
#'@export

print.inar <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(!inherits(x, "inar"))
      stop("method is only for inar objects")
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    cat("\nModel:\nINAR(",x$order.max,")\n", sep = "")
    cat("\nResiduals:\n")
    rq <- structure(quantile(x$resid,na.rm = TRUE),
                    names = c("Min","1Q","Median","3Q","Max"))
    print(rq, digits = digits, ...)
    cat("\n Coefficient(s):\n")
    print.default(format(coef(x), digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\nFit:\n")
    cat("RMS =", format(x$rms, digits = digits),
        ",  mean.error = ", format(round(x$mean.e, 2)),
        ",  var.error = ", format(round(x$var, 2)), "\n", sep = "")
    
    # cat("\nInformation Criterion:\n")
    #  cat("AIC =", format(round(x$bic, 2)),
    #      ", BIC = ", format(round(x$bic, 2)),
    #     ",  AICc = ", format(round(x$aicc, 2)), "\n", sep = "")
    
    cat("\n")
    invisible(x)
  }

#'Methods for Fitted INAR Models
#'
#'@description Methods for fitted INAR model objects.
#'
#'@usage 
#'## S3 method for class 'inar'
#'coef(object, ...)
#'## S3 method for class 'inar'
#' residuals(object, ...)
#'## S3 method for class 'inar'
#' fitted(object, ...)
#'## S3 method for class 'inar'
#' print(x, digits = max(3, getOption("digits") - 3), ...)
#' ## S3 method for class 'inar'
#' plot(x, ask = interactive(), ...)
#'
#' 
#'@param object, x an object of class "inar"; usually, a result of a call to inar.
#'@param digits see \code{\link{printCoefmat}}.
#'@param ark Should the plot method work interactively? See \code{\link{interactive}}.
#'
#'@param ...
#'further arguments passed to or from other methods.
#'
#'@return
#'For coef, a numeric vector; 
#'for residuals and fitted a univariate time series; 
#'for plot and print, the fitted INAR model object.
#'
#'@export

summary.inar <-
  function(object, ...)
  {
    if(!inherits(object, "inar"))
      stop("method is only for arma objects")
    return(summary(object$resid))
  }

#'Methods for Fitted INAR Models
#'
#'@description Methods for fitted INAR model objects.
#'
#'@usage 
#'## S3 method for class 'inar'
#'coef(object, ...)
#'## S3 method for class 'inar'
#' residuals(object, ...)
#'## S3 method for class 'inar'
#' fitted(object, ...)
#'## S3 method for class 'inar'
#' print(x, digits = max(3, getOption("digits") - 3), ...)
#' ## S3 method for class 'inar'
#' plot(x, ask = interactive(), ...)
#'
#' 
#'@param object, x an object of class "inar"; usually, a result of a call to inar.
#'@param digits see \code{\link{printCoefmat}}.
#'@param ark Should the plot method work interactively? See \code{\link{interactive}}.
#'
#'@param ...
#'further arguments passed to or from other methods.
#'
#'@return
#'For coef, a numeric vector; 
#'for residuals and fitted a univariate time series; 
#'for plot and print, the fitted INAR model object.
#'
#'@export

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