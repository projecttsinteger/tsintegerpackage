otest <- function(x, conf.level = 0.05)
  {
  alpha <- conf.level
  DNAME <- deparse(substitute(data))
  ans <- NULL
  METHOD <- "Test for overdispersion"
  n         <- length(x)
  alpha.hat <- acf(data)$acf[2]
  z         <- qnorm(1 - conf.level)
  r         <- sqrt((2*(1 + alpha.hat^2))/(n*(1 - alpha.hat^2)))
  s         <- sqrt((n*(1 - alpha.hat^2))/(2*(1 + alpha.hat^2)))
  STATISTIC   <- 1 + z*r
  PARAMETER <- ((n-1)/(n))*var(x)/mean(x)
  PVAL <- 1 - pnorm(s*(PARAMETER - 1))
  ans$statistic <- STATISTIC
  ans$p.value <- PVAL
  ans$method <- METHOD
  ans$parameter <- PARAMETER
  ans$data.name <- DNAME
  ans$coef <- cbind(alpha,PARAMETER,STATISTIC,PVAL)
  dimnames(ans$coef) <- list( "values" , c("alpha", "ID", "Statistic","Pr(>|z|)"))
  class(ans) <- "otest"
  return(ans)
}

print.otest <- function(x, digits = max(3, getOption("digits") - 2),signif.stars = getOption("show.signif.stars"), ...)
{
  if(!inherits(x, "otest"))
  stop("method is only for otest objects")
  cat("\n\t", x$method, "\n\n")
  cat("data: ", x$data.name, "\n")
  cat("\nResult(s):\n")
  printCoefmat(x$coef, digits = digits,signif.stars = signif.stars, ...)
  cat("\n")
  invisible(x)
}




  