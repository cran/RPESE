# Compute the standard error for the xts object
SE.xts = function(x, se.fun, myfun, myfun.IF,
                  prewhiten=FALSE, cleanOutliers=FALSE, fitting.method=c("Exponential", "Gamma")[1],
                  ...){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(se.fun(x = x, myfun = myfun, myfun.IF = myfun.IF,
                  prewhiten=prewhiten, cleanOutliers=cleanOutliers, fitting.method=fitting.method,
                  ...))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, se.fun, myfun = myfun, myfun.IF = myfun.IF,
                 prewhiten=prewhiten, cleanOutliers=cleanOutliers, fitting.method=fitting.method,
                 ...))
  }
}

# Compute the standard error using influence function approach for a vector
SE.IF.iid = function(x, myfun.IF, ...){
  N=length(x)
  x.IF = myfun.IF(x, ...)
  x.IF.2 = x.IF^2
  tmp = mean(x.IF.2)
  return(sqrt(tmp/N))
}

# Compute the standard error of the measure by iid bootstrapping
SE.BOOT.iid = function(x, myfun, myfun.IF, prewhiten=FALSE, ..., nsim = 100){
  res = boot(data = x, statistic = function(x,i,...) myfun(x[i],...), R = nsim, ... = ...)
  return(sd(res$t))
}

# Compute the standard error of the measure by tsboot()
SE.BOOT.cor = function(x, myfun, myfun.IF, prewhiten=FALSE, ..., nsim = 1000,
                       sim = "fixed", l = round(length(x)/5)){
  res = tsboot(tseries = x, statistic = function(x,...) myfun(x,...), R = nsim,
               sim = sim, l = l,...)
  return(sd(res$t))
}

# Compute the standard error using GLM-EN approach for serially correlated data using RPEGLMEN
SE.IF.cor = function(x, myfun.IF, return.coeffs = FALSE, d.GLM.EN = 5, alpha.EN = 0.5, keep = 1,
                     standardize = FALSE,
                     prewhiten=FALSE, cleanOutliers=FALSE, fitting.method=c("Exponential", "Gamma")[1],
                     ...){
  d = d.GLM.EN
  data.IF = myfun.IF(x, prewhiten=prewhiten, cleanOutliers=cleanOutliers, ...)
  if(prewhiten){
    ar.coeffs <- as.numeric(arima(x=x, order=c(1,0,0), include.mean=TRUE)$coef[1])
  } else{
    ar.coeffs <- NULL
  }
  tmp = SE.glmnet_exp(data.IF, standardize = standardize,
                      return.coeffs = return.coeffs, d = d, alpha.EN = alpha.EN, keep = keep,
                      fitting.method = fitting.method,
                      prewhiten = prewhiten, ar.coeffs = ar.coeffs,
                      ...)
  if(return.coeffs){
    coeffs = tmp[[2]]
    tmp = tmp[[1]]
    return(list(sqrt(tmp), coeffs))
  }
  return(sqrt(tmp))
}
