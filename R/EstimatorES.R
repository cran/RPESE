#'
#' @import RPEIF
#' @import PerformanceAnalytics
#' @import RPEGLMEN
#' @importFrom stats cov quantile SSD arima lm coef ar na.omit fft
#' @importFrom boot boot tsboot
#' @importFrom xts is.xts
#' @importFrom zoo zoo
#'
#' @title Wrapper Function for Standard Errors Estimates Functions
#'
#' @description \code{EstimatorSE} computes the standard error for specified risk and performance measures.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param estimator.fun Risk or performance measure to compute estimates of standard errors.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of:
#' \code{"IFiid"}, \code{"IFcor"}, \code{"IFcorAdapt"}, \code{"IFcorPW"},
#' \code{"BOOTiid"}, \code{"BOOTcor"}, or \code{"none"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
#' @param a First adaptive method parameter.
#' @param b Second adaptive method parameter.
#' @param ... Additional parameters.
#'
#' @return A vector standard error estimates.
#'
#' @export
#'
#' @author Xin Chen, \email{chenx26@uw.edu}
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @examples
#' # Loading data from PerformanceAnalytics
#' data(edhec, package = "PerformanceAnalytics")
#' class(edhec)
#' # Changing the data colnames
#' names(edhec) = c("CA", "CTA", "DIS", "EM", "EMN",
#'                  "ED", "FIA", "GM", "LS", "MA",
#'                  "RV", "SS", "FOF")
#' # Computing the standard errors for
#' # the three influence functions based approaches
#' EstimatorSE(edhec[,"CA"], se.method=c("IFcor"),
#'             cleanOutliers=FALSE,
#'             fitting.method=c("Exponential", "Gamma")[1])
#'
EstimatorSE <- function(data,
                        estimator.fun = c("Mean","SD","VaR","ES","SR","SoR","ESratio","VaRratio",
                                          "SoR","LPM","OmegaRatio","SemiSD","RachevRatio"),
                        se.method = c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor"),
                        cleanOutliers=FALSE,
                        fitting.method=c("Exponential", "Gamma")[1],
                        a=0.3, b=0.7,
                        ...){

  estimator.fun = estimator.fun[1]
  se.method = se.method[1]

  # Available estimator functions
  estimators.available <- c("Mean","SD","VaR","ES","SR","SoR",
                            "ESratio", "VaRratio", "SoR", "LPM", "OmegaRatio", "SemiSD", "RachevRatio")
  # Checking if the specified risk measure is available
  if(!(estimator.fun %in% estimators.available))
    stop("The specified estimator function is not available.")

  # Available SE methods
  se.available <- c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")
  # Checking if the standard error method is available
  if(!(se.method %in% se.available))
    stop("The specified standard error method is not available.")

  # Available fitting methods
  fitting.available <- c("Exponential", "Gamma")
  # Checking if the standard error method is available
  if(!(fitting.method %in% fitting.available))
    stop("The specified fitting method is not available.")

  myfun = switch(estimator.fun,
                 Mean = mean,
                 SD = sd,
                 VaR = VaR.hist,
                 ES = ES.hist,
                 SR = SR,
                 SoR = SoR,
                 ESratio = ESratio,
                 VaRratio = VaRratio,
                 LPM = LPM,
                 OmegaRatio = OmegaRatio,
                 SemiSD = SemiSD,
                 RachevRatio = RachevRatio,
                 stop("The estimator.fun specified is not implemented yet, please contact Anthony Christidis (anthony.christidis@stat.ubc.ca) or submit an issue at the github repository")
  )
  myfun.IF = switch (estimator.fun,
                     Mean = IF.mean,
                     SD = IF.SD,
                     VaR = IF.VaR,
                     ES = IF.ES,
                     SR = IF.SR,
                     SoR = IF.SoR,
                     ESratio = IF.ESratio,
                     VaRratio = IF.VaRratio,
                     LPM = IF.LPM,
                     OmegaRatio = IF.Omega,
                     SemiSD = IF.SSD,
                     RachevRatio = IF.RachR,
                     stop("The estimator.fun specified is not implemented yet, please contact Anthony Christidis (anthony.christidis@stat.ubc.ca) or submit an issue at the github repository")
  )

  res <- switch(se.method,
    none = NULL,
    IFiid = SE.xts(data, SE.IF.iid, myfun, myfun.IF, ...),
    IFcor = SE.xts(data, SE.IF.cor, myfun, myfun.IF,
                   prewhiten=FALSE, cleanOutliers=cleanOutliers, fitting.method=fitting.method,
                   ...),
    IFcorPW = SE.xts(data, SE.IF.cor, myfun, myfun.IF,
                   prewhiten=TRUE, cleanOutliers=cleanOutliers, fitting.method=fitting.method,
                   ...),
    IFcorAdapt = list(cor=SE.xts(data, SE.IF.cor, myfun, myfun.IF, prewhiten=FALSE,
                                 cleanOutliers=cleanOutliers, fitting.method=fitting.method,
                                 ...),
                      corPW=SE.xts(data, SE.IF.cor, myfun, myfun.IF, prewhiten=TRUE,
                                   cleanOutliers=cleanOutliers, fitting.method=fitting.method,
                                   ...)),
    BOOTiid = SE.xts(data, SE.BOOT.iid, myfun, myfun.IF, ...),
    BOOTcor = SE.xts(data, SE.BOOT.cor, myfun, myfun.IF,...)
  )

  # Adaptive method computation of weighted estimates
  if(se.method=="IFcorAdapt"){
    if(is.vector(data) || sum(ncol(data)>1)==0){
      ar1.param = arima(x=data, order=c(1,0,0), include.mean=TRUE)[[1]][1]
      if(0<=ar1.param & ar1.param<a)
        w = 0 else if(a<=ar1.param & ar1.param<=b)
          w = (ar1.param - a)/(b - a) else
            w = 1
      res$corAdapt = (1-w)*res$cor + w*res$corPW
      names(res$corAdapt) = colnames(data)
    } else{
      for(my.col in 1:ncol(data)){
        temp.data = data[, my.col]
        ar1.param = arima(x=temp.data, order=c(1,0,0), include.mean=TRUE)[[1]][1]
        if(0<=ar1.param & ar1.param<a)
          w = 0 else if(a<=ar1.param & ar1.param<=b)
            w = (ar1.param - a)/(b - a) else
              w = 1
        res$corAdapt[my.col] = (1-w)*res$cor[my.col] + w*res$corPW[my.col]
      }
      names(res$corAdapt) = colnames(data)
    }
    return(res$corAdapt)
  } else
    return(res)
}










