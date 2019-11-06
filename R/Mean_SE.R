#' @title Standard Error Estimate for Mean of Returns
#'
#' @description \code{Mean.SE} computes the standard error of the mean of the returns.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One or a combination of:
#' \code{"IFiid"} (default), \code{"IFcor"}, \code{"IFcorPW"}, \code{"IFcorAdapt"} (default),
#' \code{"BOOTiid"} or \code{"BOOTcor"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
#' @param ... Additional parameters.
#'
#' @return A vector or a list depending on se.method
#'
#' @export
#'
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
#' # the two influence functions based approaches
#' Mean.SE(edhec, se.method=c("IFiid","IFcorAdapt"),
#'         cleanOutliers=FALSE,
#'         fitting.method=c("Exponential", "Gamma")[1])
#'
Mean.SE <- function(data, se.method=c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")[c(1,3)],
                    cleanOutliers=FALSE, fitting.method=c("Exponential", "Gamma")[1],
                    ...){
  data = checkData(data)
  myMean = t(apply(data, 2, mean))
  rownames(myMean) = "Mean"
  if(is.null(se.method)){
    return(myMean)
  } else{
    res=list(Mean=myMean)
    # for each of the method specified in se.method, compute the standard error
    for(mymethod in se.method){
      res[[mymethod]]=EstimatorSE(data, estimator.fun = "Mean",
                                  se.method=mymethod,
                                  cleanOutliers=cleanOutliers,
                                  fitting.method=fitting.method,
                                  ...)
    }
    return(res)
  }
}
