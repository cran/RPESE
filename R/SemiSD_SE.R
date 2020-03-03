#' @title Standard Error Estimate for Semi-Standared Deviation (SemiSD) of Returns
#'
#' @description \code{SemiSD.SE} computes the standard error of the SSD of the returns.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param rf Risk-free interest rate.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One or a combination of:
#' \code{"IFiid"} (default), \code{"IFcor"} (default), \code{"IFcorPW"}, \code{"IFcorAdapt"},
#' \code{"BOOTiid"}, \code{"BOOTcor"}, or \code{"none"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
#' @param d.GLM.EN Order of the polynomial for the Exponential or Gamma fitting. Default polynomial order of 5.
#' @param freq.include Frequency domain inclusion criteria. Must be one of "All" (default), "Decimate" or "Truncate."
#' @param freq.par Percentage of the frequency used if \code{"freq.include"} is "Decimate" or "Truncate." Default is 0.5.
#' @param ... Additional parameters.
#'
#' @return A vector or a list depending on \code{se.method}.
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
#' SemiSD.SE(edhec, se.method=c("IFiid","IFcor"),
#'           cleanOutliers=FALSE,
#'           fitting.method=c("Exponential", "Gamma")[1])
#'
SemiSD.SE <- function(data, rf=0,
                      se.method=c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")[1:2],
                      cleanOutliers=FALSE, fitting.method=c("Exponential", "Gamma")[1], d.GLM.EN = 5,
                      freq.include=c("All", "Decimate", "Truncate")[1], freq.par=0.5,
                      ...){
  data = checkData(data)
  mySSD = t(apply(data, 2, SemiSD, rf, ...))
  rownames(mySSD) = "SSD"
  if(is.null(se.method) & length(se.method)==1){
    return(mySSD)
  } else {
    res=list(SSD=mySSD)
    # for each of the method specified in se.method, compute the standard error
    for(mymethod in se.method){
      res[[mymethod]]=EstimatorSE(data, estimator.fun = "SemiSD",
                                  rf=rf,
                                  se.method = mymethod,
                                  cleanOutliers=cleanOutliers,
                                  fitting.method=fitting.method, d.GLM.EN=d.GLM.EN,
                                  freq.include=freq.include, freq.par=freq.par,
                                  ...)
    }
    return(res)
  }
}
