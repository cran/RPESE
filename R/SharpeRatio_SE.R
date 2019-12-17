#' @title Standard Error Estimate for Sharpe Ratio of Returns
#'
#' @description \code{SharpeRatio.SE} computes the standard error of the Sharpe ratio of the returns.
#'
#' @details
#' The Sharpe ratio is simply the return per unit of risk (represented by
#' variability).  In the classic case, the unit of risk is the standard
#' deviation of the returns.
#'
#' \deqn{\frac{\overline{(R_{a}-R_{f})}}{\sqrt{\sigma_{(R_{a}-R_{f})}}}}
#'
#' William Sharpe now recommends \code{\link{InformationRatio}} preferentially
#' to the original Sharpe Ratio.
#'
#' The higher the Sharpe ratio, the better the combined performance of "risk"
#' and return.
#'
#' As noted, the traditional Sharpe Ratio is a risk-adjusted measure of return
#' that uses standard deviation to represent risk.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param Rf Risk free rate.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One or a combination of:
#' \code{"IFiid"} (default), \code{"IFcor"}, \code{"IFcorPW"}, \code{"IFcorAdapt"} (default),
#' \code{"BOOTiid"} or \code{"BOOTcor"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
#' @param freq.include Frequency domain inclusion criteria. Must be one of "All" (default), "Decimate" or "Truncate."
#' @param freq.par Percentage of the frequency used if \code{"freq.include"} is "Decimate" or "Truncate." Default is 0.5.
#' @param ... Additional parameters.
#'
#' @return A vector or a list depending on \code{se.method}.
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
#' # the two influence functions based approaches
#' SharpeRatio.SE(edhec, se.method=c("IFiid","IFcorAdapt"),
#'                cleanOutliers=FALSE,
#'                fitting.method=c("Exponential", "Gamma")[1])
#'
SharpeRatio.SE <- function (data, Rf = 0,
                            se.method=c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")[c(1,3)],
                            cleanOutliers=FALSE, fitting.method=c("Exponential", "Gamma")[1],
                            freq.include=c("All", "Decimate", "Truncate")[1], freq.par=0.5,
                            ...){

    # Forcing the following parameters
    FUN = FUN=c("StdDev", "VaR","ES")[1]
    weights=NULL
    annualize = FALSE

    mySR = t(apply(data, 2, SR, rf = Rf, ...))
    rownames(mySR) = "SharpeRatio"

    if(length(FUN)==1 & FUN=="StdDev" & is.null(weights) & annualize == FALSE & !is.null(se.method)){

      data = checkData(data)

      if(!is.null(dim(Rf)))
        Rf = checkData(Rf)

      res=list(SR=mySR)
      # for each of the method specified in se.method, compute the standard error
      for(mymethod in se.method){
        res[[mymethod]]=EstimatorSE(data, estimator.fun = "SR", rf = Rf,
                                    se.method = mymethod,
                                    cleanOutliers=cleanOutliers,
                                    fitting.method=fitting.method,
                                    freq.include=freq.include, freq.par=freq.par,
                                    ...)
      }
      return(res)
    } else {
      return(mySR)
    }
  }
