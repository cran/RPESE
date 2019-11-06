#' @title Standard Error Estimate for Sortino Ratio of Returns
#'
#' @description \code{SortinoRatio.SE} computes the standard error of the Sortino ratio of the returns.
#'
#' @details
#' Sortino proposed an improvement on the Sharpe Ratio to better account for
#' skill and excess performance by using only downside semivariance as the
#' measure of risk.
#'
#' Sortino contends that risk should be measured in terms of not meeting the
#' investment goal.  This gives rise to the notion of \dQuote{Minimum
#' Acceptable Return} or MAR.  All of Sortino's proposed measures include the
#' MAR, and are more sensitive to downside or extreme risks than measures that
#' use volatility(standard deviation of returns) as the measure of risk.
#'
#' Choosing the MAR carefully is very important, especially when comparing
#' disparate investment choices.  If the MAR is too low, it will not adequately
#' capture the risks that concern the investor, and if the MAR is too high, it
#' will unfavorably portray what may otherwise be a sound investment.  When
#' comparing multiple investments, some papers recommend using the risk free
#' rate as the MAR.  Practitioners may wish to choose one MAR for consistency,
#' several standardized MAR values for reporting a range of scenarios, or a MAR
#' customized to the objective of the investor.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param MAR Minimum Acceptable Return for threshold.
#' @param threshold Parameter to determine whether we use a "mean" or "const" threshold.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One or a combination of:
#' \code{"IFiid"} (default), \code{"IFcor"}, \code{"IFcorPW"}, \code{"IFcorAdapt"} (default),
#' \code{"BOOTiid"} or \code{"BOOTcor"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
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
#' SortinoRatio.SE(edhec[,"CA"], se.method=c("IFiid","IFcorAdapt"),
#'                 cleanOutliers=FALSE,
#'                 fitting.method=c("Exponential", "Gamma")[1])
#'
SortinoRatio.SE <- function (data, MAR = 0, threshold = c("mean", "const")[1],
                             se.method=c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")[c(1,3)],
                             cleanOutliers=FALSE, fitting.method=c("Exponential", "Gamma")[1],
                             ...)
  { # @author Brian G. Peterson and Xin Chen
    # modified from function by Sankalp Upadhyay <sankalp.upadhyay [at] gmail [dot] com> with permission

    # Description:
    # Sortino proposed to better account for skill and excess peRformance
    # by using only downside semivariance as the measure of risk.

    # data     return vector
    # MAR   minimum acceptable return
    # Function:

    # Forcing the following parameters
    weights = NULL

    mySoR = t(apply(data, 2, SoR, threshold = threshold, ...))
    rownames(mySoR) = "SortinoRatio"

    #if we have a weights vector, use it
    if(!is.null(weights)){
      data=Return.portfolio(data,weights,...)
    }

    if(!is.null(se.method)){
      res=list(SoR=mySoR)
      # for each of the method specified in se.method, compute the standard error
      for(mymethod in se.method){
        res[[mymethod]]=EstimatorSE(data, estimator.fun = "SoR", const = MAR, threshold=threshold, MAR = MAR,
                                    se.method = mymethod,
                                    cleanOutliers=cleanOutliers,
                                    fitting.method=fitting.method,
                                    ...)
      }
      return(res)
    } else {
      return(mySoR)
    }

  }

