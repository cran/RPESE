#' @title Standard Error Estimate for Value-at-Risk (VaR) of Returns
#'
#' @description \code{VaR.SE} computes the standard error of the value-at-risk of the returns.
#'
#' @note
#' The option to \code{invert} the VaR measure should appease both
#' academics and practitioners.  The mathematical definition of VaR as the
#' negative value of a quantile will (usually) produce a positive number.
#' Practitioners will argue that VaR denotes a loss, and should be internally
#' consistent with the quantile (a negative number).  For tables and charts,
#' different preferences may apply for clarity and compactness.  As such, we
#' provide the option, and set the default to TRUE to keep the return
#' consistent with prior versions of PerformanceAnalytics, but make no value
#' judgment on which approach is preferable.
#'
#' The prototype of the univariate Cornish Fisher VaR function was completed by
#' Prof. Diethelm Wuertz.  All corrections to the calculation and error
#' handling are the fault of Brian Peterson.
#'
#' @section
#' Background: This function provides several estimation methods for
#' the Value at Risk (typically written as VaR) of a return series and the
#' Component VaR of a portfolio. Take care to capitalize VaR in the commonly
#' accepted manner, to avoid confusion with var (variance) and VAR (vector
#' auto-regression).  VaR is an industry standard for measuring downside risk.
#' For a return series, VaR is defined as the high quantile (e.g. ~a 95% or 99%
#' quantile) of the negative value of the returns. This quantile needs to be
#' estimated.  With a sufficiently large data set, you may choose to utilize
#' the empirical quantile calculated using \code{\link{quantile}}.  More
#' efficient estimates of VaR are obtained if a (correct) assumption is made on
#' the return distribution, such as the normal distribution.  If your return
#' series is skewed and/or has excess kurtosis, Cornish-Fisher estimates of VaR
#' can be more appropriate.  For the VaR of a portfolio, it is also of interest
#' to decompose total portfolio VaR into the risk contributions of each of the
#' portfolio components.  For the above mentioned VaR estimators, such a
#' decomposition is possible in a financially meaningful way.
#'
#' @references Boudt, Kris, Peterson, Brian, and Christophe Croux. 2008.
#' Estimation and decomposition of downside risk for portfolios with non-normal
#' returns. 2008. The Journal of Risk, vol. 11, 79-103.
#'
#' Cont, Rama, Deguest, Romain and Giacomo Scandolo. Robustness and sensitivity
#' analysis of risk measurement procedures. Financial Engineering Report No.
#' 2007-06, Columbia University Center for Financial Engineering.
#'
#' Denton M. and Jayaraman, J.D. Incremental, Marginal, and Component VaR.
#' Sunguard. 2004.
#'
#' Epperlein, E., Smillie, A. Cracking VaR with kernels. RISK, 2006, vol.  19,
#' 70-74.
#'
#' Gourieroux, Christian, Laurent, Jean-Paul and Olivier Scaillet.  Sensitivity
#' analysis of value at risk. Journal of Empirical Finance, 2000, Vol. 7,
#' 225-245.
#'
#' Keel, Simon and Ardia, David. Generalized marginal risk. Aeris CAPITAL
#' discussion paper.
#'
#' Laurent Favre and Jose-Antonio Galeano. Mean-Modified Value-at-Risk
#' Optimization with Hedge Funds. Journal of Alternative Investment, Fall 2002,
#' v 5.
#'
#' Martellini, Lionel, and Volker Ziemann.  Improved Forecasts of Higher-Order
#' Comoments and Implications for Portfolio Selection. 2007. EDHEC Risk and
#' Asset Management Research Centre working paper.
#'
#' Zangari, Peter. A VaR Methodology for Portfolios that include Options. 1996.
#' RiskMetrics Monitor, First Quarter, 4-12.
#'
#' Rockafellar, Terry and Uryasev, Stanislav. Optimization of Conditional VaR.
#' The Journal of Risk, 2000, vol. 2, 21-41.
#'
#' Dowd, Kevin. Measuring Market Risk, John Wiley and Sons, 2010.
#'
#' Jorian, Phillippe. Value at Risk, the new benchmark for managing financial risk.
#' 3rd Edition, McGraw Hill, 2006.
#'
#' Hallerback, John. "Decomposing Portfolio Value-at-Risk: A General Analysis",
#' 2003. The Journal of Risk vol 5/2.
#'
#' Yamai and Yoshiba (2002). "Comparative Analyses of Expected Shortfall and
#'    Value-at-Risk: Their Estimation Error, Decomposition, and Optimization",
#'    Bank of Japan.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param p Confidence level for calculation. Default is p=0.95.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One or a combination of:
#' \code{"IFiid"} (default), \code{"IFcor"} (default), \code{"IFcorPW"}, \code{"IFcorAdapt"},
#' \code{"BOOTiid"} or \code{"BOOTcor"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
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
#' VaR.SE(edhec, se.method=c("IFiid","IFcor"),
#'        cleanOutliers=FALSE,
#'        fitting.method=c("Exponential", "Gamma")[1])
#'
VaR.SE <- function(data=NULL , p=0.95,
                   se.method=c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")[1:2],
                   cleanOutliers=FALSE, fitting.method=c("Exponential", "Gamma")[1],
                   ...)
  { # @author Brian G. Peterson and Xin Chen

    # Descripion:

    # wrapper for univariate and multivariate VaR functions.

    # Setup:
    #if(exists(modified)({if( modified == TRUE) { method="modified" }}
    #if(method == TRUE or is.null(method) ) { method="modified" }

    # Forcing the following parameters
    method=c("modified","gaussian","historical", "kernel")[3]
    clean=c("none","boudt","geltner")
    portfolio_method=c("single","component","marginal")
    weights=NULL; mu=NULL; sigma=NULL; m3=NULL; m4=NULL; invert=TRUE
    clean = clean[1]
    portfolio_method = portfolio_method[1]

    myVaR = VaR(R = data , p=p, ..., method=method,
                       clean=clean,  portfolio_method=portfolio_method,
                       weights=weights, mu=mu, sigma=sigma, m3=m3, m4=m4, invert=invert)

    ## when VaR is computed using single and historical, compute standard error

    if(portfolio_method == "single" & is.null(weights) & method == "historical"){
    if(!is.null(data)){
      data <- checkData(data, method="xts", ...)
      columns=colnames(data)
      if (!is.null(weights) & portfolio_method != "single") {
        if ( length(weights) != ncol(data)) {
          stop("number of items in weights not equal to number of columns in data")
        }
      }
      # weights = checkData(weights, method="matrix", ...) #is this necessary?
      # TODO check for date overlap with data and weights
      if(clean!="none" & is.null(mu)){ # the assumption here is that if you've passed in any moments, we'll leave data alone
        data = as.matrix(PerformanceAnalytics::Return.clean(data, method=clean))
      }
      if(portfolio_method != "single"){
        # get the moments ready
        if (is.null(mu)) { mu =  apply(data,2,'mean' ) }
        if (is.null(sigma)) { sigma = cov(data) }
        if(method=="modified"){
          if (is.null(m3)) {m3 = M3.MM(data)}
          if (is.null(m4)) {m4 = M4.MM(data)}
        }
      }
    } else {
      #data is null, check for moments
      if(is.null(mu)) stop("Nothing to do! You must pass either data or the moments mu, sigma, etc.")
      if ( length(weights) != length(mu)) {
        stop("number of items in weights not equal to number of items in the mean vector")
      }
    }
    if(is.null(se.method)){
        return(myVaR)
      } else {
        res=list(VaR=myVaR)
        # for each of the method specified in se.method, compute the standard error
        for(mymethod in se.method){
          res[[mymethod]]=EstimatorSE(data, estimator.fun = "VaR", alpha=1-p,
                                      se.method = mymethod,
                                      cleanOutliers=cleanOutliers,
                                      fitting.method=fitting.method,
                                      ...)
        }
        return(res)
      }
    }
  } # end VaR wrapper function

