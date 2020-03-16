#'
#' @import RPEIF
#'
#' @title Standard Error Estimate for Expected Shortfall (ES) of Returns
#'
#' @description \code{ES.SE} computes the standard error of the expected shortfall of the returns.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param p Confidence level for calculation. Default value is p=0.95.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One or a combination of:
#' \code{"IFiid"} (default), \code{"IFcor"} (default), \code{"IFcorPW"}, \code{"IFcorAdapt"},
#' \code{"BOOTiid"} or \code{"BOOTcor"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
#' @param d.GLM.EN Order of the polynomial for the Exponential or Gamma fitting. Default polynomial order of 5.
#' @param freq.include Frequency domain inclusion criteria. Must be one of "All" (default), "Decimate" or "Truncate."
#' @param freq.par Percentage of the frequency used if \code{"freq.include"} is "Decimate" or "Truncate." Default is 0.5.
#' @param corOut Return correlation of the returns or the influence function transformed returns. Must be one of "retCor", "retIFCor" or "none" (default).
#' @param ... Additional parameters.
#'
#' @return A vector or a list depending on \code{se.method}.
#'
#' @export
#'
#' @note The option to \code{invert} the ES measure should appease both
#' academics and practitioners.  The mathematical definition of ES as the
#' negative value of extreme losses will (usually) produce a positive number.
#' Practitioners will argue that ES denotes a loss, and should be internally
#' consistent with the quantile (a negative number).  For tables and charts,
#' different preferences may apply for clarity and compactness.  As such, we
#' provide the option, and set the default to TRUE to keep the return
#' consistent with prior versions of PerformanceAnalytics, but make no value
#' judgement on which approach is preferable.
#'
#' @section Background: This function provides several estimation methods for
#' the Expected Shortfall (ES) (also called Expected Tail Loss (ETL)
#' or Conditional Value at Risk (CVaR)) of a return series and the Component ES
#' (ETL/CVaR) of a portfolio.
#'
#' At a preset probability level denoted \eqn{c}, which typically is between 1
#' and 5 per cent, the ES of a return series is the negative value of the
#' expected value of the return when the return is less than its
#' \eqn{c}-quantile.  Unlike value-at-risk, conditional value-at-risk has all
#' the properties a risk measure should have to be coherent and is a convex
#' function of the portfolio weights (Pflug, 2000).  With a sufficiently large
#' data set, you may choose to estimate ES with the sample average of all
#' returns that are below the \eqn{c} empirical quantile. More efficient
#' estimates of VaR are obtained if a (correct) assumption is made on the
#' return distribution, such as the normal distribution. If your return series
#' is skewed and/or has excess kurtosis, Cornish-Fisher estimates of ES can be
#' more appropriate. For the ES of a portfolio, it is also of interest to
#' decompose total portfolio ES into the risk contributions of each of the
#' portfolio components. For the above mentioned ES estimators, such a
#' decomposition is possible in a financially meaningful way.
#'
#' @author Xin Chen, \email{chenx26@uw.edu}
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @references Boudt, Kris, Peterson, Brian, and Christophe Croux. 2008.
#' Estimation and decomposition of downside risk for portfolios with non-normal
#' returns. 2008. The Journal of Risk, vol. 11, 79-103.
#'
#' Cont, Rama, Deguest, Romain and Giacomo Scandolo. Robustness and sensitivity
#' analysis of risk measurement procedures. Financial Engineering Report No.
#' 2007-06, Columbia University Center for Financial Engineering.
#'
#' Laurent Favre and Jose-Antonio Galeano. Mean-Modified Value-at-Risk
#' Optimization with Hedge Funds. Journal of Alternative Investment, Fall 2002,
#' v 5.
#'
#' Martellini, Lionel, and Volker Ziemann.  Improved Forecasts of Higher-Order
#' Comoments and Implications for Portfolio Selection. 2007. EDHEC Risk and
#' Asset Management Research Centre working paper.
#'
#' Pflug, G. Ch.  Some remarks on the value-at-risk and the conditional
#' value-at-risk. In S. Uryasev, ed., Probabilistic Constrained Optimization:
#' Methodology and Applications, Dordrecht: Kluwer, 2000, 272-281.
#'
#' Scaillet, Olivier. Nonparametric estimation and sensitivity analysis of
#' expected shortfall. Mathematical Finance, 2002, vol. 14, 74-86.
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
#' ES.SE(edhec, se.method=c("IFiid","IFcor"),
#'       cleanOutliers=FALSE,
#'       fitting.method=c("Exponential", "Gamma")[1])
#'
ES.SE <- function(data, p=0.95,
                  se.method=c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")[1:2],
                  cleanOutliers=FALSE, fitting.method=c("Exponential", "Gamma")[1], d.GLM.EN=5,
                  freq.include=c("All", "Decimate", "Truncate")[1], freq.par=0.5,
                  corOut = c("none", "retCor","retIFCor", "retIFCorPW")[1],
                  ...)
{ # @author Brian G. Peterson and Xin Chen

  # We force the following parameters
  method=c("historical","gaussian","modified")[1]
  clean=c("none","boudt", "geltner")[1]
  portfolio_method=c("single","component")[1]
  weights=NULL; mu=NULL; sigma=NULL; m3=NULL; m4=NULL
  invert=TRUE; operational=TRUE
  # Descripion:

  # wrapper for univariate and multivariate ES functions.

  # Setup:
  #if(exists(modified)({if( modified == TRUE) { method="modified" }}
  #if(method == TRUE or is.null(method) ) { method="modified" }
  myES=ES(R=data , p=p, ...,
          method=method,
          clean=clean,
          portfolio_method=portfolio_method,
          weights=weights, mu=mu, sigma=sigma, m3=m3, m4=m4,
          invert=invert, operational=operational)

  # se.method=se.method[1]
  # SE=NULL
  method = method[1]
  clean = clean[1]
  portfolio_method = portfolio_method[1]
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
          if (is.null(m3)) {m3 = M3.MM(data,mu=mu)}
          if (is.null(m4)) {m4 = M4.MM(data,mu=mu)}
        }
      }
    } else {
      #data is null, check for moments
      if(is.null(mu)) stop("Nothing to do! You must pass either data or the moments mu, sigma, etc.")
      if (length(weights) != length(mu)) {
        stop("number of items in weights not equal to number of items in the mean vector")
      }
    }
    if(is.null(se.method)){
      return(myES)
    } else {
      res=list(ES=myES)
      # for each of the method specified in se.method, compute the standard error
      for(mymethod in se.method){
        res[[mymethod]]=EstimatorSE(data, estimator.fun = "ES", alpha.ES = 1-p,
                                    se.method=mymethod,
                                    cleanOutliers=cleanOutliers,
                                    fitting.method=fitting.method, d.GLM.EN=d.GLM.EN,
                                    freq.include=freq.include, freq.par=freq.par,
                                    ...)
      }

      # Adding the correlations to the list
      res <- Add_Correlations(res=res, data=data, cleanOutliers=cleanOutliers, corOut=corOut, IF.func=IF.ES, ...)

      # Returning the output
      return(res)
    }
  }
} # end ES.SE wrapper function

