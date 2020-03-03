#' @title Formatted Output for Standard Errors Functions in RPESE
#'
#' @description \code{printSE} returns a formatted output from standard error functions from RPESE.
#'
#' @param SE.data Standard error estimates output from RPESE functions.
#' @param round.digit Number of digits for rounding.
#' @param round.out Round data (TRUE) with round.digit number of digits. Default is TRUE.
#'
#' @return A data frame with formatted output from standard error functions from \code{RPESE}.
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
#' ES.out <- ES.SE(edhec, se.method=c("IFiid","IFcor"),
#'                 cleanOutliers=FALSE,
#'                 fitting.method=c("Exponential", "Gamma")[1])
#' # Print the output
#' printSE(ES.out)
#'
printSE <- function(SE.data, round.digit = 3, round.out = TRUE){
  N = length(SE.data)
  # if(N != 2) {
  #   cat("the SE.dataults do not contain standard errors!\n")
  #   return()
  # }
  list.names = names(SE.data)
  SE.data.df = data.frame(t(SE.data[[1]]))
  for(i in 2:length(list.names)){
    SE.data.df = cbind(SE.data.df, SE.data[[i]]$se)
  }
  colnames(SE.data.df) = list.names
  rownames(SE.data.df) = colnames(SE.data[[1]])
#  SE.data.df = round(SE.data.df, digits = round.digit)
  if(round.out){
    return(round(SE.data.df, digits = round.digit))
  }
  return(SE.data.df)
  # SE.data.df[2] = paste("(",SE.data.df[,2],")",sep="")
  # SE.data.df[,-1] = apply(as.data.frame(SE.data.df[,-1]),2,function(x) paste("(",x,")",sep=""))

}
