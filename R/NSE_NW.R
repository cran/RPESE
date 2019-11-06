# Newey-West NSE estimators.
nse.nw <- function(x,prewhite = FALSE) {
  out = sandwich::lrvar(x = x, type = "Newey-West", prewhite = prewhite, adjust = TRUE)
  out = unname(out)
  return(out)
}
