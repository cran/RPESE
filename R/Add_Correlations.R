#'
#' @import RPEIF
#' @importFrom stats acf
#'
Add_Correlations <- function(res, data, cleanOutliers, corOut, IF.func, ...){

  # Robust filtering
  if(cleanOutliers)
    data <- apply(data, 2, function(x, ...) return(RPEIF:::robust.cleaning(x, ...)))

  # Computing the correlations for corOut
  for(cor.id in corOut){
    if(cor.id=="retCor"){
      res$retCor$out <- apply(data, 2, function(x){return((acf(x, plot=F)$acf[2]))})
    }
    if(cor.id=="retIFCor"){
      IF.ret <- apply(data, 2, function(x, prewhiten, cleanOutliers, ...)
        return(IF.func(x, prewhiten=prewhiten, cleanOutliers=cleanOutliers, ...)),
        prewhiten=FALSE, cleanOutliers=cleanOutliers)
      res$retIFCor$out <- apply(IF.ret, 2, function(x){return((acf(x, plot=F)$acf[2]))})
    }
    if(cor.id=="retIFCorPW"){
      IF.retPW <- apply(data, 2, function(x, prewhiten, cleanOutliers, ...)
        return(IF.func(x, prewhiten=prewhiten, cleanOutliers=cleanOutliers, ...)),
        prewhiten=TRUE, cleanOutliers=cleanOutliers)
      res$retIFCorPW$out <- apply(IF.retPW, 2, function(x){return((acf(x, plot=F)$acf[2]))})
    }
  }

  # Return the output
  return(res)
}


