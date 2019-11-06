# Compute the periodogram as defined in H&W 1981
myperiodogram <- function(data,..., max.freq=0.5, twosided = FALSE, keep = 1){
  ## data.fft=myfft(data) This is very slow
  data.fft=fft(data)
  N=length(data)
  #   tmp=1:N
  #   inset=tmp[(1:N)<floor(N/2)]
  tmp = Mod(data.fft[2:floor(N/2)])^2/N
  tmp = sapply(tmp, function(x) max(0.00001,x))
  freq = ((1:(floor(N/2)-1))/N)
  tmp = tmp[1:floor(length(tmp) * keep)]
  freq = freq[1:floor(length(freq) * keep)]
  if (twosided){
    tmp = c(rev(tmp), tmp)
    freq = c(-rev(freq), freq)
  }
  return(list(spec=tmp,freq=freq))
}

SE.glmnet_exp <- function(data, ...,
                          d=7, alpha.EN=0.5,
                          keep=1,
                          standardize = FALSE,
                          return.coeffs = FALSE,
                          prewhiten = FALSE, ar.coeffs=NULL,
                          fitting.method = c("Exponential", "Gamma")[1]){
  ##### perform prewhitening

  N=length(data)
  # Step 1: compute the periodograms
  my.periodogram=myperiodogram(data, ...,  keep = keep)
  my.freq=my.periodogram$freq
  my.periodogram=my.periodogram$spec

  # remove values of frequency 0 as it does not contain information about the variance
  # my.freq=my.freq[-1]
  # my.periodogram=my.periodogram[-1]

  # implement cut-off
  nfreq=length(my.freq)
  # my.freq=my.freq[1:floor(nfreq*keep)]
  # my.periodogram=my.periodogram[1:floor(nfreq*keep)]

  # standardize

  # create 1, x, x^2, ..., x^d
  x.mat=rep(1,length(my.freq))
  for(col.iter in 1:d){
    x.mat=cbind(x.mat,my.freq^col.iter)
  }

  # b0 = rnorm(d + 1)

  # standardize x.mat
  mean_vec = apply(x.mat[,-1], 2, mean)
  sd_vec = apply(x.mat[,-1], 2, sd)
  if (standardize){
  for(i in 2:ncol(x.mat)){
    tmp = x.mat[,i]
    x.mat[,i] = (tmp - mean(tmp)) / sd(tmp)
  }
  }

  # Testing whether the fitting.method vector is valid
  if(!(fitting.method %in% c("Exponential", "Gamma")))
    stop("The specified fitting method is not available.")

  # fit the glmnet Gamma model
  if(fitting.method=="Gamma")
    res = RPEGLMEN::fit.glmGammaNet(x.mat, my.periodogram, alpha.EN = alpha.EN, ...) else if(fitting.method=="Exponential")
      res = RPEGLMEN::glmnet_exp(x.mat, my.periodogram, ..., alpha.EN = alpha.EN)


  # Step 3: return the estimated variance, and coeffs if return.coeffs = TRUE
  if(return.coeffs){
    if(standardize){
      variance = exp(sum(res * c(1, -mean_vec / sd_vec)))/N
      if(prewhiten)
        variance = variance / ( 1 - sum(ar.coeffs))^2
      coeffs = res
      return(list(variance, coeffs))
    }
    variance = exp(res[1])/N
    if(prewhiten)
      variance = variance / ( 1 - sum(ar.coeffs))^2
    coeffs = res
    return(list(variance, coeffs))
  }


  if (standardize){
    variance = exp(sum(res * c(1, -mean_vec / sd_vec)))/N
    if(prewhiten)
      variance = variance / ( 1 - sum(ar.coeffs))^2
    return(variance)
  }
  variance = exp(res[1])/N
  if(prewhiten)
    variance = variance / ( 1 - sum(ar.coeffs))^2
  return(variance)
}




