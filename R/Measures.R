################### This file includes all the functions to compute nonparametric sample estimators of risk/performance measures

# Compute sample Sharpe Ratio from data
SR=function(data,...,rf=0){
  mu.hat=mean(data)
  sigma.hat=sd(data)
  return((mu.hat-rf)/sigma.hat)
}


# Compute sample value-at-risk
VaR.hist=function(data,...,alpha=0.05){
  return(-quantile(data,alpha))
}

# Compute sample Expected Shortfall
ES.hist=function(data,...,alpha.ES=0.05){
  return(-mean(data[data<=quantile(data,alpha.ES)]))
}

# Compute sample Robust Expected Shortfall
RES=function(data,...,alpha=0.01,beta=0.1){
  return(-mean(data[(data<=quantile(data,beta))&(data>=quantile(data,alpha))]))
}

# Compute sample Sortino Ratio with mean threshold
SoR = function(data, ..., rf = 0, MAR = 0, threshold=c("mean", "const")[1]){

  if(threshold=="mean"){
    mu.hat = mean(data)
    sigma.minus.hat = sqrt(mean((data-mu.hat)^2*(data<=mu.hat)))
    SoR.hat = (mu.hat-rf)/sigma.minus.hat
    return(SoR.hat)
  } else if (threshold=="const"){
    mu.hat = mean(data)
    sigma.minus.hat = sqrt(mean((data-MAR)^2*(data<=MAR)))
    SoR.const.hat = (mu.hat-MAR)/sigma.minus.hat
    return(SoR.const.hat)
  }
}

# Compute sample Sortino Ratio with mean threshold
SoR.mean = function(data, ..., rf = 0){
  mu.hat = mean(data)
  sigma.minus.hat = sqrt(mean((data-mu.hat)^2*(data<=mu.hat)))
  SoR.hat = (mu.hat-rf)/sigma.minus.hat
  return(SoR.hat)
}

# Compute sample Sortino Ratio with constant threshold
SoR.const = function(data, ..., MAR = 0){
  mu.hat = mean(data)
  sigma.minus.hat = sqrt(mean((data-MAR)^2*(data<=MAR)))
  SoR.const.hat = (mu.hat-MAR)/sigma.minus.hat
  return(SoR.const.hat)
}


# Compute sample ESratio
ESratio = function(data, ..., alpha = 0.05, rf = 0){
  mu.hat=mean(data)
  return((mu.hat - rf)/ES.hist(data, alpha = alpha))
}

# Compute sample VaRratio
VaRratio = function(data, ..., alpha = 0.05, rf = 0){
  mu.hat=mean(data)
  return((mu.hat - rf)/VaR.hist(data, alpha = alpha))
}

# Compute sample lower partial moment
LPM = function(data, ..., const = 0, k = 1){
  N = length(data)
  return(1/N*sum((const-data[data<=const])^k)^(1/k))
}

# Compute sample Omega Ratio
OmegaRatio = function(data, ..., const = 0){
  N = length(data)
  OmegaPlus = sum(data[data>=const]-const)/N
  OmegaMinus = sum(const-data[data<=const])/N
  return(OmegaPlus/OmegaMinus)
}

# Compute Semi-Standard Deviation
SemiSD = function(data, ..., rf=0){

    mySSD = PerformanceAnalytics::SemiDeviation(data-rf)
    return(mySSD)
}

# Compute Rachev Ratio
RachevRatio = function(data, ..., alpha=0.1, beta=0.1, rf=0){

  # Computing the mean of the data
  mu.hat <- mean(data)
  # Computing the SD of the data
  sigma.hat <- mean((data-mu.hat)^2)

  # Computing the VaR of the data (lower tail)
  VaR.hat.lower <- -quantile(data, alpha)
  # Storing the negative value of the VaR based on the desired alpha (lower tail)
  quantile.lower <- -VaR.hat.lower
  # Computing the ES of the data (lower tail)
  ES.lower <- -mean(data[data<=-VaR.hat.lower])

  # Computing the VaR of the data (upper tail)
  n.upper <- floor((1-beta)*length(data))
  sorted.returns <- sort(as.numeric(data))
  VaR.hat.upper <- sorted.returns[n.upper]
  # Storing the negative value of the VaR based on the desired alpha (upper tail)
  quantile.upper <- VaR.hat.upper
  # Computing the ES of the data (upper tail)
  ES.upper <- mean(data[data>=quantile.upper])

  # Computing Rachev Ratio
  myRachevRatio <- ES.upper/ES.lower

  # Return Rachev Ratio
  return(myRachevRatio)
}
