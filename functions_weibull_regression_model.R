###########################
# MODELO DE RISCO WEIBULL #
###########################

#### FUNÇÃO DE RISCO BASE ####

HR_base = function(y, X, beta, shape){
  
  scale = exp(X%*%beta)
  hr = dweibull(y, shape, scale)/pweibull(y, shape, scale, lower.tail = F)
  
  return(hr)
  
}

#### CUMULADA FUNÇÃO RISCO BASE ####

CHR_base = function(y, X, beta, shape){
  
  scale = exp(X%*%beta)
  chr = -pweibull(y, shape, scale, lower.tail = F, log.p = T)
  
  return(chr)
  
}

#### FUNÇÃO DE LOG-VEROSSIMILHANÇA ####

log.lik = function(par){
  
  beta  = par[1:ncol(X)]
  shape = par[ncol(X)+1]
  scale = exp(X%*%beta)
  
  loglik = sum( 
    cens*log(HR_base(y, X, beta, shape)) + pweibull(y, shape = shape, scale = scale, lower.tail = F, log.p = T)
  )
  
  return(loglik)
}

log.lik02 = function(par){
  
  beta  = par[1:ncol(X)]
  shape = par[ncol(X)+1]
  scale = exp(X%*%beta)
  
  loglik = sum( 
    cens*log(HR_base(y, X, beta, shape)) + pweibull(y, shape = shape, scale = scale, lower.tail = F, log.p = T)
  )
  
  return(-loglik)
}

#### FUNÇÃO SUMMARY ####

estimates = function(values){
  
  loglik      = -values$value
  estimativas = values$par
  stderror    = sqrt(diag(solve(values$hessian)))
  teste       = estimativas/stderror
  pvalue      = ifelse(teste > 0, 2*pnorm(-teste, lower.tail = T), 2*pnorm(teste, lower.tail = T))
  
  estimates = round(data.frame(estimativas, stderror, teste, pvalue), 4)
  
  result    = list(estimates = estimates, loglik = loglik)
  
  return(result)
  
}
