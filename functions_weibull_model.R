###########################
# MODELO DE RISCO WEIBULL #
###########################

#### FUNÇÃO DE RISCO BASE ####

HR_base = function(y, shape, scale){
  
  hr = dweibull(y, shape, scale)/pweibull(y, shape, scale, lower.tail = F)
  
  return(hr)
  
}

#### CUMULADA FUNÇÃO RISCO BASE ####

CHR_base = function(y, shape, scale){
  
  chr = -pweibull(y, shape, scale, lower.tail = F, log.p = T)
  
  return(chr)
  
}

#### FUNÇÃO DE LOG-VEROSSIMILHANÇA ####

log.lik = function(par){
  
  shape   = par[1]
  scale   = par[2]
  
  loglik = sum( 
    cens*log(HR_base(y, shape, scale)) + pweibull(y, shape = shape, scale = scale, lower.tail = F, log.p = T)
  )
  
  return(loglik)
}

log.lik02 = function(par){
  
  shape   = par[1]
  scale   = par[2]
  
  loglik = sum( 
    cens*log(HR_base(y, shape, scale)) + pweibull(y, shape = shape, scale = scale, lower.tail = F, log.p = T)
  )
  
  return(-loglik)
}

#### FUNÇÃO SUMMARY ####

estimates = function(fn, lower = c(-Inf, -Inf), upper = c(Inf, Inf)){
  
  # Chutes Iniciais
  
  fitweibull = fitdistr(y, "weibull")
  fitkaplan  = survfit(Surv(y, cens) ~ 1)
  guessshape = round(as.numeric(fitweibull$estimate[1]), 2)
  guessscale = round(as.numeric(fitweibull$estimate[2]), 2)
  par0       = c(guessshape, guessscale)
  
  values = optim(par = par0, fn, lower = lower, upper = upper, hessian = T, method = 'L-BFGS-B')
  
  loglik      = -values$value
  estimativas = values$par
  stderror    = sqrt(diag(solve(values$hessian)))
  teste       = estimativas/stderror
  
  estimates = data.frame(estimativas, stderror, teste)
  
  result    = list(estimates, loglik)
  
  return(result)
  
}
