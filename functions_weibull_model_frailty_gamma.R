###########################
# MODELO DE RISCO WEIBULL #
# COM FRAGILIDADE GAMMA   #
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

#### FDP NÃO CONDICIONAL FRAGILIDADE ####

FDP_Fragilidade_Gamma = function(y, shape, scale, sigma2){
  
  theta = exp(sigma2)
  pdf = (HR_base(y, shape, scale))/((1 + theta*CHR_base(y, shape, scale))^(1+1/theta))
  
  return(as.vector(pdf))
  
}

#### SOBREVIVÊNCIA NÃO CONDICIONAL FRAGILIDADE GAMA ####

SF_Fragilidade_Gamma = function(y, shape, scale, sigma2){
  
  theta = exp(sigma2)
  sf = (1 + theta*CHR_base(y, shape, scale))^(-1/theta)
  
  return(as.vector(sf))
  
}

#### FUNÇÃO DE LOG-VEROSSIMILHANÇA ####

log.lik = function(par){
  
  shape   = par[1]
  scale   = par[2]
  sigma2  = par[3]
  
  loglik = sum( 
    cens*log(FDP_Fragilidade_Gamma(y, shape, scale, sigma2)) + 
      (1-cens)*log(SF_Fragilidade_Gamma(y, shape, scale, sigma2))
  )
  
  return(loglik)
}

log.lik02 = function(par){
  
  shape   = par[1]
  scale   = par[2]
  sigma2  = par[3]
  
  loglik = sum( 
    cens*log(FDP_Fragilidade_Gamma(y, shape, scale, sigma2)) + 
      (1-cens)*log(SF_Fragilidade_Gamma(y, shape, scale, sigma2))
  )
  
  return(-loglik)
}

#### SUMMARY ####

estimates = function(values){
  
  loglik      = -values$value
  frail       = exp(values$par[3])
  estimativas = c(values$par[1:2], frail)
  varerror    = diag(solve(values$hessian))
  varerror[3]  = varerror[3]*exp(values$par[3])^2
  stderror    = sqrt(varerror)
  teste       = estimativas/stderror
  pvalue      = ifelse(teste > 0, 2*pnorm(-teste, lower.tail = T), 2*pnorm(teste, lower.tail = T))
  
  result = round(data.frame(estimativas, stderror, teste, pvalue), 4)
  rownames(result) = c('shape', 'scale', 'frail')
  
  return(list(estimates = result, Loglik = loglik))
  
}