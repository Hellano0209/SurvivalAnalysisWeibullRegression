###########################
# MODELO DE RISCO WEIBULL #
# COM FRAÇÃO DE CURA      #
###########################

logit = function(x){
  
  y = (1+exp(-x))^(-1)
  return(y)
  
}

dlogit = function(x){
  
  aux1 = exp(-x)
  aux2 = (1+exp(-x))^2
  
  return(aux1/aux2)
  
}

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

#### FDP COM FRAÇÃO DE CURA ####

FDP_Pop = function(y, shape, scale, p0){
  
  split = logit(p0)
  pdf = split*dweibull(y, shape, scale)
  
  return(as.vector(pdf))
  
}

#### SOBREVIVÊNCIA COM FRAÇÃO DE CURA ####

SF_Pop = function(y, shape, scale, p0){
  
  split = logit(p0)
  sf = (1-split) + split*pweibull(y, shape, scale, lower.tail = F)
  
  return(as.vector(sf))
  
}

#### FUNÇÃO DE LOG-VEROSSIMILHANÇA ####

log.lik = function(par){
  
  shape   = par[1]
  scale   = par[2]
  p0      = par[3]
  
  loglik = sum( 
    cens*log(FDP_Pop(y, shape, scale, p0)) + (1-cens)*log(SF_Pop(y, shape, scale, p0))
  )
  
  return(loglik)
}

log.lik02 = function(par){
  
  shape   = par[1]
  scale   = par[2]
  p0      = par[3]
  
  loglik = sum( 
    cens*log(FDP_Pop(y, shape, scale, p0)) + (1-cens)*log(SF_Pop(y, shape, scale, p0))
  )
  
  return(-loglik)
}


#### FUNÇÃO SUMMARY ####

estimates = function(values){
  
  loglik      = -values$value
  split       = logit(values$par[3])
  estimativas = c(values$par[1:2], split)
  varerror    = diag(solve(values$hessian))
  varerror[3] = varerror[3]*dlogit(values$par[3])^2
  stderror    = sqrt(varerror)
  teste       = estimativas/stderror
  pvalue      = ifelse(teste > 0, 2*pnorm(-teste, lower.tail = T), 2*pnorm(teste, lower.tail = T))
  
  estimates = round(data.frame(estimativas, stderror, teste, pvalue), 4)
  rownames(estimates) = c("shape", "scale", "split")
  
  result    = list(estimates = estimates, loglik = loglik)
  
  return(result)
  
}
