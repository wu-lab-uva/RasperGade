#' Probability density function of trait change under pulsed evolution
#'
#' This function returns the probability density of having a trait change x given branch length t
#'
#' @param x the trait change between two nodes
#' @param t the branch length bewteen two nodes
#' @param lambda the rate at which jumps happen
#' @param size the variance of trait change introduced per jump
#' @param sigma the rate of BM (variance of gradual trait change per unit branch length)
#' @param epsilon the variance of time-independent variation
#' @param margin the total probability mass of the number of jumps not covered, defaults to 1e-6
#' @return The numeric value of probability density
#' @seealso [pPEpoisnorm()]
#' @export dPEpoisnorm
#' @examples
#' dPEpoisnorm(x=0,t=1,lambda=1,size=1,sigma=1,epsilon=0,margin=1e-6)
dPEpoisnorm = function(x,t,lambda,size,sigma,epsilon=0,margin=1e-6){
  En = lambda*t
  if(En>1e4){
    p = dnorm(x = x,mean = 0,sd = sqrt(size*En+sigma*t+epsilon))
  }else{
    k = qpois(p = c(margin/2,1-margin/2),lambda = En)
    k= floor(k[1]):ceiling(k[2])
    probs = dpois(x = k,lambda = En)
    p = sum(dnorm(x = x,mean = 0,sd = sqrt(size*k+sigma*t+epsilon))*probs)
  }
  return(p)
}

#' Cumulative distribution function of trait change under pulsed evolution
#'
#' This function returns the cumulative probability of having a trait change q given branch length t
#'
#' @param q the trait change between two nodes
#' @param t the branch length bewteen two nodes
#' @param lambda the rate at which jumps happen
#' @param size the variance of trait change introduced per jump
#' @param sigma the rate of BM (variance of gradual trait change per unit branch length)
#' @param epsilon the variance of time-independent variation
#' @param margin the total probability mass of the number of jumps not covered, defaults to 1e-6
#' @param lower.tail whether
#' @return The numeric value of probability density
#' @seealso [dPEpoisnorm()]
#' @export
#' @examples
#' pPEpoisnorm(q=0,t=1,lambda=1,size=1,sigma=1,epsilon=0)
pPEpoisnorm = function(q,t,lambda,size,sigma,epsilon=0,margin=1e-6,lower.tail=TRUE){
  En = lambda*t
  if(En>1e4){
    p = pnorm(q = q,mean = 0,sd = sqrt(size*En+sigma*t+epsilon),lower.tail = lower.tail)
  }else{
    k = qpois(p = c(margin/2,1-margin/2),lambda = En)
    k= floor(k[1]):ceiling(k[2])
    probs = dpois(x = k,lambda = En)
    p = sum(pnorm(q = q,mean = 0,sd = sqrt(size*k+sigma*t+epsilon),lower.tail = lower.tail)*probs)
  }
  return(p)
}

#' @export
rPEpoisnorm = function(n,t,lambda,size,sigma,epsilon=0){
  En = lambda*t
  k = rpois(n = n,lambda = En)
  sds = sqrt(k*size+sigma*t+epsilon)
  x = rnorm(n = n,mean = 0,sd = sds)
  return(x)
}

#' @export
dMixNormal = function(x,mean,sd,probs){
  dd = sapply(x,function(xx){
    sum(dnorm(x = xx,mean = mean,sd = sd)*probs)
  })
  return(dd)
}

#' @export
pMixNormal = function(q,mean,sd,probs){
  pp = sapply(q,function(qq){
    sum(pnorm(q = qq,mean = mean,sd = sd)*probs)
  })
  return(pp)
}

#' @export
rMixNormal = function(n,mean,sd,probs){
  rr = sapply(1:n,function(i){
    sample(x=rnorm(n = length(probs),mean = mean,sd = sd),size = 1,replace = TRUE,prob = probs)
  })
  return(rr)
}
