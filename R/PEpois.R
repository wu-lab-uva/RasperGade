#' @title  Probability functions of trait change under pulsed evolution with 1 compound Poisson process
#'
#' @description  Density function, cumulative distribution function and random generation for the trait change under Pulsed evolution with 1 compound Poisson process
#'
#' @param x the trait change between two nodes
#' @param q the trait change between two nodes
#' @param n the number of random values to return
#' @param t the branch length between two nodes
#' @param lambda the rate at which jumps happen
#' @param size the variance of trait change introduced per jump
#' @param sigma the rate of BM (variance of gradual trait change per unit branch length)
#' @param epsilon the variance of time-independent variation
#' @param margin the total probability mass of the number of jumps not covered, defaults to 1e-6
#' @param lower.tail logical, if TRUE (default), probability of trait change smaller or equal to q is returned. If false, the probability of being greater is returned.
#' @return `dPEpoisnorm` returns the probability density
#' @return `pPEpoisnorm` returns the cumulative distribution
#' @return `rPEpoisnorm` returns a vector of random numbers following the specified distribution
#' @examples
#' # generate some random trait changes under pulsed evolution
#' dtrait = rPEpoisnorm(n=1e3,t=1,lambda=1,size=1,sigma=1,epsilon=0.1)
#'
#' # test the empirical distribution against the theoretical one
#' ks.test(dtrait,pPEpoisnorm,t=1,lambda=1,size=1,sigma=1,epsilon=0.1)
#' @export
#' @rdname PEpoisnorm
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


#' @export
#' @rdname PEpoisnorm
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
#' @rdname PEpoisnorm
rPEpoisnorm = function(n,t,lambda,size,sigma,epsilon=0){
  En = lambda*t
  k = rpois(n = n,lambda = En)
  sds = sqrt(k*size+sigma*t+epsilon)
  x = rnorm(n = n,mean = 0,sd = sds)
  return(x)
}
