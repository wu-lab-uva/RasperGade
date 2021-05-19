#' @title  Probability functions of mixture of normal distributions
#'
#' @description  Density function, cumulative distribution function and random generation for the mixture of normal distributions
#'
#' @param x the quantile for which density will be determined
#' @param q the quantile for which cumulative distribution will be determined
#' @param n the number of random values to return
#' @param mean the vector of means of mixed normal distributions
#' @param sd the vector of standard deviations of mixed normal distributions
#' @param probs the vector of normalized probability weight of mixed normal distributions
#' @param lower.tail logical, if TRUE (default), probability of X being smaller than or equal to q is returned. If false, the probability of being greater is returned.
#' @return `dMixNormal` returns the probability density
#' @return `pMixNormal` returns the cumulative distribution
#' @return `rMixNormal` returns a vector of random numbers following the specified distribution
#' @examples
#' # generate some random numbers from the mixture of two normal distributions
#' # with equal variance and weight but different means
#' norm.mix = rMixNorm(n=1e3,mean=c(-1,1),sd=c(1,1),probs=c(0.5,0.5))
#'
#' # test the empirical distribution against the theoretical one
#' ks.test(norm.mix,pMixNormal,mean=c(-1,1),sd=c(1,1),probs=c(0.5,0.5))
#' @export
#' @rdname MixNormal
dMixNormal = function(x,mean,sd,probs){
  dd = sapply(x,function(xx){
    sum(dnorm(x = xx,mean = mean,sd = sd)*probs)
  })
  return(dd)
}

#' @export
#' @rdname MixNormal
pMixNormal = function(q,mean,sd,probs,lower.tail=TRUE){
  pp = sapply(q,function(qq){
    sum(pnorm(q = qq,mean = mean,sd = sd,lower.tail = lower.tail)*probs)
  })
  return(pp)
}

#' @export
#' @rdname MixNormal
rMixNormal = function(n,mean,sd,probs){
  rr = sapply(1:n,function(i){
    this.x = rnorm(n = length(probs),mean = mean,sd = sd)
    if(length(this.x)>1){
      return(sample(x=this.x,size = 1,replace = TRUE,prob = probs))
    }else{
      return(this.x)
    }
  })
  return(rr)
}
