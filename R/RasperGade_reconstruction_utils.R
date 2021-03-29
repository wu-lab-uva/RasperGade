#' @title  Calculate the product of normal probability density functions
#'
#' @description  Calculate the mean, variance and the scaling factor of the product of normal pdfs (a Gaussian function)
#'
#' @param mean the vector of means of mixed normal distributions
#' @param var the vector of standard deviations of mixed normal distributions
#' @details The function is implemented based on \href{http://www.lucamartino.altervista.org/2003-003.pdf}{Bromiley 2014}
#' @return A numeric vector of length 3 (mean, variance and scale)
#' @export
productNormalDensity = function(mean,var){
  if(length(mean)!=length(var)) stop("Mean and variance for every normal density function should be provided.")
  sig2 = 1/sum(1/var)
  mu = sum(mean/var)*sig2
  S = sqrt(2*pi*sig2)/prod(sqrt(2*pi*var))*exp(-0.5*(sum(mean^2/var)-mu^2/sig2))
  return(c(mean=mu,var=sig2,scale = S))
}

#' @title  Calculate four central moments of a mixture of normal distributions
#'
#' @description  Calculate the mean, variance, skewness and excess kurtosis of the mixture of normal distributions
#'
#' @param means the vector of means of mixed normal distributions
#' @param vars the vector of standard deviations of mixed normal distributions
#' @param probs the vector of probability weight of normal distributions
#' @return A numeric vector of length 4 (mean, variance, skewness and kurtosis)
#' @export
compoundNormalMoments = function(probs,means,vars){
  probs = probs/sum(probs)
  compound.mean = sum(probs*means)
  compound.var = sum(probs*vars)+sum(probs*(means^2))-compound.mean^2
  compound.skewness = (sum(probs*(means-compound.mean)^3)+sum(3*probs*(means-compound.mean)*vars))/
    sqrt(compound.var)^3
  compound.kurtosis=(sum(probs*(means-compound.mean)^4)+
                       sum(6*probs*(means-compound.mean)^2*vars)+
                       sum(probs*3*vars^2))/
    compound.var^2-3
  return(c(mean=compound.mean,var=compound.var,skewness = compound.skewness,kurtosis = compound.kurtosis))
}
