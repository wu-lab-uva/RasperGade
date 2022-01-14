#' @title  Get confidence interval for ancestral or hidden states
#'
#' @description  Calculate the confidence interval from the calculated posterior probability
#'
#' @param error a data frame listing the mean, variance and weight of normal distributions
#' @param alpha the desired confidence level
#' @export
#' @rdname calculate_CI
calculate_CI_from_error_distribution = function(error,alpha=0.05){
  lower.bound = qnorm(p = alpha/2,mean = error$x,sd = sqrt(error$var),lower.tail = TRUE)
  upper.bound = qnorm(p = alpha/2,mean = error$x,sd = sqrt(error$var),lower.tail = FALSE)
  eqa = function(q,lower.tail=TRUE){
    pMixNormal(q = q,lower.tail = lower.tail,
               mean = error$x,sd = sqrt(error$var),probs = error$probs)
  }
  lower.root = uniroot(f = eqa,interval = 1)
  return(c(lower=lower.root,upper=upper.root))
}