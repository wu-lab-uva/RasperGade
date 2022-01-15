#' @title  Get confidence interval for ancestral or hidden states
#'
#' @description  Calculate the confidence interval from the calculated posterior probability
#'
#' @param error a data frame listing the mean, variance and weight of normal distributions
#' @param alpha the desired confidence level;default to 0.05
#' @param tol the precision of calculated CI that is passed to `uniroot` with the same default
#' @details `calculate_CI_from_error_distribution` calls `uniroot` function to find the CI at the desired confidence level.
#' @details For a mixture of normal distribution, the quantile at cumulative probability p Q_p will always be within the interval
#' @details [min(Q_p_i),max(Q_p_i)] where Q_p_i stands for the quantiles at cumulative probability p for each of the normal component.
#' @export
#' @rdname calculate_CI
calculate_CI_from_error_distribution = function(error,alpha=0.05,tol=.Machine$double.eps^0.25){
  lower.bound = qnorm(p = alpha/2,mean = error$x,sd = sqrt(error$var),lower.tail = TRUE)
  upper.bound = qnorm(p = alpha/2,mean = error$x,sd = sqrt(error$var),lower.tail = FALSE)
  eqa = function(q,lower.tail=TRUE){
    pMixNormal(q = q,lower.tail = lower.tail,
               mean = error$x,sd = sqrt(error$var),probs = error$probs)-alpha/2
  }
  lower.root = uniroot(f = eqa,lower.tail=TRUE,tol = tol,
                       lower = min(lower.bound)-0.1,upper = max(lower.bound)+0.1)$root
  upper.root = uniroot(f = eqa,lower.tail=FALSE,tol = tol,
                       lower = min(upper.bound)-0.1,upper = max(upper.bound)+0.1)$root
  return(c(lower=lower.root,upper=upper.root))
}