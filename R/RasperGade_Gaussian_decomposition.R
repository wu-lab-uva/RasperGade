#' @title  Mix error distributions
#'
#' @description Mix error distributions by their weights
#'
#' @param errors a list of error distributions (a list of data frames)
#' @param weights a numerical vector of probability weights
#' @return The merged error distribution (a data frame)
#' @export
#' @rdname mix_errors_by_weight
mix_errors_by_weight = function(errors,weights){
  weights = weights/sum(weights)
  new.error = do.call(rbind,lapply(1:length(errors),function(i){
    this.error = errors[[i]]
    this.error$probs = this.error$probs*weights[i]
    return(this.error)
  }))
  return(new.error)
}

#' @title Calculate the mean and variance of an error distribution
#'
#' @description The function [calculate_error_mean_and_var] calculates the mean and variance of an error distribution
#'
#' @param error an error distributions (a data frame)
#' @return A named vector of two containing the mean and variance of the error distribution
#' @export
#' @rdname calculate_error_mean_and_var
calculate_error_mean_and_var = function(error){
  this.moments = unname(compoundNormalMoments(probs = error$probs,means = error$x,vars = error$var))
  return(this.moments[1:2])
}