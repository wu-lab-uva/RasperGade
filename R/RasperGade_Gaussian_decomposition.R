#' @title  Mix error distributions
#'
#' @description Mix error distributions by weights
#'
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

#' @title  Calculate error mean value
#'
#' @description Calculate error mean value
#'
#' @export
#' @rdname calculate_error_mean
calculate_error_mean = function(error){
  sum(error$x*error$probs)
}