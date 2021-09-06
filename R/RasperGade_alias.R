#' @export
#' @rdname LOO_CV_with_PE
crossValidationWithPE = function(FMR,add.epsilon=TRUE,laplace=FALSE,numApprox=1,margin=1e-6,numCores=1,asymptotic=5){
  warning("crossValidationWithPE is replaced by LOO_CV_with_PE and will be removed in future updates.",immediate. = TRUE)
  LOO_CV_with_PE(FMR,add.epsilon,laplace,numApprox,margin,numCores,asymptotic)
}

#' @export
#' @rdname predict_hidden_state_by_pic
predTraitByPIC = function(phy, trait){
  warning("predTraitByPIC is replaced by predict_hidden_state_by_pic and will be removed in future updates.",immediate. = TRUE)
  predict_hidden_state_by_pic(phy,trait)
}
