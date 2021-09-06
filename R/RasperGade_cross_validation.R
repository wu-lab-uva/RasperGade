#' @title  Leave-one-out cross-validation under pulsed evolution model
#'
#' @description Predict hidden states under the pulsed evolution model in a leave-one-out cross-validation
#'
#' @param FMR the returned data structure from function `fullMarginalReconstructionWithPE`
#' @param add.epsilon logical, if true, time-independent variation is added to the variance in the data frame
#' @param laplace logical, if true, Laplace distribution is used for time-independent variation
#' @param numApprox the number of normal distributions to approximate the Laplace distribution
#' @param margin the total probability mass that the number of jumps omitted in a compound Poisson process
#' @param numCores the number of cores to run in parallel
#' @param asymptotic the threshold of expected number of jumps on a branch beyond which normal distribution is assumed
#' @return `$summary` is a data frame listing the means and variances of the hidden states
#' @return `$error` is a list of error distributions where each element is a data frame
#' @export
#' @rdname LOO_CV_with_PE
LOO_CV_with_PE = function(FMR,add.epsilon=TRUE,laplace=FALSE,numApprox=1,margin=1e-6,numCores=1,asymptotic=5){
  cat("Conducting leave-one-out cross validation...\n")
  # the internal function to predict hidden states for each tip
  cv.func = function(tip,FMR){
    this.match = which(FMR$marginal$orientation==tip)
    error.index = do.call(expand.grid,lapply(1:length(this.match),function(i){
      1:dim(FMR$error[[this.match[i]]]$approximate)[1]
    }))
    rec = do.call(rbind,lapply(1:dim(error.index)[1],function(i){
      this.rec = marginalReconstructionWithPE(x = sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$x[error.index[i,j]]
      }),
      l = FMR$phy$edge.length[FMR$marginal$edge[this.match]],
      epsilons = sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$var[error.index[i,j]]
      }),
      sigma = FMR$params$sigma,lambdas=FMR$params$lambdas,sizes = FMR$params$sizes,margin = margin,asymptotic = asymptotic)
      this.probs = prod(sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$probs[error.index[i,j]]
      }))
      this.error = this.rec$profile
      this.error$probs = this.error$probs*this.probs
      return(this.error)
    }))
    rec$probs = rec$probs/sum(rec$probs)
    rec.moments = compoundNormalMoments(probs = rec$probs,means = rec$x,vars = rec$var)
    this.hsp = data.frame(node=tip,label = FMR$phy$tip.label[tip],
                          x = unname(rec.moments[1]),
                          var = unname(rec.moments[2]),
                          stringsAsFactors = FALSE)
    return(list(error=rec,hsp=this.hsp))
  }
  # run in parallel if required
  cv.error = mclapply(1:Ntip(FMR$phy),cv.func,FMR=FMR,mc.cores = numCores)
  # parse results
  cv = do.call(rbind,lapply(cv.error,function(x){x$hsp}))
  cv.error = lapply(cv.error,function(x){x$error})
  # add time-independent variation to error profiles in needed
  if(add.epsilon){
    cv$var = cv$var +FMR$params$epsilon/2
    if(numApprox>1&laplace){
      cv.error = lapply(cv.error,function(x){
        return(data.frame(x=c(x$x,x$x),
                          var=c(x$var+FMR$params$epsilon/2*exp(-1),x$var+FMR$params$epsilon/2*exp(1)),
                          probs=c(x$probs*exp(1)/(1+exp(1)),x$probs/(1+exp(1)))))
      })
    }else{
      for(tip in 1:Ntip(FMR$phy)){
        cv.error[[tip]]$var = cv.error[[tip]]$var+FMR$params$epsilon/2
      }
    }
  }
  #
  return(list(cv=cv,error=cv.error))
}

#' @title  Leave-one-out cross-validation under BM model
#'
#' @description Predict hidden states under the BM model in a leave-one-out cross-validation
#'
#' @param trait named vector of tip trait values; names should match the tip labels; missing values are hidden states
#' @param phy phylo-class object from `ape` package
#' @param numCores number of threads to run in parallel
#' @return `$summary` is a data frame listing the means and variances of the hidden states
#' @return `$error` is a list of error distributions where each element is a data frame
#' @export
#' @rdname LOO_CV_with_BM
LOO_CV_with_BM = function(phy,trait){
  cat("Conducting leave-one-out cross validation...\n")
  cv = lapply(1:Ntip(phy),function(i){predict_hidden_state_by_pic(phy=phy,trait = trait[phy$tip.label[-i]])})
  return(list(summary=do.call(rbind,lapply(cv,function(x){x$summary})),
              error=do.call(c,lapply(cv,function(x){x$error}))))
}
