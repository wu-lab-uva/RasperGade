#' @title  Shift the center of deviations to a probability
#' @description  Shift the center of deviations to a probability for a balanced measurement 
#' @export
#' @rdname logitSigmoid
logitSigmoid=function(x,p){
  sapply(x,function(xx){
    if(is.na(xx)|is.nan(xx)|is.nan(p)|is.na(p)) return(NA)
    logit.x = logit(x = xx,a = 1,b = -logit(p))
    trans.x = sigmoid(logit.x)
  })
}
#' @export
#' @rdname logitSigmoid
reverse.logitSigmoid=function(x,p){
  sapply(x,function(xx){
    if(is.na(xx)|is.nan(xx)|is.nan(p)|is.na(p)) return(NA)
    logit.x = logit(x = xx,a = 1,b = logit(p))
    trans.x = sigmoid(logit.x)
  })
}

#' @title  Calculate pseudo-Z-score
#' @description  Calculate pseudo-Z-score that follows the standard normal distribution
#' @export
#' @rdname pseudoZscore
pseudo.Z.score = function(trait,pred,error,epsilon=0,laplace=FALSE){
  pZ = sapply(1:length(pred),function(i){
    if(is.list(error)){
      if(laplace){
        pp = pMixNormalLaplace(q = trait[i],mean = error[[i]]$x,sd = sqrt(error[[i]]$var),
                               probs = error[[i]]$probs,laplace.var = epsilon/2)
      }else{
        pp = pMixNormal(q = trait[i],mean = error[[i]]$x,
                        sd = sqrt(error[[i]]$var+epsilon/2),probs = error[[i]]$probs)
      }
    }else{
      if(laplace){
        pp = pconv.norm.laplace(q = trait[i]-pred[i],sigma = error[i],epsilon = epsilon/2)
      }else{
        pp = pnorm(q = trait[i],mean = pred[i],sd = sqrt(error[i]+epsilon/2))
      }
    }
    return(pp)
  })
  return(list(p=pZ,q=qnorm(pZ)))
}

#' @title  Discretize result
#' @description  Extending RasperGade to integer traits
#' @export
#' @rdname discretizeResult
discretizeResult = function(res,error=NULL,epsilon=0,laplace=FALSE,numCores=1){
  new.res = mclapply(1:dim(res)[1],function(i){
    df = data.frame(node=res$node[i],label=res$label[i],x=round(res$x[i]))
    if(laplace){
      if(is.null(error)){
        df$probs = sapply(res$var,function(x){pconv.norm.laplace(q = 0.5,sigma = x,epsilon = epsilon/2)})-
          sapply(error,function(x){pconv.norm.laplace(q = -0.5,sigma = x,epsilon = epsilon/2)})
      }else{
        df$probs = pMixNormalLaplace(q = round(res$x[i])+0.5,mean = error[[i]]$x,
                                     sd = sqrt(error[[i]]$var),laplace.var = epsilon/2,probs = error[[i]]$probs)-
          pMixNormalLaplace(q = round(res$x[i])-0.5,mean = error[[i]]$x,
                            sd = sqrt(error[[i]]$var),laplace.var = epsilon/2,probs = error[[i]]$probs)
      }
    }else{
      if(is.null(error)){
        df$probs = pnorm(q = round(res$x[i])+0.5,mean = res$x[i],sd = sqrt(res$var[i]+epsilon/2))-
          pnorm(q = round(res$x[i])-0.5,mean = res$x[i],sd = sqrt(res$var[i]+epsilon/2))
      }else{
        df$probs = pMixNormal(q = round(res$x[i])+0.5,mean = error[[i]]$x,
                              sd = sqrt(error[[i]]$var+epsilon/2),probs = error[[i]]$probs)-
          pMixNormal(q = round(res$x[i])-0.5,mean = error[[i]]$x,
                     sd = sqrt(error[[i]]$var+epsilon/2),probs = error[[i]]$probs)
      }
    }
    return(df)
  },mc.cores=numCores)
  new.res = do.call(rbind,new.res)
  return(new.res)
}
#' @export
calculate_R2cv = function(trait,pred){
  trait.variance = mean((trait-mean(trait))^2)
  return(1-mean((trait - pred)^2)/trait.variance)
}
#' @export
F_beta_score = function(obs,pred,beta=1){
  tp = sum(obs&pred)
  fp = sum((!obs)&pred)
  fn = sum(obs&(!pred))
  fbeta = ((1+beta^2)*tp)/((1+beta^2)*tp+beta^2*fn+fp)
  return(fbeta)
}
#' @export
precision_recall = function(obs,pred,as.na=FALSE){
  tp = sum(obs&pred)
  fp = sum((!obs)&pred)
  fn = sum(obs&(!pred))
  precision=tp/(tp+fp)
  recall = tp/(tp+fn)
  if(as.na){
    if(is.nan(precision)) precision = NA
    if(is.nan(recall)) recall=NA
  }else{
    if(is.nan(precision)) precision = 0
    if(is.nan(recall)) recall=0
  }
  return(c(precision=precision,recall=recall))
}
