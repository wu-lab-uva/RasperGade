#' @import pracma

#' @title  Calculate RSS
#' @export
calculateBinomRSS = function(p,obs,group=NULL,bin=20,equal.bin=TRUE,detail=FALSE,...){
  if(is.null(group)){
    mid.bin = seq(0,1,length.out=bin+1)
    mid.bin = (mid.bin[-1]+mid.bin[-(bin+1)])/2
    all.bins = lapply(mid.bin,function(z){
      which((p<=(z+mid.bin[1]))&(p>(z-mid.bin[1])))
    })
  }else{
    if(equal.bin){
      group.idx = sort(group,decreasing = FALSE,index.return=TRUE)$ix
      bin.idx = round(seq(from=1,to = length(group)+1,length.out = bin+1))
      all.bins = lappply(1:bin,function(z){group.idx[bin.idx[z]:(bin.idx[z+1]-1)]})
    }else{
      mid.bin = seq(min(group[is.finite(group)]),max(group[is.finite(group)])+1e-6,length.out=bin+1)
      mid.bin = (mid.bin[-1]+mid.bin[-(bin+1)])/2
      all.bins = lapply(mid.bin,function(z){
        which((group<(z+mid.bin[1]))&(group>=(z-mid.bin[1])))
      })
    }
  }
  bin.rss = do.call(rbind,lapply(1:bin,function(z){
    idx = all.bins[[z]]
    total.count = length(idx)
    if(total.count<=0) return(data.frame(RSS=NA,emp=NA,exp=NA,count=total.count))
    exp.z = mean(p[idx])
    obs.count = sum(obs[idx])
    rss = (2*logitSigmoid(x = obs.count/total.count,p = exp.z)-1)^2
    return(data.frame(RSS=rss,emp=obs.count/total.count,exp=exp.z,count=total.count))
  }))
  this.rss = bin*mean(bin.rss,na$RSS.rm = TRUE)
  if(detail) return(list(RSS=this.rss,detail=bin.rss))
  return(this.rss)
}

#' @title  Calculate median Chi-square
#' @export
calculateHeteroscedasticity = function(x,y,bin=100,n=1000){
  idx = sort(y,index.return=TRUE)$ix
  bins=round(seq(0,length(x),length.out = bin+1))
  vars = lapply(1:bin,function(i){
    x[idx[(bins[i]+1):bins[i+1]]]
  })
  h = fligner.test(x = vars)
  return(c(h=h$statistic,p.value=h$p.value))
}

#' @title  Analyze quality of uncertainty in ancestral or hidden state prediction
#'
#' @description  Analyze quality of uncertainty in ancestral or hidden state prediction from different perspectives
#'
#' @param trait true trait values
#' @param pred point estimates of ancestral or hidden states
#' @param error uncertainty of predictions
#' @param epsilon variance of time-independent variation
#' @param laplace logical, if true, time-independent variation follows a Laplace distribution
#' @param discrete logical, if true, results should be analyzed as integers
#' @param bin number of bins in the test
#' @export
#' @rdname analyzeResidualError
analyzeResidualErrorByPPplot = function(trait,pred,error,epsilon=0,laplace=FALSE,discrete=FALSE,n=1000){
  if(is.null(epsilon)) epsilon=0
  if(discrete){
    if(is.list(error)){
      discretized.res = discretizeResult(res = data.frame(node=1,label="1",x=pred),
                                         error = error,epsilon = epsilon,laplace = laplace)
    }else{
      discretized.res = discretizeResult(res = data.frame(node=1,label="1",x=pred,var=error),
                                         error = NULL,epsilon = epsilon,laplace = laplace)
    }
    res = calculateBinomRSS(p = discretized.res$probs,obs = discretized.res$x==trait)
  }else{
    pZ = pseudo.Z.score(trait,pred,error,epsilon,laplace)
    res = ks.test(pZ$q,pnorm)$statistic
  }
  return(res)
}

#' @export
#' @rdname analyzeResidualError
analyzeResidualErrorByHeteroscedasticity = function(trait,pred,error,epsilon=0,laplace=FALSE,bin=100,n=1000){
  if(is.null(epsilon)) epsilon=0
  pZ = pseudo.Z.score(trait,pred,error,epsilon,laplace)
  if(is.list(error)){error.as.var = sapply(error,function(x){
    unname(compoundNormalMoments(probs = x$probs,means = x$x,vars = x$var)[2])
  })}else{
    error.as.var=error
  }
  rss = calculateHeteroscedasticity(x = pZ$q,y = error.as.var,bin = bin,n = n)
  return(rss)
}

#' @export
#' @rdname analyzeResidualError
analyzeResidualErrorByCI = function(trait,pred,error,epsilon=0,laplace=FALSE,bin=20,alpha=0.05){
  if(is.null(epsilon)) epsilon=0
  pZ = pseudo.Z.score(trait,pred,error,epsilon,laplace)
  if(is.list(error)){error.as.var = sapply(error,function(x){
    unname(compoundNormalMoments(probs = x$probs,means = x$x,vars = x$var)[2])
  })}else{
    error.as.var=error
  }
  rss = calculateBinomRSS(p = rep(alpha,length(error.as.var)),obs = (sapply(pZ$p,function(x){min(x,1-x)})<(alpha/2)),group = log(error.as.var),bin = bin)
  return(rss)
}
