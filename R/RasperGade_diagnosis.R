#
logitSigmoid=function(x,p){
  sapply(x,function(xx){
    if(is.na(xx)|is.nan(xx)|is.nan(p)|is.na(p)) return(NA)
    logit.x = logit(x = xx,a = 1,b = -logit(p))
    trans.x = sigmoid(logit.x)
  })
}
#
reverse.logitSigmoid=function(x,p){
  sapply(x,function(xx){
    if(is.na(xx)|is.nan(xx)|is.nan(p)|is.na(p)) return(NA)
    logit.x = logit(x = xx,a = 1,b = logit(p))
    trans.x = sigmoid(logit.x)
  })
}
#
calculateBinomRSS = function(p,obs,group=NULL,bin=20,...){
  if(is.null(group)){
    bin.rss = sapply(seq(0.005,0.995,0.01),function(z){
      idx = (p<=(z+0.005))&(p>(z-0.005))
      total.count = sum(idx)
      exp.z = z
      if(total.count>0) exp.z = mean(p[idx])
      obs.count = 0
      if(total.count>0) obs.count = sum(obs[idx])
      rss = (2*logitSigmoid(x = obs.count/total.count,p = exp.z)-1)^2
      return(rss)
    })
    this.rss = 100*mean(bin.rss,na.rm = TRUE)
  }else{
    group.idx = sort(group,decreasing = FALSE,index.return=TRUE)$ix
    bin.idx = round(seq(from=1,to = length(group)+1,length.out = bin+1))
    bin.rss = sapply(1:bin,function(z){
                                idx = group.idx[bin.idx[z]:(bin.idx[z+1]-1)]
                                total.count = length(idx)
                                exp.z = mean(p)
                                if(total.count>0) exp.z = mean(p[idx])
                                obs.count = 0
                                if(total.count>0) obs.count = sum(obs[idx])
                                rss = (2*logitSigmoid(x = obs.count/total.count,p = exp.z)-1)^2
                                return(rss)
                              })
    this.rss = bin*mean(bin.rss,na.rm = TRUE)
  }
  return(this.rss)
}
#
calculateHeteroscedasticity = function(x,y,bin=100,n=1000){
  idx = sort(y,index.return=TRUE)$ix
  bins=round(seq(0,length(x),length.out = bin+1))
  vars = lapply(1:bin,function(i){
    x[idx[(bins[i]+1):bins[i+1]]]
  })
  h = fligner.test(x = vars)
  return(c(h=h$statistic,p.value=h$p.value))
}
#
calculateExclusion = function(pred,error,epsilon=0,tolerance,alpha=0.05,laplace=FALSE,discrete=FALSE){
  if(discrete){
    pred = round(pred)
    tolerance = ceiling(tolerance)+0.5
  }
  if(laplace){
    if(is.list(error)){
      prob.within.range = sapply(1:length(pred),function(i){
        pMixNormalLaplace(q = pred[i]+tolerance,mean = error[[i]]$x,sd = sqrt(error[[i]]$var),
                          laplace.var = epsilon/2,probs = error[[i]]$probs)
      }) - 
        sapply(1:length(pred),function(i){
          pMixNormalLaplace(q = pred[i]-tolerance,mean = error[[i]]$x,sd = sqrt(error[[i]]$var),
                            laplace.var = epsilon/2,probs = error[[i]]$probs)
        })
    }else{
      prob.within.range = sapply(error,function(x){pconv.norm.laplace(q = tolerance,sigma = x,epsilon = epsilon/2)})-
        sapply(error,function(x){pconv.norm.laplace(q = -tolerance,sigma = x,epsilon = epsilon/2)})
    }
  }else{
    if(is.list(error)){
      prob.within.range = sapply(1:length(pred),function(i){
        pMixNormal(q = pred[i]+tolerance,mean = error[[i]]$x,sd = sqrt(error[[i]]$var+epsilon/2),probs = error[[i]]$probs)
      }) - 
        sapply(1:length(pred),function(i){
          pMixNormal(q = pred[i]-tolerance,mean = error[[i]]$x,sd = sqrt(error[[i]]$var+epsilon/2),probs = error[[i]]$probs)
        })
    }else{
      prob.within.range = pnorm(q = tolerance,mean = 0,sd = sqrt(error+epsilon/2))-
        pnorm(q = -tolerance,mean = 0,sd = sqrt(error+epsilon/2))
    }
  }
  return(list(p=prob.within.range, decision=prob.within.range < (1-alpha)))
}
#
checkExclusion = function(obs,pred,tolerance,discrete=FALSE){
  if(discrete){
    pred = round(pred)
    tolerance = ceiling(tolerance)
  }
  return(abs(obs-pred) > tolerance)
}
#
calculateExclusionErrorRate = function(obs,pred,error,epsilon=0,tolerance,alpha=0.05,laplace=FALSE,discrete=FALSE){
  observed.exclusion = checkExclusion(obs,pred,tolerance,discrete)
  predicted.exclusion = calculateExclusion(pred,error,epsilon,tolerance,alpha,laplace,discrete)
  predicted.p = predicted.exclusion$p
  predicted.exclusion = predicted.exclusion$decision
  errors = c(true.exclusion=sum(observed.exclusion&predicted.exclusion),
             false.exclusion=sum((!observed.exclusion)&predicted.exclusion),
             true.inclusion=sum((!observed.exclusion)&(!predicted.exclusion)),
             false.inclusion=sum(observed.exclusion&(!predicted.exclusion)))
  expected.errors = c(true.exclusion=sum(1-predicted.p[predicted.exclusion]),
                      false.exclusion=sum(predicted.p[predicted.exclusion]),
                      true.inclusion=sum(predicted.p[!predicted.exclusion]),
                      false.inclusion=sum(1-predicted.p[!predicted.exclusion]))
  rates=c(FP=sum((!observed.exclusion)&predicted.exclusion)/sum(!observed.exclusion),
          FN=sum(observed.exclusion&(!predicted.exclusion))/sum(observed.exclusion),
          FDR=sum((!observed.exclusion)&predicted.exclusion)/sum(predicted.exclusion),
          FRR=sum(observed.exclusion&(!predicted.exclusion))/sum(!predicted.exclusion))
  return(list(obs=observed.exclusion,pred=predicted.exclusion,error=errors,expect=expected.errors,rate=rates))
}
#
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
#
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
#
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
#
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
