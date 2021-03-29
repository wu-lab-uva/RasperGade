#' @title  Find approximation for a mixture of normal distributions using limited number of normal distributions
#'
#' @description Find a solution of using 1 or 2 normal distributions to approximate a mixture of normal distributions
#'
#' @param means the vector of means of mixed normal distributions
#' @param vars the vector of standard deviations of mixed normal distributions
#' @param probs the vector of probability weight of normal distributions
#' @return A dataframe with 3 columns, with columns representing the means, variances and probability weights of the normal distributions to approximate the full mixture
#' @export
#' @rdname findBestNormalApproximation
findBestNormalApproximation1 = function(probs,means,vars){
  true.moments = compoundNormalMoments(probs = probs,means = means,vars = vars)
  return(data.frame(mean=unname(true.moments[1]),var=unname(true.moments[2]),probs=1))
}

#' @export
#' @rdname findBestNormalApproximation
findBestNormalApproximation2 = function(probs,means,vars,numCores=1){
  if(length(probs)<2) return(findBestNormalApproximation1(probs=probs,means=means,vars=vars))
  true.moments = compoundNormalMoments(probs = probs,means = means,vars = vars)
  sample.range = sort(unique(means),decreasing = FALSE)
  sample.range = sort(c(sample.range,(sample.range[-1]+sample.range[-length(sample.range)])/2),decreasing = FALSE)
  sample.likelihood = c(0,sapply(sample.range,function(x){
    sum(dnorm(x = x,mean = means,sd = sqrt(vars))*probs)
  }),0)
  sample.peak = which(sapply(2:(length(sample.range)+1),function(i){
    (sample.likelihood[i]>=sample.likelihood[i+1])&(sample.likelihood[i]>=sample.likelihood[i-1])
  }))
  if(length(sample.peak)<2){
    means.order = sort(abs(means-sample.range[sample.peak]),index.return=TRUE)$ix
    if(length(probs)>1e3){
      sample.break = round(seq(from=1,to = length(means.order)-1,length.out = 1000))
    }else{
      sample.break = 1:(length(means.order)-1)
    }
    moments.diff = mclapply(sample.break,function(i){
      this.idx = means.order[1:i]
      moments.1 = compoundNormalMoments(probs = probs[this.idx]/sum(probs[this.idx]),means = means[this.idx],vars = vars[this.idx])
      moments.2 = compoundNormalMoments(probs = probs[-(this.idx)]/sum(probs[-(this.idx)]),means = means[-(this.idx)],vars = vars[-(this.idx)])
      new.df = data.frame(mean=unname(c(moments.1[1],moments.2[1])),var=unname(c(moments.1[2],moments.2[2])),
                          probs=c(sum(probs[this.idx]),sum(probs[-(this.idx)])))
      new.likelihood = c(0,sapply(sample.range,function(x){
        sum(dnorm(x = x,mean = new.df$mean,sd = sqrt(new.df$var))*new.df$probs)
      }),0)
      rss = sum((sample.likelihood-new.likelihood)^2)
      return(list(diff=rss,df=new.df))
    },mc.cores = numCores)
  }else{
    peak.combn = combn(x = sample.peak,m = 2,simplify = FALSE)
    moments.diff = mclapply(peak.combn,function(i){
      new.df = do.call(rbind,lapply(i,function(j){
        this.group = which(((means-mean(sample.range[i]))*(sample.range[j]-mean(sample.range[i])))>0)
        half.group = which(((means-mean(sample.range[i]))*(sample.range[j]-mean(sample.range[i])))==0)
        this.moment = compoundNormalMoments(probs = c(probs[this.group],probs[half.group]/2),
                                            means = means[c(this.group,half.group)],vars = vars[c(this.group,half.group)])
        data.frame(mean=unname(this.moment[1]),var=unname(this.moment[2]),
                   probs=sum(c(probs[this.group],probs[half.group]/2)))
      }))
      new.likelihood = c(0,sapply(sample.range,function(x){
        sum(dnorm(x = x,mean = new.df$mean,sd = sqrt(new.df$var))*new.df$probs)
      }),0)
      rss = sum((sample.likelihood-new.likelihood)^2)
      return(list(diff=rss,df=new.df))
    },mc.cores = numCores)
  }
  best.diff = which.min(sapply(moments.diff,function(x){x$diff}))
  return(moments.diff[[best.diff]]$df)
}
