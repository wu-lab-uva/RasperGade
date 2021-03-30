#
library(RasperGade)
#
all.data = readRDS("Data/Landis.data.RDS")[readRDS("Landis.PE.index.RDS")]
raw.model = readRDS("Result/Empirical/BodySize.models.RDS")
all.stat = normalize.AIC(AIC = raw.model$stat,scale = raw.model$scale,
                         sample.size = sum(sapply(all.data,function(x){Nnode(x$phy)})))
all.model = lapply(raw.model$params[c(1,3)],function(x){x/c(1,rep(raw.model$scale,3))})
bm.rate = var(do.call(c,lapply(all.data,function(x){pic(x = x$dat[x$phy$tip.label],phy = x$phy)})))
bm.pic = do.call(rbind,lapply(all.data,function(dat){
  pic(x = dat$dat[dat$phy$tip.label],phy = dat$phy,scaled = FALSE,var.contrasts = TRUE)
}))
bm.stat = -2*(sum(dnorm(x = bm.pic[,1],sd = sqrt(bm.rate*bm.pic[,2]),log = TRUE)))+2
#
if(file.exists("Result/Empirical/Landis.CV.RDS")){
  all.asr = readRDS("Result/Empirical/Landis.CV.RDS")
}else{
  all.asr = list()
  for(clade in names(all.data)){
    print(clade)
    tree = all.data[[clade]]$phy
    trait = all.data[[clade]]$dat[all.data[[clade]]$phy$tip.label]
    asr.res = list(list(FMR=list(phy=tree,params=list(sigma=bm.rate,epsilon=0)),trait=trait,
                        CV=do.call(rbind,lapply(1:Ntip(tree),function(i){
                          predTraitByPIC(phy = tree,trait = trait[-i],rate=bm.rate)}))
    ))
    for(mdl in all.model){
      t0 = Sys.time()
      FMR = fullMarginalReconstructionWithPE(phy = tree,x = trait,
                                             params = mdl,
                                             laplace = FALSE,approximate = 1,numCores = 1,asymptotic = 5)
      CV = crossValidationWithPE(FMR = FMR,add.epsilon = FALSE,numCores = 1,laplace = FALSE,numApprox = 1)
      t1 = Sys.time()
      asr.res[[length(asr.res)+1]] = list(FMR=FMR,CV=CV,trait = trait,t = c(t0,t1))
    }
    all.asr[[clade]] = asr.res
  }
  saveRDS(all.asr,"Result/Empirical/Landis.CV.RDS")
}
# read in expected statistics
stat.expectation = readRDS("Result/stat.expectation.RDS")
#
cv.all.q = do.call(rbind,lapply(all.asr,function(x){
  sapply(x,function(y){
    if(is.data.frame(y$CV)){
      pseudo.Z.score(trait = x[[2]]$trait[x[[2]]$FMR$phy$tip.label],pred = y$CV$x,
                                   error = y$CV$var,epsilon = 0,laplace = FALSE)$q
    }else{
      pseudo.Z.score(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                                   error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE)$q
    }
  })
}))
cv.d = apply(cv.all.q,2,function(x){
  ks.test(x,pnorm)$statistic
})
cv.d/mean(stat.expectation$cv.D)

# Figure 1 part
BM.Felsenstein.pZ = lapply(list(1:length(all.asr)),function(idx){
  do.call(rbind,lapply(all.asr[idx],function(x){
    Z=(x[[1]]$trait-x[[1]]$CV$x)/sqrt(x[[1]]$CV$var)
    data.frame(p=pnorm(Z),Z=Z,var = x[[1]]$CV$var,se = sqrt(x[[1]]$CV$var))
  }))
})
saveRDS(BM.Felsenstein.pZ,"Landis.BM.plot.data.RDS")
calculateHeteroscedasticity(x = BM.Felsenstein.pZ[[1]]$Z,y = BM.Felsenstein.pZ[[1]]$se,bin = 20)[1]/mean(stat.expectation$cv.h)
# Figure 2 part
PEBME.pZ = lapply(list(1:length(all.asr)),function(idx){
  do.call(rbind,lapply(all.asr[idx],function(x){
    pZ=sapply(1:dim(x[[3]]$CV$cv)[1],function(i){
      pMixNormal(q = x[[3]]$trait[i],mean = x[[3]]$CV$error[[i]]$x,
                 sd = sqrt(x[[3]]$CV$error[[i]]$var+x[[3]]$FMR$params$epsilon/2),
                 probs = x[[3]]$CV$error[[i]]$probs)
    })
    data.frame(p=pZ,Z=qnorm(pZ),var = x[[1]]$CV$var,se = sqrt(x[[1]]$CV$var))
  }))
})
calculateHeteroscedasticity(x = PEBME.pZ[[1]]$Z,y = PEBME.pZ[[1]]$se,bin = 20)[1]/mean(stat.expectation$cv.h)

saveRDS(PEBME.pZ,"Result/Empirical/Landis.PEBME.plot.data.RDS")
#
save.image("Landis.summary.RData")


