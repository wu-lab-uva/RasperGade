#
setwd("F:/Dropbox/Scripts/16S_Copy_Number/Version3/")
source("fitPE_functions.R")
setwd("F:/Dropbox/Scripts/16S_Copy_Number/Version4/")
source("RasperGade_reconstruction.R")
source("RasperGade_diagnosis.R")
setwd("F:/Dropbox/Scripts/16S_Copy_Number/20210302/")
#
job.names = c("BacNormal/bac.16SGCN.jitter.normal",
              "BacNormal/bac.gGC.rec.normal",
              "BacNormal/bac.GS.rec.normal",
              "BacNormal/bac.N.rec.normal",
              "BacNormal/bac.rGC.rec.normal")
names(job.names) = sapply(job.names,function(x){strsplit(x,".",fixed = TRUE)[[1]][2]})
trait.label = c("16S rRNA GCN","Genomic GC%","Genome size","N-ARSC","rRNA GC%","Body size")
#
all.data = lapply(job.names,function(f){
  readRDS(paste0(f,".data.RDS"))
})
#
all.model = lapply(job.names,function(f){
  res = readRDS(paste0(f,".models.RDS"))
  mdls = lapply(res$params[c(1,3)],function(mdl){
    params=summarize.model.parameter(mdl,scale=res$scale)$params
    if(!is.null(res$jitter)) params["epsilon"] = params["epsilon"]-2*res$jitter
    return(params)
  })
  names(mdls) = c("BMe","PEBMe")
  return(mdls)
})
#
tiv.contribution = sapply(1:5,function(i){
  all.model[[i]]$PEBMe[4]/var(all.data[[i]]$dat)
})
names(tiv.contribution) = names(all.data)
#
sample.size = sapply(1:length(all.data),function(i){
  dat = all.data[[i]]
  tpic = findTerminalPIC(phy = dat$phy,x = dat$dat[dat$phy$tip.label],remove.zero = FALSE)
  if(i<=1) return(Nnode(dat$phy))
  return(Nnode(dat$phy)-sum(tpic$pic==0))
})
names(sample.size) = names(all.data)
#
all.bm.fit = lapply(1:length(all.data),function(i){
  dat = all.data[[i]]
  bm.pic = pic(x = dat$dat[dat$phy$tip.label],phy = dat$phy,scaled = FALSE,var.contrasts = TRUE)
  bm.rate = var(pic(x = dat$dat[dat$phy$tip.label],phy = dat$phy))
  tpic = findTerminalPIC(phy = dat$phy,x = dat$dat[dat$phy$tip.label],remove.zero = FALSE)
  if(i>1&(sum(tpic$pic==0)>0)){
    bm.aic = -2*(sum(dnorm(x = bm.pic[-(tpic$node[tpic$pic==0]-Ntip(dat$phy)),1],
                           sd = sqrt(bm.rate*bm.pic[-(tpic$node[tpic$pic==0]-Ntip(dat$phy)),2]),log = TRUE)))+2
  }else{
    bm.aic = -2*(sum(dnorm(x = bm.pic[,1],
                           sd = sqrt(bm.rate*bm.pic[,2]),log = TRUE)))+2
  }
  return(data.frame(bm.aic=bm.aic,sigma=bm.rate))
})
#
all.scale = sapply(job.names,function(f){
  res = readRDS(paste0(f,".models.RDS"))
  return(res$scale)
})
all.stat = lapply(job.names,function(f){
  res = readRDS(paste0(f,".models.RDS"))
  aics = res$stat[c(1,3)]
  names(aics) = c("BMe","PEBMe")
  return(aics)
})
all.tabs1.df = do.call(rbind,lapply(1:5,function(i){
  norm.aic = normalize.AIC(AIC = all.stat[[i]],scale = all.scale[i],sample.size = sample.size[i])
  data.frame(trait=names(job.names)[i],
             bm.aic = unname(all.bm.fit[[i]]$bm.aic),bme.aic = unname(norm.aic[1]),pebme.aic = unname(norm.aic[2]),
             bm.sigma = unname(all.bm.fit[[i]]$sigma),
             bme.sigma=unname(all.model[[i]]$BMe[3]),bme.epsilon=unname(all.model[[i]]$BMe[4]),
             pebme.lambda = unname(all.model[[i]]$PEBMe[1]),pebme.size = unname(all.model[[i]]$PEBMe[2]),
             pebme.sigma = unname(all.model[[i]]$PEBMe[3]),pebme.epsilon = unname(all.model[[i]]$PEBMe[4])
             )
}))
write.table(x = all.tabs1.df,file = "bacteria.parameters.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
isDiscrete = c(TRUE,FALSE,FALSE,FALSE,FALSE)
# reconstruct ancestral states and conduct cross-validation 
if(file.exists("Bacteria.CV.RDS")){
  all.asr = readRDS("Bacteria.CV.RDS")
}else{
  all.asr = list()
  for(traitname in names(all.data)){
    print(traitname)
    tree = all.data[[traitname]]$phy
    trait = all.data[[traitname]]$dat[all.data[[traitname]]$phy$tip.label]
    t0 = Sys.time()
    bm.rate = var(pic(x = trait,phy = tree))
    asr.res = list(list(FMR=list(phy=tree,params=list(sigma=bm.rate,epsilon=0)),trait=trait,
                        CV=do.call(rbind,lapply(1:Ntip(tree),function(i){
                          predTraitByPIC(phy = tree,trait = trait[-i],rate=bm.rate)}))
    ))
    t1 = Sys.time()
    asr.res[[1]]$t = c(t0,t1)
    for(mdl in all.model[[traitname]]){
      t0 = Sys.time()
      names(mdl) = c("lambda","size","sigma","epsilon")
      FMR = fullMarginalReconstructionWithPE(phy = tree,x = trait,
                                             params = mdl,
                                             laplace = FALSE,approximate = 1,numCores = 1,asymptotic = 5)
      CV = crossValidationWithPE(FMR = FMR,add.epsilon = FALSE,numCores = 1,laplace = FALSE,numApprox = 1)
      t1 = Sys.time()
      asr.res[[length(asr.res)+1]] = list(FMR=FMR,CV=CV,trait = trait,t = c(t0,t1))
    }
    all.asr[[traitname]] = asr.res
  }
  saveRDS(all.asr,"Bacteria.CV.RDS")
}
# read in expected statistics
stat.expectation = readRDS("stat.expectation.RDS")
#
library(ggplot2)
library(ggpubr)
#
ref.size = 0.6
line.size = 0.6
# Figure 2
BM.Felsenstein.pZ = c(lapply(all.asr[-1],function(x){
  Z = (x[[1]]$trait[x[[1]]$FMR$phy$tip.label]-x[[1]]$CV$x)/sqrt(x[[1]]$CV$var)
  pZ.df = data.frame(p=pnorm(Z),Z=Z,var = x[[1]]$CV$var,se = sqrt(x[[1]]$CV$var))
}),
list(do.call(rbind,readRDS("Landis.BM.plot.data.RDS"))))
#
BM.Felsenstein.Z.plot = lapply(1:length(BM.Felsenstein.pZ),function(i){
  df = BM.Felsenstein.pZ[[i]]
  ggplot()+geom_point(mapping = aes(x=se,y = Z),data = df,alpha=0.5,shape=".")+
    geom_hline(yintercept = 2,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    geom_hline(yintercept = -2,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(ylim = c(-10,10),
                    xlim=c(10^(floor(2*log10(min(df$se)))/2),10^(ceiling(2*log10(max(df$se)))/2)))+
    scale_x_continuous(trans = "log10",breaks=c(1e-4,0.001,0.01,0.1,1,10),
                       labels = c(expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
    scale_y_continuous(breaks = seq(-10,10,2))+
    ggtitle(trait.label[i+1])+xlab("Predicted SE")+ylab("Z-score")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
BM.Felsenstein.pp.plot = lapply(1:length(BM.Felsenstein.pZ),function(i){
  df = BM.Felsenstein.pZ[[i]]
  this.df = data.frame(empirical = (1:dim(df)[1])/(dim(df)[1]+1),theoretical=sort(df$p,decreasing = FALSE))
  ggplot()+geom_line(mapping = aes(x=theoretical,y=empirical),data=this.df,alpha=0.8,size=line.size)+
    geom_abline(slope = 1,intercept = 0,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(xlim = c(0,1),ylim=c(0,1),expand = FALSE)+
    scale_x_continuous(breaks=seq(0,1,0.25),
                       labels = seq(0,1,0.25),expand = c(0,0.05))+
    scale_y_continuous(breaks=seq(0,1,0.25),
                       labels = seq(0,1,0.25),expand = c(0,0.05))+
    xlab("Theoretical probability")+ylab("Empirical probability")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
BM.Felsenstein.qq.plot = lapply(1:length(BM.Felsenstein.pZ),function(i){
  df = BM.Felsenstein.pZ[[i]]
  this.df = data.frame(empirical = sort(df$Z,decreasing = FALSE),theoretical= qnorm((1:dim(df)[1])/(dim(df)[1]+1)))
  ggplot()+geom_line(mapping = aes(x=theoretical,y=empirical),data=this.df,alpha=0.8,size=line.size)+
    geom_abline(slope = 1,intercept = 0,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    #coord_cartesian(ylim=c(-10,10),xlim=c(-3.25,3.25))+
    xlab("Theoretical quantile")+ylab("Empirical quantile")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
BM.Felsenstein.alpha.plot = lapply(1:length(BM.Felsenstein.pZ),function(i){
  df = BM.Felsenstein.pZ[[i]]
  bins.se = exp(seq(log(min(df$se*2)),log(max(df$se/2)),length.out = 20))
  se.idx = lapply(1:20,function(j){
    which((df$se<=(bins.se[j]*2))&(df$se>=(bins.se[j]/2)))
  })
  valid.idx = sapply(se.idx,length)>=100
  this.df = data.frame(se=bins.se[valid.idx],sample.size = sapply(se.idx[valid.idx],length),
                       cp = sapply(se.idx[valid.idx],function(x){
                         logitSigmoid(x = 1-(sum(sapply(df$p[x],function(y){min(y,1-y)<(0.05/2)})))/(length(x)),
                                      p = 0.95)
                       }),
                       upCI = sapply(se.idx[valid.idx],function(x){
                         logitSigmoid(x = 1 - binom.test(x = sum(sapply(df$p[x],function(y){min(y,1-y)<(0.05/2)})),
                                                         n = length(x),p = 0.05,conf.level = 0.95)$conf.int[1],
                                          p = 0.95)
                       }),
                       lowCI = sapply(se.idx[valid.idx],function(x){
                         logitSigmoid(x = 1 - binom.test(x = sum(sapply(df$p[x],function(y){min(y,1-y)<(0.05/2)})),
                                                         n = length(x),p = 0.05,conf.level = 0.95)$conf.int[2],
                                          p = 0.95)
                       }))
  print(range(reverse.logitSigmoid(this.df$cp,p = 0.95)))
  ggplot()+geom_line(mapping = aes(x=se,y=cp),data=this.df,alpha=0.8,size=line.size)+
    geom_line(mapping = aes(x=se,y=upCI),data=this.df,alpha=0.6,size=line.size,linetype="dashed")+
    geom_line(mapping = aes(x=se,y=lowCI),data=this.df,alpha=0.6,size=line.size,linetype="dashed")+
    geom_hline(yintercept = 0.5,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(ylim=c(0,1),xlim=c(10^floor(log10(min(this.df$se))),10^ceiling(log10(max(this.df$se)))))+
    scale_x_continuous(trans = "log10",breaks=c(0.001,0.01,0.1,1,10),
                       labels = c(expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
    scale_y_continuous(breaks = logitSigmoid(x = c(0.5,0.75,0.9,0.95,0.97,0.99,1),p = 0.95),
                       labels = c(0.5,0.75,0.9,0.95,0.97,0.99,1))+
    xlab("Predicted SE")+ylab("Coverage probability")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
BM.Felsenstein.combined.plot = ggarrange(plotlist = do.call(c,lapply(c(4,1,2,3,5),function(i){
  list(BM.Felsenstein.Z.plot[[i]],BM.Felsenstein.qq.plot[[i]],BM.Felsenstein.alpha.plot[[i]])
})),ncol = 3,nrow = 5,labels = "AUTO",align = "hv")
#
ggsave(filename = "Fig2.pdf",device = "pdf",plot = BM.Felsenstein.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 5.5,units = "in")
ggsave(filename = "Fig2.png",device = "png",plot = BM.Felsenstein.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 5.5,units = "in",dpi = "print",type="cairo")
# Figure 4
PEBME.pZ = c(lapply(all.asr[-1],function(x){
  Z = pseudo.Z.score(trait=x[[3]]$trait[x[[3]]$FMR$phy$tip.label],pred = x[[3]]$CV$cv$x,
                     error = x[[3]]$CV$error,epsilon = x[[3]]$FMR$params$epsilon,laplace = FALSE)
  pZ.df = data.frame(p=Z$p,Z=Z$q,var = x[[3]]$CV$cv$var,se = sqrt(x[[3]]$CV$cv$var))
}),
list(do.call(rbind,readRDS("Landis.PEBME.plot.data.RDS"))))
#
PEBME.Z.plot = lapply(1:length(PEBME.pZ),function(i){
  df = PEBME.pZ[[i]]
  ggplot()+geom_point(mapping = aes(x=se,y = Z),data = df,alpha=0.5,shape=".")+
    geom_hline(yintercept = 2,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    geom_hline(yintercept = -2,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(ylim = c(-6,6),
                    xlim=c(10^(floor(2*log10(min(df$se)))/2),10^(ceiling(2*log10(max(df$se)))/2)))+
    scale_x_continuous(trans = "log10",breaks=c(1e-4,0.001,0.01,0.1,1,10),
                       labels = c(expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
    scale_y_continuous(breaks = seq(-10,10,2))+
    ggtitle(trait.label[i+1])+xlab("Predicted SE")+ylab("Pseudo-Z-score")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
PEBME.pp.plot = lapply(1:length(PEBME.pZ),function(i){
  df = PEBME.pZ[[i]]
  this.df = data.frame(empirical = (1:dim(df)[1])/(dim(df)[1]+1),theoretical=sort(df$p,decreasing = FALSE))
  ggplot()+geom_line(mapping = aes(x=theoretical,y=empirical),data=this.df,alpha=0.8,size=line.size)+
    geom_abline(slope = 1,intercept = 0,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(xlim = c(0,1),ylim=c(0,1),expand = FALSE)+
    scale_x_continuous(breaks=seq(0,1,0.25),
                       labels = seq(0,1,0.25),expand = c(0,0.02))+
    scale_y_continuous(breaks=seq(0,1,0.25),
                       labels = seq(0,1,0.25),expand = c(0,0.02))+
    xlab("Theoretical probability")+ylab("Empirical probability")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
PEBME.qq.plot = lapply(1:length(PEBME.pZ),function(i){
  df = PEBME.pZ[[i]]
  this.df = data.frame(empirical = sort(df$Z,decreasing = FALSE),theoretical= qnorm((1:dim(df)[1])/(dim(df)[1]+1)))
  ggplot()+geom_line(mapping = aes(x=theoretical,y=empirical),data=this.df,alpha=0.8,size=line.size)+
    geom_abline(slope = 1,intercept = 0,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    #coord_cartesian(ylim=c(-10,10),xlim=c(-3.25,3.25))+
    xlab("Theoretical quantile")+ylab("Empirical quantile")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
PEBME.alpha.plot = lapply(1:length(PEBME.pZ),function(i){
  df = PEBME.pZ[[i]]
  bins.se = exp(seq(log(min(df$se*2)),log(max(df$se/2)),length.out = 20))
  se.idx = lapply(1:20,function(j){
    which((df$se<=(bins.se[j]*2))&(df$se>=(bins.se[j]/2)))
  })
  valid.idx = sapply(se.idx,length)>=100
  this.df = data.frame(se=bins.se[valid.idx],sample.size = sapply(se.idx[valid.idx],length),
                       cp = sapply(se.idx[valid.idx],function(x){
                         logitSigmoid(x = 1-(sum(sapply(df$p[x],function(y){min(y,1-y)<(0.05/2)})))/(length(x)),
                                      p = 0.95)
                       }),
                       upCI = sapply(se.idx[valid.idx],function(x){
                         logitSigmoid(x = 1 - binom.test(x = sum(sapply(df$p[x],function(y){min(y,1-y)<(0.05/2)})),
                                                         n = length(x),p = 0.05,conf.level = 0.95)$conf.int[1],
                                      p = 0.95)
                       }),
                       lowCI = sapply(se.idx[valid.idx],function(x){
                         logitSigmoid(x = 1 - binom.test(x = sum(sapply(df$p[x],function(y){min(y,1-y)<(0.05/2)})),
                                                         n = length(x),p = 0.05,conf.level = 0.95)$conf.int[2],
                                      p = 0.95)
                       }))
  print(range(reverse.logitSigmoid(this.df$cp,p = 0.95)))
  ggplot()+geom_line(mapping = aes(x=se,y=cp),data=this.df,alpha=0.8,size=line.size)+
    geom_line(mapping = aes(x=se,y=upCI),data=this.df,alpha=0.6,size=line.size,linetype="dashed")+
    geom_line(mapping = aes(x=se,y=lowCI),data=this.df,alpha=0.6,size=line.size,linetype="dashed")+
    geom_hline(yintercept = 0.5,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(ylim=c(0,1),xlim=c(10^floor(log10(min(this.df$se))),10^ceiling(log10(max(this.df$se)))))+
    scale_x_continuous(trans = "log10",breaks=c(0.001,0.01,0.1,1,10),
                       labels = c(expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
    scale_y_continuous(breaks = logitSigmoid(x = c(0.5,0.75,0.9,0.95,0.97,0.99,1),p = 0.95),
                       labels = c(0.5,0.75,0.9,0.95,0.97,0.99,1))+
    xlab("Predicted SE")+ylab("Coverage probability")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
PEBME.combined.plot = ggarrange(plotlist = do.call(c,lapply(c(4,1,2,3,5),function(i){
  list(PEBME.Z.plot[[i]],PEBME.qq.plot[[i]],PEBME.alpha.plot[[i]])
})),ncol = 3,nrow = 5,labels = "AUTO",align = "hv")
#
ggsave(filename = "Fig4.pdf",device = "pdf",plot = PEBME.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 5.5,units = "in")
ggsave(filename = "Fig4.png",device = "png",plot = PEBME.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 5.5,units = "in",dpi = "print",type="cairo")
#
tab2.left = do.call(rbind,lapply(BM.Felsenstein.pZ[c(4,1,2,3,5)],function(df){
  this.D = ks.test(df$Z,pnorm)$statistic/mean(stat.expectation$cv.D)
  this.h = unname(calculateHeteroscedasticity(x = df$Z,y = df$se,bin = 20)[1])/mean(stat.expectation$cv.h)
  this.rss = calculateBinomRSS(p = rep(0.95,dim(df)[1]),obs = ((df$p<=0.975)&(df$p>=0.025)), group = df$se,bin = 20)/
    mean(stat.expectation$cv.rss)
  data.frame(D=this.D,h=this.h,rss=this.rss)
}))
#
tab2.right = do.call(rbind,lapply(PEBME.pZ[c(4,1,2,3,5)],function(df){
  this.D = ks.test(df$Z,pnorm)$statistic/mean(stat.expectation$cv.D)
  this.h = unname(calculateHeteroscedasticity(x = df$Z,y = df$se,bin = 20)[1])/mean(stat.expectation$cv.h)
  this.rss = calculateBinomRSS(p = rep(0.95,dim(df)[1]),obs = ((df$p<=0.975)&(df$p>=0.025)), group = df$se,bin = 20)/
    mean(stat.expectation$cv.rss)
  data.frame(D=this.D,h=this.h,rss=this.rss)
}))
#
gcn.RSS = sapply(all.asr$`16SGCN`,function(res){
  if(is.data.frame(res$CV)){
    analyzeResidualErrorByPPplot(trait = res$trait,pred = res$CV$x,error = res$CV$var,
                                 epsilon = res$FMR$params$epsilon,laplace = FALSE,discrete = TRUE)
  }else{
    analyzeResidualErrorByPPplot(trait = res$trait,pred = res$CV$cv$x,error = res$CV$error,
                                 epsilon = res$FMR$params$epsilon,laplace = FALSE,discrete = TRUE)
  }
})
gcn.error.rate = lapply(all.asr$`16SGCN`,function(res){
  if(is.data.frame(res$CV)){
    calculateExclusionErrorRate(obs = res$trait,pred = res$CV$x,error = res$CV$var,
                                epsilon = res$FMR$params$epsilon,tolerance = 0,alpha = 0.05,laplace = FALSE,discrete = TRUE)
  }else{
    calculateExclusionErrorRate(obs = res$trait,pred = res$CV$cv$x,error = res$CV$error,
                                epsilon = res$FMR$params$epsilon,tolerance = 0,alpha = 0.05,laplace = FALSE,discrete = TRUE)
  }
})
sapply(gcn.error.rate,function(x){
  sum(x$error[c(2,4)])/sum(x$error)
})
#
save.image("Bacteria.summary.RData")
