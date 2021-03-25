#
source("fitPE_functions.R")
source("RasperGade_reconstruction.R")
source("RasperGade_diagnosis.R")
library(ggplot2)
library(ggpubr)
# find representative from simulated data
epsilons=c(0,1,10)
contributions = rev(c(0,0.1,0.3,0.5,0.7,0.9,1))
replicates = 1:20
batch = 1:5
pars.epsilon = expand.grid(epsilon=epsilons,lambda=5,contribution=0.9,rep=replicates,batch=batch,KEEP.OUT.ATTRS = FALSE)
pars.epsilon$id = rep(1:(dim(pars.epsilon)[1]/length(unique(pars.epsilon$batch))),length(unique(pars.epsilon$batch)))
pars.contribution = expand.grid(epsilon=10,lambda=5,contribution=contributions,rep=replicates,batch=batch,KEEP.OUT.ATTRS = FALSE)
pars.contribution$id = rep(1:(dim(pars.contribution)[1]/length(unique(pars.contribution$batch))),length(unique(pars.contribution$batch)))
pars.bm = expand.grid(rate=1,batch=1,id=1:100)
# set up file path
rep.file.idx = list(1:100,
                    which(pars.contribution[,3]==0),
                    which(pars.epsilon[,1]==0),
                    which(pars.contribution[,3]==0.9))
rep.folder = list("simulation_BM",
                  "simulation_contribution",
                  "simulation_Epsilon",
                  "simulation_contribution")
rep.prefix = list("bac.simulation.BM",
                  "bac.simulation.simulation.PE",
                  "bac.simulation.simulation.epsilon",
                  "bac.simulation.simulation.PE")
rep.pars = list(pars.bm,pars.contribution,pars.epsilon,pars.contribution)
trait.label = c(expression("BM"),expression("BM+"*epsilon),expression("BM+PE"),expression("BM+PE+"*epsilon))
# read in expected statistics
stat.expectation = readRDS("stat.expectation.RDS")
# read in data
all.asr = lapply(1:length(rep.file.idx),function(i){
  this.res = lapply(rep.file.idx[[i]],function(idx){
     f = paste0(rep.folder[[i]],"_",rep.pars[[i]]$batch[idx],"/",rep.prefix[[i]],".",rep.pars[[i]]$id[idx],".RDS")
     asr.res = readRDS(f)
     this.stats = lapply(asr.res,function(y){
       if(is.data.frame(y$CV)){
         cv.error = unname(mean(abs(y$trait[y$FMR$phy$tip.label]-y$CV$x)))
         cv.bias = unname(mean(y$trait[y$FMR$phy$tip.label]-y$CV$x))
         cv.pz = pseudo.Z.score(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                              error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE)
         asr.error = mean(abs(y$trait[-(1:Ntip(y$FMR$phy))]-y$ASR$x))
         asr.bias = mean(y$trait[-(1:Ntip(y$FMR$phy))]-y$ASR$x)
         asr.pz = pseudo.Z.score(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                              error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE)
         cv.rss = analyzeResidualErrorByCI(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                                    error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)
         cv.h = unname(analyzeResidualErrorByHeteroscedasticity(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                                                    error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1])
         asr.rss = analyzeResidualErrorByCI(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                                  error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)
         asr.h = unname(analyzeResidualErrorByHeteroscedasticity(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                                                  error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1])
         cv.d = ks.test(cv.pz$q,pnorm)$statistic
         asr.d = ks.test(asr.pz$q,pnorm)$statistic
       }else{
        cv.error = unname(mean(abs(y$trait[y$FMR$phy$tip.label]-y$CV$cv$x)))
        cv.bias = unname(mean(y$trait[y$FMR$phy$tip.label]-y$CV$cv$x))
        cv.pz = pseudo.Z.score(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                               error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE)
        asr.error = mean(abs(y$trait[-(1:Ntip(y$FMR$phy))]-y$ASR$ace$x))
        asr.bias = mean(y$trait[-(1:Ntip(y$FMR$phy))]-y$ASR$ace$x)
        asr.pz = pseudo.Z.score(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                                error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE)
        cv.rss = analyzeResidualErrorByCI(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                                          error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)
        cv.h = unname(analyzeResidualErrorByHeteroscedasticity(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                                                               error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1])
        asr.rss = analyzeResidualErrorByCI(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                                           error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)
        asr.h = unname(analyzeResidualErrorByHeteroscedasticity(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                                                                error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1])
        cv.d = ks.test(cv.pz$q,pnorm)$statistic
        asr.d = ks.test(asr.pz$q,pnorm)$statistic
       }
       return(list(cv.error=cv.error,cv.bias=cv.bias,
                   asr.error=asr.error,asr.bias=asr.bias,
                   cv.rss=cv.rss,cv.h=cv.h,
                   asr.rss=asr.rss,asr.h=asr.h,
                   cv.d=cv.d,asr.d=asr.d))
     })
  })
  return(this.res)
})
## Cross-validation statistics
tab1 = c(list(data.frame(error = mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.error})})),
                       error.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.error})})),
                       bias = mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.bias})})),
                       bias.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.bias})})),
                       D = mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.d})}))/mean(stat.expectation$cv.D),
                       D.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.d})}))/mean(stat.expectation$cv.D),
                       h = mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.h})}))/mean(stat.expectation$cv.h),
                       h.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.h})}))/mean(stat.expectation$cv.h),
                       rss= mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.rss})}))/mean(stat.expectation$cv.rss),
                       rss.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.rss})}))/mean(stat.expectation$cv.rss)
                       )),
            lapply(all.asr[-1],function(list.by.model){
              data.frame(error = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.error})}),1,mean),
                         error.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.error})}),1,std_err),
                         bias = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.bias})}),1,mean),
                         bias.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.bias})}),1,std_err),
                         D = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.d})}),1,mean)/mean(stat.expectation$cv.D),
                         D.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.d})}),1,std_err)/mean(stat.expectation$cv.D),
                         h = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.h})}),1,mean)/mean(stat.expectation$cv.h),
                         h.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.h})}),1,std_err)/mean(stat.expectation$cv.h),
                         rss= apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.rss})}),1,mean)/mean(stat.expectation$cv.rss),
                         rss.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.rss})}),1,std_err)/mean(stat.expectation$cv.rss)
                         )
              }))
tab1.p.value = c(list(t.test(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$cv.bias})}))$p.value),
                 lapply(all.asr[-1],function(list.by.model){
                   data.frame(error = t.test(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.error})})[2,]-
                                    sapply(list.by.model,function(x){sapply(x,function(y){y$cv.error})})[1,])$p.value,
                              bias = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.bias})}),1,
                                           function(z){t.test(z)$p.value}),
                              D = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.d})}),1,
                                        function(z){t.test(z,stat.expectation$cv.D,alternative = "greater")$p.value}),
                              h = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.h})}),1,
                                        function(z){t.test(z,stat.expectation$cv.h,alternative = "greater")$p.value}),
                              rss= apply(sapply(list.by.model,function(x){sapply(x,function(y){y$cv.rss})}),1,
                                         function(z){t.test(z,stat.expectation$cv.rss,alternative = "greater")$p.value}))
}))
## Global ancestral state statistics
tabs1 = c(list(data.frame(error = mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.error})})),
                         error.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.error})})),
                         bias = mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.bias})})),
                         bias.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.bias})})),
                         D = mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.d})}))/mean(stat.expectation$asr.D),
                         D.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.d})}))/mean(stat.expectation$asr.D),
                         h = mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.h})}))/mean(stat.expectation$asr.h),
                         h.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.h})}))/mean(stat.expectation$asr.h),
                         rss= mean(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.rss})}))/mean(stat.expectation$asr.rss),
                         rss.se = std_err(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.rss})}))/mean(stat.expectation$asr.rss)
)),
lapply(all.asr[-1],function(list.by.model){
  data.frame(error = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.error})}),1,mean),
             error.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.error})}),1,std_err),
             bias = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.bias})}),1,mean),
             bias.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.bias})}),1,std_err),
             D = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.d})}),1,mean)/mean(stat.expectation$asr.D),
             D.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.d})}),1,std_err)/mean(stat.expectation$asr.D),
             h = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.h})}),1,mean)/mean(stat.expectation$asr.h),
             h.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.h})}),1,std_err)/mean(stat.expectation$asr.h),
             rss= apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.rss})}),1,mean)/mean(stat.expectation$asr.rss),
             rss.se = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.rss})}),1,std_err)/mean(stat.expectation$asr.rss)
  )
}))
tabs1.p.value = c(list(t.test(sapply(all.asr[[1]],function(x){sapply(x,function(y){y$asr.bias})}))$p.value),
                 lapply(all.asr[-1],function(list.by.model){
                   data.frame(error = t.test(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.error})})[2,]-
                                               sapply(list.by.model,function(x){sapply(x,function(y){y$asr.error})})[1,])$p.value,
                              bias = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.bias})}),1,
                                           function(z){t.test(z)$p.value}),
                              D = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.d})}),1,
                                        function(z){t.test(z,stat.expectation$asr.D,alternative = "greater")$p.value}),
                              h = apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.h})}),1,
                                        function(z){t.test(z,stat.expectation$asr.h,alternative = "greater")$p.value}),
                              rss= apply(sapply(list.by.model,function(x){sapply(x,function(y){y$asr.rss})}),1,
                                         function(z){t.test(z,stat.expectation$asr.rss,alternative = "greater")$p.value}))
                 }))
#
ref.size = 0.6
line.size = 0.6
# Figure 1 data
BM.Felsenstein.pZ = lapply(1:length(rep.file.idx),function(i){
  f = paste0(rep.folder[[i]],"_",rep.pars[[i]]$batch[rep.file.idx[[i]][1]],"/",rep.prefix[[i]],".",rep.pars[[i]]$id[rep.file.idx[[i]][1]],".RDS")
  x = readRDS(f)
  Z = (x[[1]]$trait[x[[1]]$FMR$phy$tip.label]-x[[1]]$CV$x)/sqrt(x[[1]]$CV$var)
  pZ.df = data.frame(p=pnorm(Z),Z=Z,var = x[[1]]$CV$var,se = sqrt(x[[1]]$CV$var))
})
# Figure 1 plots
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
    ggtitle(trait.label[i])+xlab("Predicted SE")+ylab("Z-score")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
BM.Felsenstein.pp.plot = lapply(1:length(BM.Felsenstein.pZ),function(i){
  df = BM.Felsenstein.pZ[[i]]
  this.df = data.frame(empirical = (1:dim(df)[1])/(dim(df)[1]+1),theoretical=sort(df$p,decreasing = FALSE))
  ggplot()+geom_line(mapping = aes(x=theoretical,y=empirical),data=this.df,alpha=0.8,size=line.size)+
    geom_abline(slope = 1,intercept = 0,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(xlim = c(0,1),ylim=c(0,1))+
    scale_x_continuous(breaks=seq(0,1,0.25),
                       labels = seq(0,1,0.25),expand = c(0,0.02))+
    scale_y_continuous(breaks=seq(0,1,0.25),
                       labels = seq(0,1,0.25),expand = c(0,0.02))+
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
    scale_y_continuous(breaks = logitSigmoid(x = c(0.3,0.75,0.9,0.95,0.97,0.99,1),p = 0.95),
                       labels = c(0.3,0.75,0.9,0.95,0.97,0.99,1))+
    xlab("Predicted SE")+ylab("Coverage probability")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
BM.Felsenstein.combined.plot = ggarrange(plotlist = do.call(c,lapply(1:length(BM.Felsenstein.pZ),function(i){
  list(BM.Felsenstein.Z.plot[[i]],BM.Felsenstein.qq.plot[[i]],BM.Felsenstein.alpha.plot[[i]])
})),
ncol = 3,nrow = 4,labels = "AUTO",align = "hv")
#
ggsave(filename = "Fig1.pdf",device = "pdf",plot = BM.Felsenstein.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 4.4,units = "in")
ggsave(filename = "Fig1.png",device = "png",plot = BM.Felsenstein.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 4.4,units = "in",dpi = "print",type="cairo")
#
# Figure 3 data
PEBME.pZ = lapply(2:length(rep.file.idx),function(i){
  f = paste0(rep.folder[[i]],"_",rep.pars[[i]]$batch[rep.file.idx[[i]][1]],"/",rep.prefix[[i]],".",rep.pars[[i]]$id[rep.file.idx[[i]][1]],".RDS")
  x = readRDS(f)
  pZ = pseudo.Z.score(trait = x[[2]]$trait[x[[2]]$FMR$phy$tip.label],pred = x[[2]]$CV$cv$x,error = x[[2]]$CV$error,
                      epsilon = x[[2]]$FMR$params$epsilon,laplace = FALSE)
  pZ.df = data.frame(p=pZ$p,Z=pZ$q,var = x[[2]]$CV$cv$var,se = sqrt(x[[2]]$CV$cv$var))
})
# Figure 3 plots
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
    coord_cartesian(xlim = c(0,1),ylim=c(0,1))+
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
    scale_y_continuous(breaks = logitSigmoid(x = c(0.3,0.75,0.9,0.95,0.97,0.99,1),p = 0.95),
                       labels = c(0.3,0.75,0.9,0.95,0.97,0.99,1))+
    xlab("Predicted SE")+ylab("Coverage probability")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
PEBME.combined.plot = ggarrange(plotlist = do.call(c,lapply(1:length(PEBME.pZ),function(i){
  list(PEBME.Z.plot[[i]],PEBME.qq.plot[[i]],PEBME.alpha.plot[[i]])
})),
ncol = 3,nrow = 3,labels = "AUTO",align = "hv")
#
ggsave(filename = "Fig3.pdf",device = "pdf",plot = PEBME.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 3.25,units = "in")
ggsave(filename = "Fig3.png",device = "png",plot = PEBME.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 3.25,units = "in",dpi = "print",type="cairo")
#

# Figure S1 data
ASR.BM.Felsenstein.pZ = lapply(1:length(rep.file.idx),function(i){
  f = paste0(rep.folder[[i]],"_",rep.pars[[i]]$batch[rep.file.idx[[i]][1]],"/",rep.prefix[[i]],".",rep.pars[[i]]$id[rep.file.idx[[i]][1]],".RDS")
  x = readRDS(f)
  Z = (x[[1]]$trait[-(1:Ntip(x[[1]]$FMR$phy))]-x[[1]]$ASR$x)/sqrt(x[[1]]$ASR$var)
  pZ.df = data.frame(p=pnorm(Z),Z=Z,var = x[[1]]$ASR$var,se = sqrt(x[[1]]$ASR$var))
})
# Figure S1 plots
ASR.BM.Felsenstein.Z.plot = lapply(1:length(ASR.BM.Felsenstein.pZ),function(i){
  df = ASR.BM.Felsenstein.pZ[[i]]
  ggplot()+geom_point(mapping = aes(x=se,y = Z),data = df,alpha=0.5,shape=".")+
    geom_hline(yintercept = 2,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    geom_hline(yintercept = -2,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(ylim = c(-20,20),
                    xlim=c(10^(floor(2*log10(min(df$se)))/2),10^(ceiling(2*log10(max(df$se)))/2)))+
    scale_x_continuous(trans = "log10",breaks=c(1e-4,0.001,0.01,0.1,1,10),
                       labels = c(expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
    scale_y_continuous(breaks = seq(-20,20,4))+
    ggtitle(trait.label[i])+xlab("Predicted SE")+ylab("Z-score")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
ASR.BM.Felsenstein.pp.plot = lapply(1:length(ASR.BM.Felsenstein.pZ),function(i){
  df = ASR.BM.Felsenstein.pZ[[i]]
  this.df = data.frame(empirical = (1:dim(df)[1])/(dim(df)[1]+1),theoretical=sort(df$p,decreasing = FALSE))
  ggplot()+geom_line(mapping = aes(x=theoretical,y=empirical),data=this.df,alpha=0.8,size=line.size)+
    geom_abline(slope = 1,intercept = 0,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    coord_cartesian(xlim = c(0,1),ylim=c(0,1))+
    scale_x_continuous(breaks=seq(0,1,0.25),
                       labels = seq(0,1,0.25),expand = c(0,0.02))+
    scale_y_continuous(breaks=seq(0,1,0.25),
                       labels = seq(0,1,0.25),expand = c(0,0.02))+
    xlab("Theoretical probability")+ylab("Empirical probability")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
ASR.BM.Felsenstein.qq.plot = lapply(1:length(ASR.BM.Felsenstein.pZ),function(i){
  df = ASR.BM.Felsenstein.pZ[[i]]
  this.df = data.frame(empirical = sort(df$Z,decreasing = FALSE),theoretical= qnorm((1:dim(df)[1])/(dim(df)[1]+1)))
  ggplot()+geom_line(mapping = aes(x=theoretical,y=empirical),data=this.df,alpha=0.8,size=line.size)+
    geom_abline(slope = 1,intercept = 0,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    xlab("Theoretical quantile")+ylab("Empirical quantile")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
ASR.BM.Felsenstein.alpha.plot = lapply(1:length(ASR.BM.Felsenstein.pZ),function(i){
  df = ASR.BM.Felsenstein.pZ[[i]]
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
    scale_y_continuous(breaks = logitSigmoid(x = c(0.3,0.75,0.9,0.95,0.97,0.99,1),p = 0.95),
                       labels = c(0.3,0.75,0.9,0.95,0.97,0.99,1))+
    xlab("Predicted SE")+ylab("Coverage probability")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
ASR.BM.Felsenstein.combined.plot = ggarrange(plotlist = do.call(c,lapply(1:length(ASR.BM.Felsenstein.pZ),function(i){
  list(ASR.BM.Felsenstein.Z.plot[[i]],ASR.BM.Felsenstein.qq.plot[[i]],ASR.BM.Felsenstein.alpha.plot[[i]])
})),
ncol = 3,nrow = 4,labels = "AUTO",align = "hv")
#
ggsave(filename = "FigS1.pdf",device = "pdf",plot = ASR.BM.Felsenstein.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 4.4,units = "in")
ggsave(filename = "FigS1.png",device = "png",plot = ASR.BM.Felsenstein.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 4.4,units = "in",dpi = "print",type="cairo")
#
# Figure S3
ASR.PEBME.pZ = lapply(2:length(rep.file.idx),function(i){
  f = paste0(rep.folder[[i]],"_",rep.pars[[i]]$batch[rep.file.idx[[i]][1]],"/",rep.prefix[[i]],".",rep.pars[[i]]$id[rep.file.idx[[i]][1]],".RDS")
  x = readRDS(f)
  pZ = pseudo.Z.score(trait = x[[2]]$trait[-(1:Ntip(x[[2]]$FMR$phy))],pred = x[[2]]$ASR$ace$x,error = x[[2]]$ASR$error,
                      epsilon = x[[2]]$FMR$params$epsilon,laplace = FALSE)
  pZ.df = data.frame(p=pZ$p,Z=pZ$q,var = x[[2]]$ASR$ace$var,se = sqrt(x[[2]]$ASR$ace$var))
})
#
ASR.PEBME.Z.plot = lapply(1:length(ASR.PEBME.pZ),function(i){
  df = ASR.PEBME.pZ[[i]]
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
ASR.PEBME.pp.plot = lapply(1:length(ASR.PEBME.pZ),function(i){
  df = ASR.PEBME.pZ[[i]]
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
ASR.PEBME.qq.plot = lapply(1:length(ASR.PEBME.pZ),function(i){
  df = ASR.PEBME.pZ[[i]]
  this.df = data.frame(empirical = sort(df$Z,decreasing = FALSE),theoretical= qnorm((1:dim(df)[1])/(dim(df)[1]+1)))
  ggplot()+geom_line(mapping = aes(x=theoretical,y=empirical),data=this.df,alpha=0.8,size=line.size)+
    geom_abline(slope = 1,intercept = 0,col="red",size=ref.size,linetype="dashed",alpha=0.6)+
    xlab("Theoretical quantile")+ylab("Empirical quantile")+#ggtitle(trait.label[i+1])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
ASR.PEBME.alpha.plot = lapply(1:length(ASR.PEBME.pZ),function(i){
  df = ASR.PEBME.pZ[[i]]
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
    scale_y_continuous(breaks = logitSigmoid(x = c(0.3,0.75,0.9,0.95,0.97,0.99,1),p = 0.95),
                       labels = c(0.3,0.75,0.9,0.95,0.97,0.99,1))+
    xlab("Predicted SE")+ylab("Coverage probability")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line())
})
ASR.PEBME.combined.plot = ggarrange(plotlist = do.call(c,lapply(1:length(ASR.PEBME.pZ),function(i){
  list(ASR.PEBME.Z.plot[[i]],ASR.PEBME.qq.plot[[i]],ASR.PEBME.alpha.plot[[i]])
})),
ncol = 3,nrow = 3,labels = "AUTO",align = "hv")
#
ggsave(filename = "FigS3.pdf",device = "pdf",plot = ASR.PEBME.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 3.25,units = "in")
ggsave(filename = "FigS3.png",device = "png",plot = ASR.PEBME.combined.plot,path = "./",
       scale = 2,width = 3.25,height = 3.25,units = "in",dpi = "print",type="cairo")
#
save.image("simulation.summary.RData")
