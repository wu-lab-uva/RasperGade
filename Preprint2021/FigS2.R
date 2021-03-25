#
setwd("F:/Dropbox/Scripts/16S_Copy_Number/Version3/")
source("fitPE_functions.R")
setwd("F:/Dropbox/Scripts/16S_Copy_Number/Version4/")
source("RasperGade_reconstruction.R")
source("RasperGade_diagnosis.R")
setwd("F:/Dropbox/Scripts/16S_Copy_Number/20210302/")
library(ggplot2)
library(ggpubr)
# parse simulation results
epsilons=c(0,1,10)
contributions = rev(c(0,0.1,0.3,0.5,0.7,0.9,1))
replicates = 1:100
batch = 1
pars.epsilon = expand.grid(epsilon=epsilons,lambda=5,contribution=0.9,rep=replicates,batch=batch,KEEP.OUT.ATTRS = FALSE)
pars.epsilon$id = rep(1:(dim(pars.epsilon)[1]/length(unique(pars.epsilon$batch))),length(unique(pars.epsilon$batch)))
pars.contribution = expand.grid(epsilon=10,lambda=5,contribution=contributions,rep=replicates,batch=batch,KEEP.OUT.ATTRS = FALSE)
pars.contribution$id = rep(1:(dim(pars.contribution)[1]/length(unique(pars.contribution$batch))),length(unique(pars.contribution$batch)))

# read in expected statistics
stat.expectation = readRDS("stat.expectation.RDS")
# epsilon results
epsilon.df = do.call(rbind,lapply(1:dim(pars.epsilon)[1],function(i){
  asr.res = readRDS(paste0("simulation_Epsilon_only/","/bac.simulation.simulation.epsilon.only.",pars.epsilon$id[i],".RDS"))
  this.stats = sapply(asr.res,function(y){
    if(is.data.frame(y$CV)){
      ers = pseudo.Z.score(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                           error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE)
      c(analyzeResidualErrorByCI(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                                     error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20),
        analyzeResidualErrorByHeteroscedasticity(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                                                 error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1],
        1-sum(sapply(ers$p,function(x){min(x,1-x)<=0.025}))/length(ers$p),
        ks.test(ers$q,pnorm)$statistic)
    }else{
      ers = pseudo.Z.score(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                           error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE)
      c(analyzeResidualErrorByCI(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                                     error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20),
        analyzeResidualErrorByHeteroscedasticity(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                                                 error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1],
        1-sum(sapply(ers$p,function(x){min(x,1-x)<=0.025}))/length(ers$p),
        ks.test(ers$q,pnorm)$statistic)
    }
  })
  return(data.frame(rss=this.stats[1,],h=this.stats[2,],cp95=this.stats[3,],D=this.stats[4,],model=c("BM","True model")))
}))
epsilon.ace.df = do.call(rbind,lapply(1:dim(pars.epsilon)[1],function(i){
  asr.res = readRDS(paste0("simulation_Epsilon_only","/bac.simulation.simulation.epsilon.only.",pars.epsilon$id[i],".RDS"))
  this.stats = sapply(asr.res,function(y){
    if(is.data.frame(y$ASR)){
      ers = pseudo.Z.score(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                           error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE)
      c(analyzeResidualErrorByCI(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                                 error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20),
        analyzeResidualErrorByHeteroscedasticity(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                                                 error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1],
        1-sum(sapply(ers$p,function(x){min(x,1-x)<=0.025}))/length(ers$p),
        ks.test(ers$q,pnorm)$statistic)
    }else{
      ers = pseudo.Z.score(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                           error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE)
      c(analyzeResidualErrorByCI(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                                 error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20),
        analyzeResidualErrorByHeteroscedasticity(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                                                 error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1],
        1-sum(sapply(ers$p,function(x){min(x,1-x)<=0.025}))/length(ers$p),
        ks.test(ers$q,pnorm)$statistic)
    }
  })
  return(data.frame(rss = this.stats[1,],h=this.stats[2,],cp95=this.stats[3,],D=this.stats[4,],model=c("BM","True model")))
}))
# contribution results
contribution.df = do.call(rbind,lapply(1:dim(pars.contribution)[1],function(i){
  asr.res = readRDS(paste0("simulation_PEonly","/bac.simulation.simulation.PEonly.",pars.contribution$id[i],".RDS"))
  this.stats = sapply(asr.res,function(y){
    if(is.data.frame(y$CV)){
      ers = pseudo.Z.score(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                           error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE)
      c(analyzeResidualErrorByCI(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                                     error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20),
        analyzeResidualErrorByHeteroscedasticity(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$x,
                                                 error = y$CV$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1],
        1-sum(sapply(ers$p,function(x){min(x,1-x)<=0.025}))/length(ers$p),
        ks.test(ers$q,pnorm)$statistic)
    }else{
      ers = pseudo.Z.score(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                           error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE)
      c(analyzeResidualErrorByCI(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                                     error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20),
        analyzeResidualErrorByHeteroscedasticity(trait = y$trait[y$FMR$phy$tip.label],pred = y$CV$cv$x,
                                                 error = y$CV$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1],
        1-sum(sapply(ers$p,function(x){min(x,1-x)<=0.025}))/length(ers$p),
        ks.test(ers$q,pnorm)$statistic)
    }
  })
  return(data.frame(rss=this.stats[1,],h=this.stats[2,],cp95=this.stats[3,],D=this.stats[4,],model=c("BM","True model")))
}))
contribution.ace.df = do.call(rbind,lapply(1:dim(pars.contribution)[1],function(i){
  asr.res = readRDS(paste0("simulation_PEonly","/bac.simulation.simulation.PEonly.",pars.contribution$id[i],".RDS"))
  this.stats = sapply(asr.res,function(y){
    if(is.data.frame(y$ASR)){
      ers = pseudo.Z.score(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                           error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE)
      c(analyzeResidualErrorByCI(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                                 error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20),
        analyzeResidualErrorByHeteroscedasticity(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$x,
                                                 error = y$ASR$var,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1],
        1-sum(sapply(ers$p,function(x){min(x,1-x)<=0.025}))/length(ers$p),
        ks.test(ers$q,pnorm)$statistic)
    }else{
      ers = pseudo.Z.score(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                           error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE)
      c(analyzeResidualErrorByCI(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                                 error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20),
        analyzeResidualErrorByHeteroscedasticity(trait = y$trait[-(1:Ntip(y$FMR$phy))],pred = y$ASR$ace$x,
                                                 error = y$ASR$error,epsilon = y$FMR$params$epsilon,laplace = FALSE,bin = 20)[1],
        1-sum(sapply(ers$p,function(x){min(x,1-x)<=0.025}))/length(ers$p),
        ks.test(ers$q,pnorm)$statistic)
    }
  })
  return(data.frame(rss = this.stats[1,],h=this.stats[2,],cp95=this.stats[3,],D=this.stats[4,],model=c("BM","True model")))
}))

# Figure S2 row 1
plot.alpha=1
dot.size = 0.5
{
  epsilon.summary = epsilon.df
  epsilon.summary$epsilon = rep(pars.epsilon$epsilon,each=2)
  epsilon.summary$epsilon[epsilon.summary$epsilon==0] = 0.1
  epsilon.summary$epsilon = epsilon.summary$epsilon*exp(rnorm(n = dim(epsilon.summary)[1],sd = 0.1))
  epsilon.summary$h = epsilon.summary$h/mean(stat.expectation$cv.h)
  epsilon.summary$D = epsilon.summary$D/mean(stat.expectation$cv.D)
  epsilon.summary$rss = epsilon.summary$rss/mean(stat.expectation$cv.rss)
  epsilon.D.plot = ggplot()+geom_point(mapping = aes(x=epsilon,y=D,color=model),data=epsilon.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(trans="log10",breaks=c(0.1,1,10),labels = c("0",expression(4%*%10^-4),expression(4%*%10^-3)))+
    scale_y_continuous(trans="sqrt",breaks = c(1,2.5,5,7.5,10))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab(expression(epsilon))+ylab("D")+ggtitle("Time-independent variation","(hidden states)" )+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
  epsilon.h.plot = ggplot()+geom_point(mapping = aes(x=epsilon,y=h,color=model),data=epsilon.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(trans="log10",breaks=c(0.1,1,10),labels = c("0",expression(4%*%10^-4),expression(4%*%10^-3)))+
    scale_y_continuous(trans="sqrt",breaks = c(1,10,25,50,80,120))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab(expression(epsilon))+ylab(expression(median*" "*chi^2))+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
  epsilon.er.plot = ggplot()+geom_point(mapping = aes(x=epsilon,y=rss,color=model),data=epsilon.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(trans="log10",breaks=c(0.1,1,10),labels = c("0",expression(4%*%10^-4),expression(4%*%10^-3)))+
    scale_y_continuous(trans="sqrt",breaks = c(1,10,20,30,40))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab(expression(epsilon))+ylab(expression(RSS[CP]))+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))  
}

# row 2
{
  epsilon.ace.summary = epsilon.ace.df
  epsilon.ace.summary$epsilon = rep(pars.epsilon$epsilon,each=2)
  epsilon.ace.summary$epsilon[epsilon.ace.summary$epsilon==0] = 0.1
  epsilon.ace.summary$epsilon = epsilon.ace.summary$epsilon*exp(rnorm(n = dim(epsilon.ace.summary)[1],sd = 0.1))
  epsilon.ace.summary$h = epsilon.ace.summary$h/mean(stat.expectation$asr.h)
  epsilon.ace.summary$D = epsilon.ace.summary$D/mean(stat.expectation$asr.D)
  epsilon.ace.summary$rss = epsilon.ace.summary$rss/mean(stat.expectation$asr.rss)
  epsilon.ace.D.plot = ggplot()+geom_point(mapping = aes(x=epsilon,y=D,color=model),data=epsilon.ace.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(trans="log10",breaks=c(0.1,1,10),labels = c("0",expression(4%*%10^-4),expression(4%*%10^-3)))+
    scale_y_continuous(trans="sqrt",breaks = c(1,2,3,4))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab(expression(epsilon))+ylab("D")+ggtitle("Time-independent variation","(ancestral states)")+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
  epsilon.ace.h.plot = ggplot()+geom_point(mapping = aes(x=epsilon,y=h,color=model),data=epsilon.ace.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(trans="log10",breaks=c(0.1,1,10),labels = c("0",expression(4%*%10^-4),expression(4%*%10^-3)))+
    scale_y_continuous(trans="sqrt",breaks = c(1,10,25,50,80))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab(expression(epsilon))+ylab(expression(median*" "*chi^2))+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
  epsilon.ace.er.plot = ggplot()+geom_point(mapping = aes(x=epsilon,y=rss,color=model),data=epsilon.ace.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(trans="log10",breaks=c(0.1,1,10),labels = c("0",expression(4%*%10^-4),expression(4%*%10^-3)))+
    scale_y_continuous(trans="sqrt",breaks = c(1,10,20))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab(expression(epsilon))+ylab(expression(RSS[CP]))+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
}

# row 3
{
  contribution.summary = contribution.df
  contribution.summary$contribution = rep(pars.contribution$contribution,each=2)
  contribution.summary$contribution = contribution.summary$contribution+rnorm(n = dim(contribution.summary)[1],sd = 0.01)
  contribution.summary$h = contribution.summary$h/mean(stat.expectation$cv.h)
  contribution.summary$D = contribution.summary$D/mean(stat.expectation$cv.D)
  contribution.summary$rss = contribution.summary$rss/mean(stat.expectation$cv.rss)
  contribution.D.plot = ggplot()+geom_point(mapping = aes(x=contribution,y=D,color=model),data=contribution.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))+
    scale_y_continuous(trans="sqrt",breaks = c(1,10,20,30,40))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab("PE contribution")+ylab("D")+ggtitle("Pulsed evolution","(hidden states)")+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
  contribution.h.plot = ggplot()+geom_point(mapping = aes(x=contribution,y=h,color=model),data=contribution.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))+
    scale_y_continuous(trans="sqrt",breaks=c(1,10,25,50,80,120))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab("PE contribution")+ylab(expression(median*" "*chi^2))+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
  contribution.er.plot = ggplot()+geom_point(mapping = aes(x=contribution,y=rss,color=model),data=contribution.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))+
    scale_y_continuous(trans="sqrt",breaks = c(1,10,20,30))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab("PE contribution")+ylab(expression(RSS[CP]))+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
}

# row 4
{
  contribution.ace.summary = contribution.ace.df
  contribution.ace.summary$contribution = rep(pars.contribution$contribution,each=2)
  contribution.ace.summary$contribution = contribution.ace.summary$contribution+rnorm(n = dim(contribution.ace.summary)[1],sd = 0.01)
  contribution.ace.summary$h = contribution.ace.summary$h/mean(stat.expectation$asr.h)
  contribution.ace.summary$D = contribution.ace.summary$D/mean(stat.expectation$asr.D)
  contribution.ace.summary$rss = contribution.ace.summary$rss/mean(stat.expectation$asr.rss)
  contribution.ace.D.plot = ggplot()+geom_point(mapping = aes(x=contribution,y=D,color=model),data=contribution.ace.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))+
    scale_y_continuous(trans="sqrt",breaks =c(1,5,10,15))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab("PE contribution")+ylab("D")+ggtitle("Pulsed evolution","(ancestral states)")+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
  contribution.ace.h.plot = ggplot()+geom_point(mapping = aes(x=contribution,y=h,color=model),data=contribution.ace.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))+
    scale_y_continuous(trans="sqrt",breaks = c(1,10,25,50,80))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab("PE contribution")+ylab(expression(median*" "*chi^2))+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
  contribution.ace.er.plot = ggplot()+geom_point(mapping = aes(x=contribution,y=rss,color=model),data=contribution.ace.summary,alpha=plot.alpha,size=dot.size)+
    geom_hline(yintercept = 1,linetype="dashed",color="red")+
    scale_x_continuous(breaks=seq(0,1,0.25),labels=c("0%","25%","50%","75%","100%"))+
    scale_y_continuous(trans="sqrt",breaks = c(1,4,8,12,16))+
    scale_color_manual(values = c(BM="black","True model"="#69b3a2"))+
    xlab("PE contribution")+ylab(expression(RSS[CP]))+
    guides(colour = guide_legend(override.aes = list(size=1.5)))+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(hjust = 1))
}
#
combined.trend.plot = ggarrange(plotlist = list(epsilon.D.plot,epsilon.h.plot,epsilon.er.plot,
                                                epsilon.ace.D.plot,epsilon.ace.h.plot,epsilon.ace.er.plot,
                                                contribution.D.plot,contribution.h.plot,contribution.er.plot,
                                                contribution.ace.D.plot,contribution.ace.h.plot,contribution.ace.er.plot),
                                ncol = 3,nrow = 4,align = "hv",labels = "AUTO",common.legend = TRUE,legend = "bottom",label.x = c(0,0,-0.05,0,0,-0.05))
#
ggsave(filename = "FigS2.pdf",device = "pdf",plot = combined.trend.plot,path = "./",
       scale = 3,width = 3.25,height = 4.5,units = "in")
ggsave(filename = "FigS2.png",device = "png",plot = combined.trend.plot,path = "./",
       scale = 3,width = 3.25,height = 4.5,units = "in",dpi = "print",type="cairo")
#
epsilon.improve = sapply(seq(1,dim(epsilon.summary)[1],2),function(i){1-(unlist(epsilon.summary[i+1,c(1,2,4)])-1)/(unlist(epsilon.summary[i,c(1,2,4)])-1)})
contribution.improve = sapply(seq(1,dim(contribution.summary)[1],2),function(i){1-(unlist(contribution.summary[i+1,c(1,2,4)])-1)/(unlist(contribution.summary[i,c(1,2,4)])-1)})
epsilon.ace.improve = sapply(seq(1,dim(epsilon.ace.summary)[1],2),function(i){1-(unlist(epsilon.ace.summary[i+1,c(1,2,4)])-1)/(unlist(epsilon.ace.summary[i,c(1,2,4)])-1)})
contribution.ace.improve = sapply(seq(1,dim(contribution.ace.summary)[1],2),function(i){1-(unlist(contribution.ace.summary[i+1,c(1,2,4)])-1)/(unlist(contribution.ace.summary[i,c(1,2,4)])-1)})
#
save.image("trend.summary.RData")


