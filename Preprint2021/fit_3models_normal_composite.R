# setup the environment for model fitting
source("fitPE_functions.R")

# read in key arguments for the fitting
arg=commandArgs(TRUE)
# interpret the arguments
# get data
all.data = list(euka=readRDS(arg[1])[readRDS(arg[2])])
all.trait = unlist(lapply(all.data,function(x){lapply(x,function(y){y$dat})}))
all.clade.taxa = names(unlist(lapply(all.data,function(x){lapply(x,function(y){""})})))
#
pseudo.tree = rtree(length(all.clade.taxa),tip.label = paste0("root_",all.clade.taxa),br = 1e4)
pseudo.tree = makeNodeLabel(phy = pseudo.tree,prefix = "PSEUDO_NODE_")
tree = pseudo.tree
for(clade in names(all.data)){
  clade.data = all.data[[clade]]
  for(taxa in names(clade.data)){
    taxa.tree = clade.data[[taxa]]$phy
    taxa.tree$tip.label = paste0(clade,".",taxa,".",taxa.tree$tip.label)
    tree = bind.tree(x = tree,y = taxa.tree,where = which(tree$tip.label==paste0("root_",clade,".",taxa)))
  }
}
trait = all.trait
#
valid.node = which(is.na(tree$node.label))
#
# set up a name for saving the result of this analysis
job.name = arg[3]
if(is.na(job.name)) stop("Job name is required and cannot be omitted!")
# how many cores to use in the fitting
numCores=ceiling(as.numeric(arg[4]))
if(is.na(numCores)) numCores=1
# how far randomized starting points can go from the initial estimate
repRange = as.numeric(arg[5])
if(is.na(repRange)) repRange =10
# set manual fold-constraint between sizes
size.ratio = as.numeric(arg[6])
if(is.na(size.ratio)) size.ratio=10
#
ignore.zero = (as.numeric(arg[7])>0)
if(is.na(ignore.zero)) ignore.zero = FALSE
#
AIC.cutoff = as.numeric(arg[8])
if(is.na(AIC.cutoff)) AIC.cutoff=2
# what's the tolerance for coarse-grind model fitting?
reltol = as.numeric(arg[9])
if(is.na(reltol)) reltol = 1e-4
# set manual scaling coefficient
scale.coef = as.numeric(arg[10])
# organize data for easy transfer to other analysis
this.data = list(phy=tree,dat=trait[tree$tip.label])
if(any(is.na(this.data$dat))) stop("NAs are not supported in tip values!")
saveRDS(this.data,sprintf("%s.data.RDS",job.name))
# estimating scaling coefficient, epsilon and sigma term so the optimizer have a reasonable starting point
raw.contrast = pic(x = trait[tree$tip.label],phy=tree,scaled = FALSE,var.contrasts = TRUE)
if(is.na(scale.coef)) scale.coef = 2/unname(var(raw.contrast[raw.contrast[,2]<=quantile(raw.contrast[,2],probs=0.05),1]))
sprintf("The scaling coefficient is %.3f",scale.coef)
est.sigma = var(pic(x = trait[tree$tip.label],phy=tree))*scale.coef
median.l = median(tree$edge.length)
#
terminal.pic = findTerminalPIC(phy = tree,x = sqrt(scale.coef)*trait[tree$tip.label],remove.zero = TRUE)
term.pebme = fitPE(x = sapply(terminal.pic$x,function(x){x[1]-x[2]}),l = sapply(terminal.pic$l,sum),numCores = numCores,laplace = FALSE,
                   start.value = list(lambda=log(1/median.l),size=log(0.9*est.sigma*median.l),sigma=log(0.1*est.sigma),epsilon=log(1)))
# fitting the models
model.functions = c(bme="fitPE.phy",pebme ="fitPE.phy")
fixed.params = list(bme=list(lambda=log(0),size=log(0)),pebme =list())
initial.parameters = list(bme=c(lambda=0,size=0,sigma=est.sigma,epsilon=1),pebme = term.pebme$params)
ini.multiplier = matrix(data = c(1,1,repRange,repRange,repRange,1/repRange,1/repRange,repRange,1/repRange,1/repRange),nrow = 5,ncol = 2,byrow = TRUE)
#
res = list(bme=list(),pebme =list())
rep.idx = c(1,5)
for(i in rev(1:2)){
  print(names(res)[i])
  for(rr in 1:rep.idx[i]){
    print(rr)
    ini.func = "initialize.parameters.SP"
    start.params = initialize.parameters.SP(initial.parameters[[i]]*
                                              c(ini.multiplier[rr,1],ini.multiplier[rr,2],1,1))
    if(any(sapply(start.params,is.nan))) next
    try({
      mdl= fit.model.byNM.2step(fit.func = model.functions[[i]],
                                start.params = list(start.value=start.params),
                                params1 =list(fixed=fixed.params[[i]],size.ratio=size.ratio,ignore.zero=ignore.zero,laplace = FALSE,
                                              phy=tree,x=sqrt(scale.coef)*trait[tree$tip.label],numCores=numCores,reltol=reltol),
                                params2 =list(fixed=fixed.params[[i]],size.ratio=size.ratio,ignore.zero=ignore.zero,laplace = FALSE,
                                              phy=tree,x=sqrt(scale.coef)*trait[tree$tip.label],numCores=numCores),
                                AIC.func = "get.AIC",
                                params.func = "get.params",
                                ini.func = ini.func,
                                AIC.cutoff = AIC.cutoff)
      res[[i]][[rr]] = mdl
    })
    #
    saveRDS(list(scale=scale.coef,res = res),
            sprintf("%s.models.RDS",job.name))
  }
}
#
model.summary = lapply(res,function(x){
  best.rep = which.min(sapply(x,function(y){if(is.null(y)){
    return(NA)
  }else{
    return(y$low$AIC)
  }}))
  best.AIC = x[[best.rep]]$low$AIC
  best.params = x[[best.rep]]$low$params
  return(list(AIC=best.AIC,params=best.params))
})
#
saveRDS(list(scale=scale.coef,jitter=jitter.size,jitter.trait = trait,res = res,stat=sapply(model.summary,get.AIC),params=lapply(model.summary,get.params)),
        sprintf("%s.models.RDS",job.name))
#
warnings()
