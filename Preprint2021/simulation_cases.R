library(RasperGade)
# read in key arguments for the fitting
arg=commandArgs(TRUE)
# interpret the arguments
# read in phylogeny and trait from 1 RDS file
tree.name = arg[1]
if(is.na(tree.name)) stop("Tree is required and cannot be omitted!")
tree = read.tree(tree.name)
# how many cores to use in the fitting
numCores=as.numeric(arg[2])
# save file prefix
save.name = arg[3]
# test on normal distribution only
laplace=FALSE
#
epsilons= as.numeric(arg[4])
lambdas = as.numeric(arg[5])
contributions = as.numeric(arg[6])
if(any(is.na(c(epsilons,lambdas,contributions)))) stop("Model parameters cannot be NA")
replicates = 1:100
pars = expand.grid(epsilons,lambdas,contributions,replicates)
total.rate = 1
all.asr = mclapply(1:dim(pars)[1],function(i){
  model.params = c(lambda = pars[i,2],size = total.rate*pars[i,3]/pars[i,2],sigma = total.rate*(1-pars[i,3]),epsilon = total.rate*4e-4*pars[i,1])
  print(model.params)
  trait = simulatePET.phy(phy = tree,lambda = model.params[1],size = model.params[2],sigma = model.params[3],epsilon = 0)
  trait = trait +rnorm(n = length(trait))*sqrt(model.params[4]/2)
  #
  t0 = Sys.time()
  bm.rate = var(pic(x = trait[tree$tip.label],phy = tree))
  asr.res = list(list(FMR=list(phy=tree,params=list(sigma=bm.rate,epsilon=0)),trait=trait,
                      CV=do.call(rbind,lapply(1:Ntip(tree),function(i){
                        predTraitByPIC(phy = tree,trait = trait[tree$tip.label[-i]],rate=bm.rate)})),
                      ASR=continuousAceByPIC(x = trait[tree$tip.label],phy = tree,rate = bm.rate)))
  t1 = Sys.time()
  asr.res[[1]]$t= c(t0,t1)
  t0 = Sys.time()
  FMR = fullMarginalReconstructionWithPE(phy = tree,x = trait[tree$tip.label],params = model.params,
                                         laplace = laplace,approximate = 1,numCores = 1,asymptotic=5)
  CV = crossValidationWithPE(FMR = FMR,add.epsilon = FALSE,numCores = 1,laplace = laplace)
  ASR = globalReconstructionWithPE(FMR = FMR,add.epsilon = FALSE,numCores = 1,laplace = laplace)
  t1 = Sys.time()
  asr.res[[2]] = list(FMR=FMR,trait=trait,CV=CV,ASR=ASR,t = c(t0,t1))
  saveRDS(asr.res,paste0(save.name,".",i,".RDS"))
  return(asr.res)
},mc.cores = numCores)
#
