# load required packages and functions
source("../R/ape_addition.R")
source("../R/PEpois.R")
source("../R/RasperGade_reconstruction.R")
source("../R/fitPE_functions.R")

# read in key arguments for the simulation
arg=commandArgs(TRUE)
# interpret the arguments
# the job name
job.name = arg[1]
if(is.na(job.name)) stop("Job name is required and cannot be omitted!")
# read in phylogeny
data.file = sprintf("%s.tre",job.name)
tree = read.tree(data.file)
# how many cores to use in the fitting
numCores=as.numeric(arg[2])
laplace=FALSE
#
epsilons=c(0,1,10)
lambdas = c(5)
contributions = rev(c(0))
replicates=1:100
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
  saveRDS(asr.res,paste0(job.name,".simulation.epsilon.only.",i,".RDS"))
  t0 = Sys.time()
  FMR = fullMarginalReconstructionWithPE(phy = tree,x = trait[tree$tip.label],params = model.params,
                                         laplace = laplace,approximate = 1,numCores = 1,asymptotic=5)
  CV = crossValidationWithPE(FMR = FMR,add.epsilon = FALSE,numCores = 1,laplace = laplace)
  ASR = globalReconstructionWithPE(FMR = FMR,add.epsilon = FALSE,numCores = 1,laplace = laplace)
  t1 = Sys.time()
  asr.res[[2]] = list(FMR=FMR,trait=trait,CV=CV,ASR=ASR,t = c(t0,t1))
  saveRDS(asr.res,paste0(job.name,".simulation.epsilon.only.",i,".RDS"))
  return(asr.res)
},mc.cores = numCores)
#


