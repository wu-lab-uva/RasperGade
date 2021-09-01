# if you haven't install RasperGade yet
library(devtools)
install_github(repo="wu-lab-uva/RasperGade")

# if you have RasperGade installed
library(RasperGade)

# below is a simple demo of how RasperGade is run
## result may change slightly

# read the demo tree and trait values
tree = read.tree("demo.tree")
trait = read.table("demo.trait.txt",sep = "\t",stringsAsFactors = FALSE)[,2]
names(trait) = read.table("demo.trait.txt",sep = "\t",stringsAsFactors = FALSE)[,1]
# 'fit' a BM model using PIC
bm.rate = var(pic(x = trait[tree$tip.label],phy = tree))
print(bm.rate)
## [1] 87.47837

# fit the pulsed evolution model
## takes several seconds to complete
pe.model = fitPE.phy(phy = tree,x = trait[tree$tip.label],laplace = FALSE)
print(pe.model$params)
##   lambda      size     sigma   epsilon 
## 1.172651 58.925763  6.760760  2.166998 

# conducting leave-one-out cross-validation under different models
## takes about a minute to complete
bm.cv = do.call(rbind,lapply(1:Ntip(tree),function(i){predTraitByPIC(phy = tree,trait = trait[tree$tip.label[-i]],rate = bm.rate)}))
pe.fmr = fullMarginalReconstructionWithPE(phy = tree,x = trait[tree$tip.label],params = pe.model$params,laplace = FALSE)
pe.cv = crossValidationWithPE(FMR = pe.fmr,add.epsilon = FALSE,laplace = FALSE)

# get the residual errors
bm.error = trait[tree$tip.label] - bm.cv$x
pe.error = trait[tree$tip.label] - pe.cv$cv$x

# compare the error and bias
## bias
print(c(BM=mean(bm.error),PE=mean(pe.error)))
###         BM          PE 
### 0.01810704 -0.08808421 

## error
print(c(BM=mean(abs(bm.error)),PE=mean(abs(pe.error))))
###       BM       PE 
### 3.258241 3.156028 

# Which model provide more normal Z-scores? Use D statistic in the KS-test to compare.
bm.D = analyzeResidualErrorByPPplot(trait = trait[tree$tip.label],pred = bm.cv$x,error = bm.cv$var,epsilon = 0,laplace = FALSE)
pe.D = analyzeResidualErrorByPPplot(trait = trait[tree$tip.label],pred = pe.cv$cv$x,error = pe.cv$error,epsilon = pe.fmr$params$epsilon,laplace = FALSE)
print(c(BM=unname(bm.D),PE=unname(pe.D)))
###         BM         PE 
### 0.08680844 0.02119105 
