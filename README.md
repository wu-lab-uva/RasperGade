# RasperGade: Reconstructing Ancestral State under Pulsed Evolution in R by Gaussian Decomposition
`RasperGade` is an  R package that predicts ancestral and hidden states while accounting for pulsed evolution and time-independent variation.
It has 3 major functions:

1. fitting the pulsed evolution model given extant trait value and a phylogeny

2. reconstructing ancestral states or predicting hidden states

3. evaluating the quality of predicted ancestral and hidden states

Detailed algorithm and analyses are described in

Yingnan Gao, Martin Wu, Modeling Pulsed Evolution and Time-Independent Variation Improves the Confidence Level of Ancestral and Hidden State Predictions, Systematic Biology, 2022;, syac016, https://doi.org/10.1093/sysbio/syac016

## System requirements
`RasperGade` are built on `R 3.6.3`. It requires the following R packages and their dependencies for fitting the pulsed evolution model and ancestral state reconstruction: 

`ape`, `bbmle`, `castor`,`digest`, `parallel`,`pracma`

The following R packages are suggested for plotting diagnostic graphs: 

`ggplot2`, `ggpubr`

## Installation
The required packages can be installed from CRAN by running the following lines in R:
```
  install.packages("ape")
  install.packages("bbmle")
  install.packages("castor")
  install.packages("digest")
  install.packages("parallel")
  install.packages("pracma")
```
The suggested packages can also be installed from CRAN:
```
  install.packages("ggplot2")
  install.packages("ggpubr")
```
`RasperGade` can be installed from Github using the `devtools` package:
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github(repo = "wu-lab-uva/RasperGade")
```
## Data format
`RasperGade` takes two basic input data: the phylogeny and the extant trait values, to perform a typical ancestral state reconstruction or hidden state prediction.

Phylogenies are expected to be phylo-class objects as defined by the R package `ape`.

Extant trait values are expected to be a named numeric vector, and the names should correspond to tip labels of the phylogeny. Hidden states to be predicted are either absent or marked by NA values in the named vector.

Tips should have unique labels, and if nodes have labels, they should be unique.

## Demo
A small demo is provided under the folder `inst/extdata/Demo`
Once `RasperGade` is installed, the location of the demo files can also be accessed by the following command in R
```
demo.tree.file = system.file("extdata/Demo/demo.tree",package = "RasperGade",mustWork = TRUE)
demo.trait.file = system.file("extdata/Demo/demo.trait.txt",package = "RasperGade",mustWork = TRUE)
```
To run a typical leave-one-out cross-validation analysis
```
library(RasperGade)

# read the demo tree and trait values
tree = read.tree(demo.tree.file)
trait = read.table(demo.trait.file,sep = "\t",stringsAsFactors = FALSE)[,2]
names(trait) = read.table(demo.trait.file,sep = "\t",stringsAsFactors = FALSE)[,1]

# fit the pulsed evolution model
pe.model = fitPE.phy(phy = tree,x = trait[tree$tip.label],laplace = FALSE)
print(pe.model$params)
##   lambda      size     sigma   epsilon 
## 1.172651 58.925763  6.760760  2.166998 

# conducting leave-one-out cross-validation under different models
bm.cv = LOO_CV_with_BM(tree,trait)
pe.fmr = fullMarginalReconstructionWithPE(phy = tree,x = trait[tree$tip.label],params = pe.model$params,laplace = FALSE)
pe.cv = LOO_CV_with_PE(FMR = pe.fmr)
```
To evaluate the performance of CV
```
# get the residual errors
bm.error = trait[tree$tip.label] - bm.cv$summary$x
pe.error = trait[tree$tip.label] - pe.cv$summary$x

# compare the error and bias
## bias
print(c(BM=mean(bm.error),PE=mean(pe.error)))
###         BM          PE 
### 0.01810704 -0.08808421 

## error
print(c(BM=mean(abs(bm.error)),PE=mean(abs(pe.error))))
###       BM       PE 
### 3.258241 3.156028 

# get the residual Z-scores
bm.pZ = pseudo.Z.score(trait=trait[tree$tip.label],pred=bm.cv$summary$x,error=bm.cv$error)
bm.pZ$se = sqrt(bm.cv$summary$var)
pe.pZ = pseudo.Z.score(trait=trait[tree$tip.label],pred=pe.cv$summary$x,error=pe.cv$error)
pe.pZ$se = sqrt(pe.cv$summary$var)

# calculate diagnostic statistics
diag.stats = lapply(list(BM=bm.pZ,PE=pe.pZ),function(x){
  this.test = ks.test(x$q,pnorm)
  this.D = this.test$statistic
  this.H = calculateHeteroscedasticity(x = x$q,y = x$se,bin = 20)[1]
  this.RSS = calculateBinomRSS(p = rep(0.95,dim(x)[1]),
                               obs = sapply(x$p,function(pp){
                                 min(pp,1-pp)>=0.025}),group = log(x$se),bin = 10,equal.bin = FALSE)
  data.frame(D = unname(this.D),H = unname(this.H),RSS = unname(this.RSS))
})
print(diag.stats)
```
