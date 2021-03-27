# RasperGade: Reconstructing Ancestral State under Pulsed Evolution in R by Gaussian Decomposition
`RasperGade` is an  R package that predicts ancestral and hidden states while accounting for pulsed evolution and time-independent variation.
It has 3 major functions: fitting the pulsed evolution model given extant trait value and a phylogeny, reconstructing ancestral states or predicting hidden states, and evaluating the quality of predicted ancestral and hidden states

Detailed algorithm and analyses are described in the preprint: ...

## System requirements
`RasperGade` are built on `R 3.6.3`. It depends on the following R packages and their dependencies for fitting the pulsed evolution model and ancestral state reconstruction: 

`ape`, `bbmle`, `pracma`, `doParallel`

The following R packages are suggested for plotting diagnostic graphs: 

`ggplot2`, `ggpubr`

## Installation
The required and suggested packages can be installed from CRAN by running the following lines in R:
```
  install.packages("ape")
  install.packages("bbmle")
  install.packages("pracma")
  install.packages("doParallel")
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
`RasperGade` has two basic input data: the phylogeny and the extant trait values.
Phylogenies are expected to be phylo-class objects as defined by the R package `ape`.
Extant trait values should be in a named numeric vector, and the names should correspond to tip labels of the phylogeny. Hidden states to be predicted are marked by NA values in the named vector.
Tips should have unique labels, and if nodes have labels, they should be unique.

## Demo

## Reproduction of results in the preprint
