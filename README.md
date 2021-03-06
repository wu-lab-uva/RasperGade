# RasperGade: Reconstructing Ancestral State under Pulsed Evolution in R by Gaussian Decomposition
`RasperGade` is an  R package that predicts ancestral and hidden states while accounting for pulsed evolution and time-independent variation.
It has 3 major functions:

1. fitting the pulsed evolution model given extant trait value and a phylogeny

2. reconstructing ancestral states or predicting hidden states

3. evaluating the quality of predicted ancestral and hidden states

Detailed algorithm and analyses are described in the preprint: `Modeling pulsed evolution and time-independent variation improves the confidence level of ancestral and hidden state predictions in continuous traits`

doi: https://doi.org/10.1101/2021.03.29.437517

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
A small example is provided under the folder `inst/extdata/Demo`
Once `RasperGade` is installed, the location of the example can also be accessed by the following command in R
```
demo.tree.file = system.file("extdata/Demo/demo.tree",package = "RasperGade",mustWork = TRUE)
demo.trait.file = system.file("extdata/Demo/demo.trait.txt",package = "RasperGade",mustWork = TRUE)
```