##############################################################
################## System preparation ########################
##############################################################
rm(list=ls())

### Install CRAN packages, if necessary
if (!require("Seurat", quietly = TRUE))
  install.packages("Seurat", version="5.3.0")
if (!require("MCMCpack", quietly = TRUE))
  install.packages("MCMCpack", version="1.7.1")
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2", version="3.5.2")


## For Windows users, please update RTools to the latest version.



## Install R package spARI from a local file
## (Please make sure to install all required dependencies manually in advance)
install.packages("BayesRare_1.0.tar.gz", repos = NULL, type="source")
## Or install from GitHub 
## (Automatically install dependencies)
devtools::install_github("yinqiaoyan/BayesRare", dependencies = TRUE)




