library(tidyr)
library(dplyr)
library(parallel)
library(ggplot2)

###If you want to use SKAT, globaltest or G2 aside logistic regression.
##It was used but is not fully supported at the moment
##Uncomment following...

#library(SKAT)
#library(globaltest)
#source("src/G2.R")

source("src/model/G2G.R")
source("src/model/GWAS.R")
source("src/controller/G2G.R")
source("src/controller/GWAS.R")
source("src/view/G2G.R")
source("src/view/GWAS.R")

trace <- TRUE
nb_cpu = 20

dir.create("gen-data", showWarnings = F, recursive = F, mode = "0777")
