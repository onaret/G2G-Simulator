library(tidyr)
library(dplyr)
library(parallel)

library(ggplot2)

library(SKAT)
library(globaltest)
setwd(dirname(parent.frame(2)$ofile))
source("G2.R")
source("G2G.R")
source("G2G_single.R")
source("G2G_full.R")
source("GWAS.R")

trace <- TRUE

#nb_cpu = 2
dir.create("../gen-data", showWarnings = F, recursive = F, mode = "0777")