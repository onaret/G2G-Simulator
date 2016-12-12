library(tidyr)
library(dplyr)
library(parallel)

library(ggplot2)

library(SKAT)
library(globaltest)

source("src/G2.R")
source("src/G2G.R")
source("src/G2G_single.R")
source("src/G2G_full.R")
source("src/GWAS.R")
dir.create("gen-data", showWarnings = F, recursive = F, mode = "0777")
trace <- TRUE

#nb_cpu = 2


