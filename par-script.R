#rm(list = ls())
cat("\014")  
library(rjson)
library(parallel)
library(dplyr)
#library(compiler)
#library(pryr)
#library(microbenchmark)
library(ggplot2)

source("generate-SNP-datas.R")
source("generate-viral-datas.R")
#source("analyse-generated-datas.R")
#source("system-tools.R")

#enableJIT(3)

#######With virus

###Input
threshold = 5e-05
nb_pop = 2
nb_strain = 2
sample_size = 5000
repnb = 1000
trace = TRUE

get_SNP_input <- function(S_Stratified = NA, S_Biased = NA,Y_Stratified = NA, Y_Biased = NA, Associated_Strains = NA, Associated_Populations = NA, size = size, Ordered_Bias = FALSE) {
  lvl = c("full","no","half", "P1", "P2", "A", "B")
  if(is.na(Associated_Strains) & !is.na(Associated_Populations)) Associated_Strains = "full"
  if(is.na(Associated_Populations) & !is.na(Associated_Strains)) Associated_Populations = "full"
  #  if(is.na(Associated_Strains) & !is.na(Associated_Populations)) throw("Define associated_Populations")
  # if(is.na(Associated_Populations) & !is.na(Associated_Strains)) throw("Define Associated_Strains")
  tbl_df(data.frame(`S_Stratified` = factor(rep(S_Stratified, size), levels = lvl), `S_Biased` = factor(S_Biased, levels = lvl), `Y_Stratified` = factor(Y_Stratified, levels = lvl), 
                    `Y_Biased` = factor(Y_Biased, levels = lvl), `Associated_Strains` = Associated_Strains, `Associated_Populations` = Associated_Populations, `Ordered_Bias` = Ordered_Bias,
                    `size` = size, stringsAsFactors = FALSE) )}

fs1 = get_SNP_input(S_Stratified = "full",Y_Biased = "full", size = repnb, Ordered_Bias = TRUE)
fs2_ss1 = get_SNP_input(S_Stratified = "full", Associated_Strains = "full", Associated_Populations = "full", size = repnb)
ss4 = get_SNP_input(S_Stratified = "full", Associated_Strains = "half", Associated_Populations = "half", size = repnb)

cl = makeCluster(30, type = "FORK")

res = parLapply(cl, seq(0.01, 1, length.out = 100), function(param) {
  pval = scenario_viral(sample_size, nb_strain, nb_pop, fs1, fcoeff_pop = param, vir_bias = param, tag = "fs1_strat_impact")$pval
  medianpval = apply(pval, 2, function(p) median(p,na.rm = TRUE))
  rm(pval)
  c(param, medianpval)
})

fs2_beta = parLapply(cl, seq(0.01, 1, length.out = 100), function(param) {
  pval = scenario_viral(sample_size, nb_strain, nb_pop, fs2_ss1, fcoeff_pop = 0.2, beta = param, tag = "beta_0.5")$pval
  medianpval = apply(pval, 2, function(p) median(p,na.rm = TRUE))
  rm(pval)
  c(param, medianpval)
})

param = 3000
pval =  scenario_viral(sample_size= param, nb_strain, nb_pop, ss4, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25")
pval = select(data.frame(pval$pval), ends_with("pval"))
medianpval = apply(pval, 2, function(p) median(p,na.rm = TRUE))
c(param, medianpval)

param = 4000
pval =  scenario_viral(sample_size= param, nb_strain, nb_pop, ss4, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25")
pval = select(data.frame(pval$pval), ends_with("pval"))
medianpval2 = apply(pval, 2, function(p) median(p,na.rm = TRUE))
c(param, medianpval2)


ss4_size = parLapply(cl, seq(3000, 10920, length.out = 100), function(param) {
  pval =  scenario_viral(sample_size= param, nb_strain, nb_pop, ss4, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25")
  pval = select(data.frame(pval$pval), ends_with("pval"))
  medianpval = apply(pval, 2, function(p) median(p,na.rm = TRUE))
  c(param, medianpval)
})



write.table(fs2_beta,file = "result_on_fs2beta")

write.table(ss4_size,file = "ss4_size")
